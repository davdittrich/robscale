#ifndef ROBSCALE_ROBUST_CORE_H
#define ROBSCALE_ROBUST_CORE_H

#include <cmath>
#include <cstring>
#include "selection.h"
#include "sort_net.h"
#include "robscale_config.h"

// --- Platform-specific vectorized tanh ---
#if defined(__APPLE__) && defined(__MACH__)
  #define COMPLEX vDSP_COMPLEX
  #include <Accelerate/Accelerate.h>
  #undef COMPLEX
  #define ROBSCALE_HAS_ACCELERATE 1
#endif

#if defined(ROBSCALE_HAS_SLEEF) && !defined(ROBSCALE_HAS_ACCELERATE)
  #include <sleef.h>
  #if defined(__x86_64__) || defined(_M_X64)
    #include <immintrin.h>
  #elif defined(__aarch64__) || defined(_M_ARM64)
    #include <arm_neon.h>
  #endif
#endif

namespace robscale {

// --- Low-level Kernels ---

// Bulk tanh: vectorized via Accelerate (macOS), SLEEF (Linux), or OpenMP SIMD
inline void bulk_tanh(double* inout, int n) {
  if (n <= 64) {
    for (int i = 0; i < n; ++i) inout[i] = std::tanh(inout[i]);
    return;
  }
#if defined(ROBSCALE_HAS_ACCELERATE)
  vvtanh(inout, inout, &n);
#elif defined(ROBSCALE_HAS_SLEEF)
  int i = 0;
  #if defined(__AVX512F__)
    for (; i + 8 <= n; i += 8) {
      __m512d v = _mm512_loadu_pd(inout+i);
      v = Sleef_tanhd8_u10avx512f(v);
      _mm512_storeu_pd(inout+i, v);
    }
  #elif defined(__AVX2__)
    for (; i + 4 <= n; i += 4) {
      __m256d v = _mm256_loadu_pd(inout+i);
      v = Sleef_tanhd4_u10avx2(v);
      _mm256_storeu_pd(inout+i, v);
    }
  #endif
  for (; i < n; i++) inout[i] = std::tanh(inout[i]);
#else
  #if defined(_OPENMP) || defined(ROBSCALE_HAS_OMP_SIMD)
    #pragma omp simd
  #endif
  for (int i = 0; i < n; ++i) inout[i] = std::tanh(inout[i]);
#endif
}

// ADM core: constant * mean(|x - center|)
ROBSCALE_INLINE double adm_core(const double* x, int n, double center, double constant) {
  double sum = 0.0;
#if defined(_OPENMP) || defined(ROBSCALE_HAS_OMP_SIMD)
  #pragma omp simd reduction(+:sum)
#endif
  for (int i = 0; i < n; ++i) {
    sum += std::abs(x[i] - center);
  }
  return constant * sum / n;
}

// robLoc iteration kernel
inline double rob_loc_kernel(const double* x, int n, double t, double* tmp) {
  for (int k = 0; k < 100; ++k) {
    for (int i = 0; i < n; ++i) tmp[i] = (x[i] - t) * 0.5; // Scale fixed at 1 for inner
    bulk_tanh(tmp, n);
    double sum_psi = 0.0, sum_dpsi = 0.0;
    for (int i = 0; i < n; ++i) {
      double p = tmp[i];
      sum_psi += p;
      sum_dpsi += 1.0 - p*p;
    }
    double v = 2.0 * sum_psi / sum_dpsi;
    t += v;
    if (std::abs(v) <= 1e-10) break;
  }
  return t;
}

// Median of pre-sorted array
ROBSCALE_INLINE double median_sorted(const double* x, size_t n) {
  if (n & 1) return x[n / 2];
  return (x[(n / 2) - 1] + x[n / 2]) * 0.5;
}

// Selection based median
ROBSCALE_INLINE double median_select(double* x, size_t n) {
  if (ROBSCALE_UNLIKELY(n == 0)) return 0.0;
  if (n <= 16) {
    robscale::small_sort(x, n);
    return median_sorted(x, n);
  }
  size_t h = (n - 1) / 2;
  robscale::floyd_rivest_select(x, x + h, x + n);
  if (n & 1) return x[h];
  
  double v1 = x[h];
  double v2 = x[h + 1];
  for (size_t i = h + 2; i < n; ++i) {
    if (x[i] < v2) v2 = x[i];
  }
  return (v1 + v2) * 0.5;
}

// Low-median selection
inline double lowmedian_select(double* x, size_t n) {
  if (ROBSCALE_UNLIKELY(n == 0)) return 0.0;
  size_t h = (n - 1) / 2;
  robscale::floyd_rivest_select(x, x + h, x + n);
  return x[h];
}

// MAD via selection
inline double mad_select(const double* x, int n, double med, double* dev) {
  for (int i = 0; i < n; ++i) dev[i] = std::abs(x[i] - med);
  return robscale::MAD_CONSISTENCY * median_select(dev, static_cast<size_t>(n));
}

} // namespace robscale

#endif // ROBSCALE_ROBUST_CORE_H
