#ifndef ROBSCALE_ROBUST_CORE_H
#define ROBSCALE_ROBUST_CORE_H

#include <cmath>
#include <cstring>
#include "selection.h"
#include "sort_net.h"
#include "robscale_config.h"
#include "qnsn_runtime_config.h"

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
    #ifdef ROBSCALE_HAS_AVX2_DISPATCH
    // sleef.h guards Sleef_tanhd4_u10avx2 behind #ifdef __AVX2__, which is
    // absent without global -mavx2.  The symbol exists in libsleef regardless;
    // declare it here so the target-attributed wrapper can call it.
    extern "C" __m256d Sleef_tanhd4_u10avx2(__m256d);
    #endif
  #elif defined(__aarch64__) || defined(_M_ARM64)
    #include <arm_neon.h>
  #endif
#endif

namespace robscale {

// --- Low-level Kernels ---

#if defined(ROBSCALE_HAS_SLEEF) && !defined(ROBSCALE_HAS_ACCELERATE)
  #ifdef ROBSCALE_HAS_AVX2_DISPATCH
  // AVX2 SLEEF tanh: processes 4 doubles per iteration.
  // Target attribute enables AVX2 codegen without global -mavx2.
  ROBSCALE_TARGET_AVX2
  inline void bulk_tanh_sleef_avx2(double* inout, int n) {
    int i = 0;
    for (; i + 4 <= n; i += 4) {
      __m256d v = _mm256_loadu_pd(inout + i);
      v = Sleef_tanhd4_u10avx2(v);
      _mm256_storeu_pd(inout + i, v);
    }
    for (; i < n; i++) inout[i] = std::tanh(inout[i]);
  }
  #endif
#endif

// Bulk tanh: vectorized via Accelerate (macOS), SLEEF (Linux), or OpenMP SIMD
inline void bulk_tanh(double* inout, int n) {
  if (n <= 64) {
    for (int i = 0; i < n; ++i) inout[i] = std::tanh(inout[i]);
    return;
  }
#if defined(ROBSCALE_HAS_ACCELERATE)
  vvtanh(inout, inout, &n);
#elif defined(ROBSCALE_HAS_SLEEF)
  #ifdef ROBSCALE_HAS_AVX2_DISPATCH
  if (robscale::qnsn::RuntimeConfig::get().hw.simd_level >=
      robscale::qnsn::SIMDLevel::AVX2) {
    bulk_tanh_sleef_avx2(inout, n);
    return;
  }
  #endif
  // Scalar SLEEF fallback (no AVX2 at runtime, or SLEEF on non-x86)
  for (int i = 0; i < n; i++) inout[i] = std::tanh(inout[i]);
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
  if (n <= ROBSCALE_SORT_MEDIAN_THRESHOLD) {
    robscale::small_sort(x, n);
    return median_sorted(x, n);
  }
  size_t h = (n - 1) / 2;
  robscale::floyd_rivest_select(x, x + h, x + n);
  if (n & 1) return x[h];

  double v1 = x[h];
  // After floyd_rivest_select, x[h+1..n-1] are all >= x[h].
  // min_element vectorises to MINPD and has no loop-carried branch dependency.
  double v2 = *std::min_element(x + h + 1, x + n);
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
