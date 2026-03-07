#ifndef ROBSCALE_ROBUST_CORE_H
#define ROBSCALE_ROBUST_CORE_H

#include <cmath>
#include <cstring>
#include <memory>
#include "selection.h"
#include "sort_net.h"

// --- Platform-specific vectorized tanh ---
// On macOS, include only vForce (not full Accelerate) to avoid COMPLEX
// symbol collision between vecLib/vDSP.h and R's Rinternals.h
#if defined(__APPLE__) && defined(__MACH__)
  #include <vecLib/vForce.h>
  #define ROBSCALE_HAS_ACCELERATE 1
#endif

// SLEEF detection (set via PKG_CXXFLAGS from configure)
#if defined(ROBSCALE_HAS_SLEEF)
  #include <sleef.h>
  #if defined(__x86_64__) || defined(_M_X64)
    #include <immintrin.h>
  #elif defined(__aarch64__) || defined(_M_ARM64)
    #include <arm_neon.h>
  #endif
#endif

// --- Macros ---
#if defined(__GNUC__) || defined(__clang__)
  #define ROBSCALE_LIKELY(x)   __builtin_expect(!!(x), 1)
  #define ROBSCALE_UNLIKELY(x) __builtin_expect(!!(x), 0)
#else
  #define ROBSCALE_LIKELY(x)   (x)
  #define ROBSCALE_UNLIKELY(x) (x)
#endif

#ifdef _OPENMP
  #pragma omp declare simd notinbranch
#endif
inline double vec_tanh(double x) { return std::tanh(x); }

// Bulk tanh: vectorized via Accelerate (macOS), SLEEF (Linux), or OpenMP SIMD
inline void bulk_tanh(double* inout, int n) {
#if defined(ROBSCALE_HAS_ACCELERATE)
  vvtanh(inout, inout, &n);
#elif defined(ROBSCALE_HAS_SLEEF)
  int i = 0;
  #if defined(__AVX2__)
    for (; i + 4 <= n; i += 4) {
      __m256d v = _mm256_loadu_pd(inout + i);
      v = Sleef_tanhd4_u10avx2(v);
      _mm256_storeu_pd(inout + i, v);
    }
  #elif defined(__AVX512F__)
    for (; i + 8 <= n; i += 8) {
      __m512d v = _mm512_loadu_pd(inout + i);
      v = Sleef_tanhd8_u10avx512f(v);
      _mm512_storeu_pd(inout + i, v);
    }
  #elif defined(__aarch64__)
    for (; i + 2 <= n; i += 2) {
      float64x2_t v = vld1q_f64(inout + i);
      v = Sleef_tanhd2_u10(v);
      vst1q_f64(inout + i, v);
    }
  #endif
  // Scalar fallback for remainder
  for (; i < n; i++) {
    inout[i] = std::tanh(inout[i]);
  }
#else
  #ifdef _OPENMP
    #pragma omp simd
  #endif
  for (int i = 0; i < n; ++i) {
    inout[i] = vec_tanh(inout[i]);
  }
#endif
}

// Key constants from Rousseeuw & Verboven (2002)
constexpr double RHO_SCALE_CONST      = 0.37394112142347236;
constexpr double INV_RHO_SCALE_CONST  = 1.0 / 0.37394112142347236;
constexpr double MAD_CONSISTENCY      = 1.4826;
constexpr double ADM_CONSISTENCY      = 1.2533141373155001;  // sqrt(pi/2)

// Stack buffer size (2048 doubles = 16 KB per segment; 3x arena = 48 KB)
constexpr int STACK_BUF_SIZE = 2048;

// Median of a pre-sorted array (caller must ensure n >= 1)
inline double median_sorted(const double* x, int n) {
  if (n & 1) {
    return x[n >> 1];
  } else {
    return (x[(n >> 1) - 1] + x[n >> 1]) * 0.5;
  }
}

// ADM core: constant * mean(|x - center|)
inline double adm_core(const double* x, int n, double center, double constant) {
  double sum = 0.0;
  for (int i = 0; i < n; ++i) {
    sum += std::abs(x[i] - center);
  }
  return constant * sum / n;
}

// --- Selection ---
inline double median_select(double* x, size_t n) {
  if (ROBSCALE_UNLIKELY(n == 0)) return 0.0;
  if (n <= 8) {
    robscale::small_sort(x, n);
    return median_sorted(x, n);
  }
  size_t k_idx = (n - 1) / 2;
  double* k = x + k_idx;
  robscale::floyd_rivest_select(x, k, x + n);
  if (n & 1) return *k;
  
  // For even n, we need the next element
  double v1 = *k;
  double v2 = x[k_idx + 1];
  for (size_t i = k_idx + 2; i < n; ++i) {
    if (x[i] < v2) v2 = x[i];
  }
  return (v1 + v2) * 0.5;
}

// MAD via selection: |x_i - med| into dev buffer, then median_select
inline double mad_select(const double* x, int n, double med, double* dev) {
  for (int i = 0; i < n; ++i) {
    dev[i] = std::abs(x[i] - med);
  }
  return MAD_CONSISTENCY * median_select(dev, static_cast<size_t>(n));
}

#endif // ROBSCALE_ROBUST_CORE_H
