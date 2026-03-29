#ifndef ROBSCALE_ROBUST_CORE_H
#define ROBSCALE_ROBUST_CORE_H

#include <cmath>
#include <cstring>
#include "selection.h"
#include "sort_net.h"
#include "robscale_config.h"
#include "qnsn_runtime_config.h"
#include "simd_median.h"

// --- Platform-specific vectorized tanh ---
#if defined(__APPLE__) && defined(__MACH__)
  #define COMPLEX vDSP_COMPLEX
  #include <Accelerate/Accelerate.h>
  #undef COMPLEX
  #define ROBSCALE_HAS_ACCELERATE 1
#endif

// SLEEF header — only when SLEEF is actually installed
#if defined(ROBSCALE_HAS_SLEEF) && !defined(ROBSCALE_HAS_ACCELERATE)
  #include <sleef.h>
  #if defined(__aarch64__) || defined(_M_ARM64)
    #include <arm_neon.h>
  #endif
#endif

// Vectorized tanh declarations.
// Backend hierarchy (highest priority first):
//   1. macOS Accelerate vvtanh (ROBSCALE_HAS_ACCELERATE)
//   2. glibc libmvec _ZGVdN4v_tanh — AVX2 4-wide   (ROBSCALE_HAS_GLIBC_MVEC)
//   3. SLEEF Sleef_tanhd4_u10avx2  — AVX2 4-wide   (ROBSCALE_HAS_SLEEF)
//   4. #pragma omp simd / scalar std::tanh
//
// ROBSCALE_HAS_AVX2_TANH: defined when AVX2 tanh backend available (2 or 3).
#if defined(ROBSCALE_HAS_AVX2_TANH) && !defined(ROBSCALE_HAS_ACCELERATE)
  #include <immintrin.h>
  // AVX2 4-wide backend
  // Declarations carry target("avx2,fma") so clang 21+ accepts __m256d
  // parameter/return types without global -mavx2 (ABI compatibility).
  #if defined(ROBSCALE_HAS_GLIBC_MVEC)
    extern "C" ROBSCALE_TARGET_AVX2 __m256d _ZGVdN4v_tanh(__m256d);
    #define ROBSCALE_TANH4_AVX2 _ZGVdN4v_tanh
  #elif defined(ROBSCALE_HAS_SLEEF)
    // sleef.h guards Sleef_tanhd4_u10avx2 behind #ifdef __AVX2__, which is
    // absent without global -mavx2.  The symbol exists in libsleef regardless;
    // declare it here so the target-attributed wrapper can call it.
    extern "C" ROBSCALE_TARGET_AVX2 __m256d Sleef_tanhd4_u10avx2(__m256d);
    #define ROBSCALE_TANH4_AVX2 Sleef_tanhd4_u10avx2
  #endif
#endif

namespace robscale {

// Scan for first non-finite element. Returns its index, or -1 if all finite.
// No allocation; short-circuits on first bad value.
// Callers that include Rcpp.h wrap this in a static validate_finite() helper
// that issues Rcpp::stop with the appropriate message.
ROBSCALE_HIDDEN ROBSCALE_INLINE
int find_first_nonfinite(const double* ROBSCALE_RESTRICT xp, int n) {
  for (int i = 0; i < n; ++i)
    if (ROBSCALE_UNLIKELY(!std::isfinite(xp[i]))) return i;
  return -1;
}

// --- Asymptotic Relative Efficiency (ARE) constants ---

constexpr double ARE_ROBSCALE   = 0.55;
constexpr double ARE_QN         = 0.82;
constexpr double ARE_SN         = 0.58;
constexpr double ARE_GMD        = 0.98;
constexpr double ARE_MAD        = 0.368;
constexpr double ARE_IQR        = 0.37;
constexpr double ARE_ADM        = 0.88;
constexpr double ARE_SD_C4      = 1.00;

// --- Low-level Kernels ---

// ADM core: constant * mean(|x - center|).
// Scalar version: safe on all x86 (no target attribute).
ROBSCALE_HIDDEN ROBSCALE_INLINE double adm_core_scalar(const double* ROBSCALE_RESTRICT x, int n,
                                                       double center, double constant) {
  double s0 = 0.0, s1 = 0.0, s2 = 0.0, s3 = 0.0;
  int i = 0;
  for (; i + 3 < n; i += 4) {
    s0 += std::abs(x[i]   - center);
    s1 += std::abs(x[i+1] - center);
    s2 += std::abs(x[i+2] - center);
    s3 += std::abs(x[i+3] - center);
  }
  for (; i < n; ++i) s0 += std::abs(x[i] - center);
  return constant * (s0 + s1 + s2 + s3) * (1.0 / n);
}

// ADM core (AVX2): constant * mean(|x - center|)
// ROBSCALE_TARGET_AVX2: enables 256-bit AVX2 auto-vectorisation for this
// function without requiring a global -mavx2 flag (safe on CRAN).
ROBSCALE_HIDDEN ROBSCALE_TARGET_AVX2
ROBSCALE_INLINE double adm_core_avx2(const double* ROBSCALE_RESTRICT x, int n,
                                     double center, double constant) {
  double s0 = 0.0, s1 = 0.0, s2 = 0.0, s3 = 0.0;
  int i = 0;
  for (; i + 3 < n; i += 4) {
    s0 += std::abs(x[i]   - center);
    s1 += std::abs(x[i+1] - center);
    s2 += std::abs(x[i+2] - center);
    s3 += std::abs(x[i+3] - center);
  }
  for (; i < n; ++i) s0 += std::abs(x[i] - center);
  return constant * (s0 + s1 + s2 + s3) * (1.0 / n);
}

// ADM core dispatch: AVX2 if available, scalar otherwise.
ROBSCALE_HIDDEN ROBSCALE_INLINE double adm_core(const double* ROBSCALE_RESTRICT x, int n,
                                                 double center, double constant,
                                                 bool use_avx2 = false) {
#if defined(ROBSCALE_HAS_AVX2_TANH) && !defined(ROBSCALE_HAS_ACCELERATE)
  if (use_avx2) return adm_core_avx2(x, n, center, constant);
#endif
  (void)use_avx2;
  return adm_core_scalar(x, n, center, constant);
}

// ADM core for pre-sorted input: constant * mean(|x - center|).
// Scalar version: safe on all x86.
ROBSCALE_HIDDEN ROBSCALE_INLINE double adm_core_sorted_scalar(const double* ROBSCALE_RESTRICT x, int n,
                                                              double center, double constant) {
  if (ROBSCALE_UNLIKELY(n <= 1)) return 0.0;
  int k = n / 2;
  double sum_abs = 0.0;
  for (int i = 0; i < k; ++i)             sum_abs += center - x[i];
  for (int i = k + (n & 1); i < n; ++i)   sum_abs += x[i] - center;
  return constant * sum_abs * (1.0 / n);
}

// ADM core for pre-sorted input (AVX2).
ROBSCALE_HIDDEN ROBSCALE_TARGET_AVX2
ROBSCALE_INLINE double adm_core_sorted_avx2(const double* ROBSCALE_RESTRICT x, int n,
                                            double center, double constant) {
  if (ROBSCALE_UNLIKELY(n <= 1)) return 0.0;
  int k = n / 2;
  double sum_abs = 0.0;
  for (int i = 0; i < k; ++i)             sum_abs += center - x[i];
  for (int i = k + (n & 1); i < n; ++i)   sum_abs += x[i] - center;
  return constant * sum_abs * (1.0 / n);
}

// ADM core for pre-sorted input dispatch.
ROBSCALE_HIDDEN ROBSCALE_INLINE double adm_core_sorted(const double* ROBSCALE_RESTRICT x, int n,
                                                       double center, double constant,
                                                       bool use_avx2 = false) {
#if defined(ROBSCALE_HAS_AVX2_TANH) && !defined(ROBSCALE_HAS_ACCELERATE)
  if (use_avx2) return adm_core_sorted_avx2(x, n, center, constant);
#endif
  (void)use_avx2;
  return adm_core_sorted_scalar(x, n, center, constant);
}

// Absolute-deviation kernel (in-place): w[i] = |w[i] - center|.
// Single-pointer; RESTRICT valid since no other pointer aliases w in the caller.
// Used by mad_impl_auto paths and mad_from_data.
ROBSCALE_HIDDEN ROBSCALE_INLINE void bulk_abs_diff_inplace(
    double* ROBSCALE_RESTRICT w, int n, double center) {
#if defined(_OPENMP) || defined(ROBSCALE_HAS_OMP_SIMD)
  #pragma omp simd
#endif
  for (int i = 0; i < n; ++i)
    w[i] = std::abs(w[i] - center);
}

// Absolute-deviation kernel (out-of-place): dst[i] = |src[i] - center|.
// dst and src guaranteed non-aliasing by caller; RESTRICT enables SIMD.
// Used by mad_impl_center paths and the ensemble MAD block.
ROBSCALE_HIDDEN ROBSCALE_INLINE void bulk_abs_diff(
    double* ROBSCALE_RESTRICT dst,
    const double* ROBSCALE_RESTRICT src,
    int n, double center) {
#if defined(_OPENMP) || defined(ROBSCALE_HAS_OMP_SIMD)
  #pragma omp simd
#endif
  for (int i = 0; i < n; ++i)
    dst[i] = std::abs(src[i] - center);
}

// Median of pre-sorted array. Precondition: n >= 1 (n=0 guard is defensive).
ROBSCALE_HIDDEN ROBSCALE_INLINE double median_sorted(const double* x, size_t n) {
  if (ROBSCALE_UNLIKELY(n == 0)) return 0.0;
  if (n & 1) return x[n / 2];
  return x[(n / 2) - 1] + (x[n / 2] - x[(n / 2) - 1]) * 0.5;  // overflow-safe
}

// Runtime CPUID check for SIMD median dispatch (cached static bool).
#ifdef ROBSCALE_HAS_AVX2_DISPATCH
namespace {
inline bool cpu_has_avx2() {
  static const bool val = (robscale::qnsn::RuntimeConfig::get().hw.simd_level >=
                           robscale::qnsn::SIMDLevel::AVX2);
  return val;
}
} // anon namespace
#endif

// Selection based median
ROBSCALE_HIDDEN ROBSCALE_INLINE double median_select(double* x, size_t n) {
  if (ROBSCALE_UNLIKELY(n == 0)) return 0.0;
#ifdef ROBSCALE_HAS_AVX2_DISPATCH
  if (cpu_has_avx2()) {
    switch (n) {
      case 8:  return robscale::simd::simd_median_sel_8(x);
      case 16: return robscale::simd::simd_median_sel_16(x);
      case 32: return robscale::simd::simd_median_sel_32(x);
    }
  }
#endif
  if (n <= ROBSCALE_MEDIAN_NET_THRESHOLD) {
    return robscale::median_net(x, n);
  }
  size_t h = (n - 1) / 2;
  robscale::floyd_rivest_select(x, x + h, x + n);
  if (n & 1) return x[h];

  double v1 = x[h];
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
ROBSCALE_HIDDEN inline double mad_select(const double* x, int n, double med, double* dev) {
  for (int i = 0; i < n; ++i) dev[i] = std::abs(x[i] - med);
  return robscale::MAD_CONSISTENCY * median_select(dev, static_cast<size_t>(n));
}

// GMD weighted sum: computes scale * sum_i{ w_i * x[i] } where
// x is pre-sorted ascending and w_i = 2*i - (n-1) (0-indexed).
// Dispatches to AVX2 FMA 4-wide kernel (ROBSCALE_HAS_AVX2_TANH guard reuses
// the same AVX2 dispatch infrastructure; target attribute emits avx2/fma
// instructions without global -mavx2).
// Scalar tail handles remainder; same scalar path used when AVX2 absent.
// Threshold n<8: fewer than two full AVX2 vectors, scalar is faster.
#if defined(ROBSCALE_HAS_AVX2_TANH) && !defined(ROBSCALE_HAS_ACCELERATE)
ROBSCALE_HIDDEN ROBSCALE_TARGET_AVX2
ROBSCALE_INLINE double gmd_weighted_sum_avx2(
    const double* ROBSCALE_RESTRICT x, int n, double scale) {
  // Weight w_i = 2*i - (n-1).  For 4-wide lane at offset i:
  //   {w, w+2, w+4, w+6}  where w = 2*i - (n-1).
  // After 4 elements the base increments by 8.
  const double w0  = -static_cast<double>(n - 1);
  const double inc = 8.0;
  __m256d wv  = _mm256_set_pd(w0 + 6.0, w0 + 4.0, w0 + 2.0, w0);
  const __m256d d8 = _mm256_set1_pd(inc);
  __m256d acc = _mm256_setzero_pd();
  int i = 0;
  for (; i + 4 <= n; i += 4) {
    __m256d xv = _mm256_loadu_pd(x + i);
    acc = _mm256_fmadd_pd(wv, xv, acc);
    wv  = _mm256_add_pd(wv, d8);
  }
  // Horizontal sum of acc
  __m128d lo   = _mm256_castpd256_pd128(acc);
  __m128d hi   = _mm256_extractf128_pd(acc, 1);
  __m128d s128 = _mm_add_pd(lo, hi);
  s128         = _mm_hadd_pd(s128, s128);
  double sum   = _mm_cvtsd_f64(s128);
  // Scalar tail: weights continue from where vector left off
  double wi = w0 + 2.0 * i;
  for (; i < n; ++i, wi += 2.0) sum += wi * x[i];
  return scale * sum;
}
#endif

// use_avx2: pre-hoisted AVX2 flag (from RuntimeConfig::get().hw.simd_level);
// avoids a TLS read on every call.  Callers that run inside hot loops (e.g.
// ensemble_one_replicate, cpp_single_estimator_ci_bounds) hoist once before
// the loop and pass the flag here.  Default false keeps non-hot callers working.
ROBSCALE_HIDDEN ROBSCALE_INLINE double gmd_weighted_sum(
    const double* ROBSCALE_RESTRICT x, int n, double scale,
    bool use_avx2 = false) {
#if defined(ROBSCALE_HAS_AVX2_TANH) && !defined(ROBSCALE_HAS_ACCELERATE)
  if (use_avx2 && n >= 8)
    return gmd_weighted_sum_avx2(x, n, scale);
#endif
  double sum = 0.0;
  const double w0 = -static_cast<double>(n - 1);
  for (int i = 0; i < n; ++i)
    sum += (w0 + 2.0 * i) * x[i];
  return scale * sum;
}

} // namespace robscale

#endif // ROBSCALE_ROBUST_CORE_H
