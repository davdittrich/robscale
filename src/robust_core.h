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
    // glibc libmvec _ZGVdN4v_tanh is 25-50% faster than SLEEF on glibc >= 2.35
    // systems (Zen/Skylake). Preferred when ROBSCALE_HAS_GLIBC_MVEC is set.
    // SLEEF remains the fallback on older glibc and non-glibc platforms.
    #if defined(ROBSCALE_HAS_GLIBC_MVEC)
    extern "C" __m256d _ZGVdN4v_tanh(__m256d);
    #define ROBSCALE_TANH4_AVX2 _ZGVdN4v_tanh
    #else
    #define ROBSCALE_TANH4_AVX2 Sleef_tanhd4_u10avx2
    #endif
    #endif
  #elif defined(__aarch64__) || defined(_M_ARM64)
    #include <arm_neon.h>
  #endif
#endif

namespace robscale {

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

#if defined(ROBSCALE_HAS_SLEEF) && !defined(ROBSCALE_HAS_ACCELERATE)
  #ifdef ROBSCALE_HAS_AVX2_DISPATCH
  // AVX2 tanh: processes 4 doubles per iteration via ROBSCALE_TANH4_AVX2.
  // Resolves to glibc libmvec _ZGVdN4v_tanh (preferred, 25-50% faster) when
  // ROBSCALE_HAS_GLIBC_MVEC is defined; falls back to Sleef_tanhd4_u10avx2.
  // Target attribute enables AVX2 codegen without global -mavx2.
  ROBSCALE_TARGET_AVX2
  inline void bulk_tanh_sleef_avx2(double* inout, int n) {
    int i = 0;
    for (; i + 4 <= n; i += 4) {
      __m256d v = _mm256_loadu_pd(inout + i);
      v = ROBSCALE_TANH4_AVX2(v);
      _mm256_storeu_pd(inout + i, v);
    }
    for (; i < n; i++) inout[i] = std::tanh(inout[i]);
  }
  #endif
#endif

// Bulk tanh: vectorized via Accelerate (macOS), SLEEF (Linux), or OpenMP SIMD
// Scalar fallback only for n<8: at n<8 SIMD setup overhead exceeds the gain
// (fewer than 2 full AVX2 vectors). For n>=8, use SIMD path.
ROBSCALE_HIDDEN inline void bulk_tanh(double* inout, int n) {
  if (n < 8) {
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

// OPT-L2: dispatch variant accepting a pre-hoisted AVX2 flag, avoiding a
// repeated RuntimeConfig::get() (TLS read) on each NR iteration.
// CPUID features are invariant for process lifetime; hoisting is safe.
ROBSCALE_HIDDEN inline void bulk_tanh_dispatched(double* inout, int n, bool use_avx2) {
  if (n < 8) {
    for (int i = 0; i < n; ++i) inout[i] = std::tanh(inout[i]);
    return;
  }
#if defined(ROBSCALE_HAS_ACCELERATE)
  (void)use_avx2;
  vvtanh(inout, inout, &n);
#elif defined(ROBSCALE_HAS_SLEEF)
  #ifdef ROBSCALE_HAS_AVX2_DISPATCH
  if (use_avx2) {
    bulk_tanh_sleef_avx2(inout, n);
    return;
  }
  #else
  (void)use_avx2;
  #endif
  for (int i = 0; i < n; i++) inout[i] = std::tanh(inout[i]);
#else
  (void)use_avx2;
  #if defined(_OPENMP) || defined(ROBSCALE_HAS_OMP_SIMD)
    #pragma omp simd
  #endif
  for (int i = 0; i < n; ++i) inout[i] = std::tanh(inout[i]);
#endif
}

// ADM core: constant * mean(|x - center|)
// ROBSCALE_TARGET_AVX2: enables 256-bit AVX2 auto-vectorisation for this
// function without requiring a global -mavx2 flag (safe on CRAN).
// ROBSCALE_RESTRICT: eliminates aliasing analysis so the vectoriser can
// freely reorder loads.
// 4-wide dual-accumulator unroll: breaks the single-accumulator 4-cycle
// latency chain; four independent chains fill one 256-bit AVX2 register
// (4 doubles) per fused-load+abs+add cycle.
ROBSCALE_HIDDEN ROBSCALE_TARGET_AVX2
ROBSCALE_INLINE double adm_core(const double* ROBSCALE_RESTRICT x, int n,
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

// ADM core for pre-sorted input: constant * mean(|x - center|) without abs().
// PRECONDITION: x is sorted ascending AND center = exact median of x.
// Algorithm: sum|x_i - center| = upper_sum - lower_sum (verified formula):
//   lower_sum = sum(x[0..k-1]), upper_sum = sum(x[k+(n&1)..n-1]), k = n/2.
// Two independent SIMD-friendly accumulation loops: no branches, no abs().
// center parameter retained for API symmetry and precondition visibility.
ROBSCALE_HIDDEN ROBSCALE_TARGET_AVX2
ROBSCALE_INLINE double adm_core_sorted(const double* ROBSCALE_RESTRICT x, int n,
                                       double center, double constant) {
  if (ROBSCALE_UNLIKELY(n <= 1)) return 0.0;  // n=0: avoid 1.0/0=Inf; n=1: trivially 0
  (void)center;  // only documents the precondition; not used in computation
  int k = n / 2;
  double lower_sum = 0.0, upper_sum = 0.0;
  for (int i = 0; i < k; ++i)             lower_sum += x[i];
  for (int i = k + (n & 1); i < n; ++i)   upper_sum += x[i];
  return constant * (upper_sum - lower_sum) * (1.0 / n);
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

// Median of pre-sorted array
ROBSCALE_HIDDEN ROBSCALE_INLINE double median_sorted(const double* x, size_t n) {
  if (n & 1) return x[n / 2];
  return (x[(n / 2) - 1] + x[n / 2]) * 0.5;
}

// Selection based median
ROBSCALE_HIDDEN ROBSCALE_INLINE double median_select(double* x, size_t n) {
  if (ROBSCALE_UNLIKELY(n == 0)) return 0.0;
  if (n <= ROBSCALE_SORT_MEDIAN_THRESHOLD) {
    return robscale::median_net(x, n);
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
ROBSCALE_HIDDEN inline double mad_select(const double* x, int n, double med, double* dev) {
  for (int i = 0; i < n; ++i) dev[i] = std::abs(x[i] - med);
  return robscale::MAD_CONSISTENCY * median_select(dev, static_cast<size_t>(n));
}

} // namespace robscale

#endif // ROBSCALE_ROBUST_CORE_H
