#ifndef ROBSCALE_CONFIG_H
#define ROBSCALE_CONFIG_H

#include <cstddef>
#include <cmath>
#include <memory>

/**
 * robscale_config.h
 * 
 * Centralized architectural thresholds and portability macros for the robscale package.
 * Decouples performance logic from compiler-specific pragmas.
 */


// --- Performance Thresholds ---

/**
 * Threshold for switching from selection (Floyd-Rivest) to full sorting.
 * Sorting is more branch-predictable and often faster for small arrays.
 */
#ifndef ROBSCALE_SORT_THRESHOLD
#define ROBSCALE_SORT_THRESHOLD 512
#endif

/**
 * Threshold for surgical micro-dispatchers (zero-allocation fast paths).
 * These paths use stack-allocated buffers and bypass Rcpp boundary overhead.
 */

// Fixed thresholds (not cache-sensitive; cache-sensitive thresholds are
// derived from hw.l2_per_core at runtime in RuntimeConfig::calculate_thresholds)
#define ROBSCALE_QN_EXACT_THRESHOLD    64
#define ROBSCALE_SN_STACK_THRESHOLD    2048
#define ROBSCALE_SORT_BOOST_THRESHOLD  512

/**
 * Threshold for median_select: use small_sort (NBIS) for n <= this value,
 * floyd_rivest_select for n > this value. Empirically calibrated: NBIS
 * (branchless O(n²)) beats floyd_rivest_select for small n because NBIS has
 * zero branch mispredictions and full L1 residency. Set conservatively at 32;
 * raise after benchmarking Phase 3+4 of the optimization plan.
 */
#ifndef ROBSCALE_SORT_MEDIAN_THRESHOLD
#define ROBSCALE_SORT_MEDIAN_THRESHOLD 64
#endif


// --- Portability & Optimization Macros ---

/**
 * Branch prediction hints.
 */
#if defined(__GNUC__) || defined(__clang__)
#define ROBSCALE_LIKELY(x)   __builtin_expect(!!(x), 1)
#define ROBSCALE_UNLIKELY(x) __builtin_expect(!!(x), 0)
#else
#define ROBSCALE_LIKELY(x)   (x)
#define ROBSCALE_UNLIKELY(x) (x)
#endif

/**
 * Force inline for micro-kernels.
 * Relaxed to standard inline to allow compiler better optimization for tiny-n paths.
 */
#define ROBSCALE_INLINE inline

/**
 * Buffer management for micro-samples.
 * Ensures enough space for internal calculations (e.g. absolute differences).
 */
#define ROBSCALE_MICRO_BUFFER_SIZE 128

/**
 * Prevent inlining so small-n paths get their own minimal stack frame.
 */
#if defined(_MSC_VER)
  #define ROBSCALE_NOINLINE __declspec(noinline)
#elif defined(__GNUC__) || defined(__clang__)
  #define ROBSCALE_NOINLINE __attribute__((noinline))
#else
  #define ROBSCALE_NOINLINE
#endif


// --- Common Statistical Constants ---

#if defined(_MSC_VER)
  #define ROBSCALE_RESTRICT __restrict
#elif defined(__GNUC__) || defined(__clang__)
  #define ROBSCALE_RESTRICT __restrict__
#else
  #define ROBSCALE_RESTRICT
#endif

namespace robscale {
  constexpr double MAD_CONSISTENCY      = 1.482602218505602;
  constexpr double ADM_CONSISTENCY      = 1.2533141373155001; // sqrt(pi/2), matches adm() default
  constexpr double GMD_CONSISTENCY      = 0.886226925452758;  // sqrt(pi)/2
  constexpr double IQR_CONSISTENCY      = 0.741301109252801;  // 1/(Phi^-1(0.75) - Phi^-1(0.25))
  constexpr double RHO_SCALE_CONST      = 0.37394112142347236;
  constexpr double INV_RHO_SCALE_CONST  = 1.0 / RHO_SCALE_CONST;  // N3

  /// MAD implosion threshold for M-scale and ensemble estimators.
  /// When MAD(x) <= IMPLOSION_BOUND, more than 50% of observations are tied
  /// and the scale is degenerate; triggers ADM fallback.
  /// Value 1e-4: small enough to be effectively zero for practical scales,
  /// but above floating-point noise for data with legitimate small spread.
  /// Ref: Rousseeuw & Verboven (2002), Sec. 4.3.
  constexpr double IMPLOSION_BOUND = 1e-4;

  // c4(n) consistency constant for unbiased standard deviation
  // Formula: sqrt(2/(n-1)) * Gamma(n/2) / Gamma((n-1)/2)
  inline double c4_factor(int n) {
    if (n < 2) return 1.0;
    return std::exp(0.5 * std::log(2.0 / (n - 1.0))
                    + std::lgamma(n / 2.0)
                    - std::lgamma((n - 1.0) / 2.0));
  }
}


/// Tiered stack/heap scratch buffer.
/// N_MICRO: fits in L1 (hot path for very small n, zero malloc overhead).
/// N_STACK: fits in L2 (medium n, no heap allocation).
/// Falls back to heap for n > N_STACK.
template <size_t N_MICRO = 128, size_t N_STACK = 2048>
struct StackArena {
  double buf_micro[N_MICRO];
  double buf_stack[N_STACK];
  std::unique_ptr<double[]> heap;

  double* get(size_t n) {
    if (n <= N_MICRO) return buf_micro;
    if (n <= N_STACK) return buf_stack;
    heap.reset(new double[n]);
    return heap.get();
  }
};

// --- Runtime SIMD dispatch ---
// Per-function target attributes: the compiler emits AVX2/FMA instructions
// for annotated functions without requiring -mavx2 globally.  This produces
// a portable binary that activates AVX2 at runtime on capable hardware.
// Supported by GCC 4.9+, all Clang, all Rtools MinGW-w64 GCC.
#if (defined(__x86_64__) || defined(_M_X64)) && \
    (defined(__GNUC__) || defined(__clang__))
  #define ROBSCALE_TARGET_AVX2  __attribute__((target("avx2,fma")))
  #define ROBSCALE_HAS_AVX2_DISPATCH 1
#else
  #define ROBSCALE_TARGET_AVX2
#endif

#endif // ROBSCALE_CONFIG_H
