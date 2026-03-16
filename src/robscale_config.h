#ifndef ROBSCALE_CONFIG_H
#define ROBSCALE_CONFIG_H

#include <cstddef>
#include <cmath>

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
  constexpr double INV_RHO_SCALE_CONST  = 1.0 / 0.37394112142347236;

  // c4(n) consistency constant for unbiased standard deviation
  // Formula: sqrt(2/(n-1)) * Gamma(n/2) / Gamma((n-1)/2)
  inline double c4_factor(int n) {
    if (n < 2) return 1.0;
    return std::exp(0.5 * std::log(2.0 / (n - 1.0))
                    + std::lgamma(n / 2.0)
                    - std::lgamma((n - 1.0) / 2.0));
  }
}


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
