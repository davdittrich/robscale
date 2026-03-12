#ifndef ROBSCALE_CONFIG_H
#define ROBSCALE_CONFIG_H

#include <cstddef>

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

// Performance-critical thresholds for Qn/Sn (Aligned to v0.1.5 Gold Standard)
#define ROBSCALE_QN_EXACT_THRESHOLD    64
#define ROBSCALE_SN_STACK_THRESHOLD    2048
#define ROBSCALE_SN_PARALLEL_THRESHOLD 12288
#define ROBSCALE_QN_PARALLEL_THRESHOLD 8192
#define ROBSCALE_SORT_TBB_THRESHOLD    6144
#define ROBSCALE_TBB_GRAIN_SIZE       1024

// Include tuned thresholds if auto-tuning was run (FAST=1 builds)
#ifdef ROBSCALE_HAS_TUNED_THRESHOLDS
#include "qnsn_tuned_thresholds.h"
#endif

// Use tuned value if available, otherwise default to 512
#ifdef TUNED_SORT_BOOST_THRESHOLD
#define ROBSCALE_SORT_BOOST_THRESHOLD  TUNED_SORT_BOOST_THRESHOLD
#else
#define ROBSCALE_SORT_BOOST_THRESHOLD  512
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
  constexpr double ADM_CONSISTENCY      = 1.3926; // Refined value
  constexpr double RHO_SCALE_CONST      = 0.37394112142347236;
  constexpr double INV_RHO_SCALE_CONST  = 1.0 / 0.37394112142347236;
}


#endif // ROBSCALE_CONFIG_H
