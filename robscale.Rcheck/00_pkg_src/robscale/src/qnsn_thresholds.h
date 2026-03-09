#ifndef ROBSCALE_QNSN_THRESHOLDS_H
#define ROBSCALE_QNSN_THRESHOLDS_H

#include <cstddef>

#ifdef ROBSCALE_HAS_TUNED_THRESHOLDS
#include "qnsn_tuned_thresholds.h"
#endif

// Cache-aware threshold constants for robscale::qnsn
#ifndef FASTQNSN_L2_CACHE_BYTES
#define FASTQNSN_L2_CACHE_BYTES 4194304 // 4 MB conservative default
#endif

#ifndef FASTQNSN_CACHE_LINE_BYTES
#ifdef __aarch64__
#define FASTQNSN_CACHE_LINE_BYTES 128 // ARM / Apple Silicon
#else
#define FASTQNSN_CACHE_LINE_BYTES 64 // x86-64
#endif
#endif

namespace robscale::qnsn {

// Brute-force threshold for Qn. Our JM implementation is very fast, so we crossover earlier.
#ifdef TUNED_QN_EXACT_THRESHOLD
constexpr size_t QN_EXACT_THRESHOLD = TUNED_QN_EXACT_THRESHOLD;
#else
constexpr size_t QN_EXACT_THRESHOLD = 64;
#endif

// Sn stack threshold: Max size for stack-allocated working array
constexpr size_t SN_STACK_THRESHOLD = 2048;

// Parallel thresholds (serial -> RcppParallel/TBB)
constexpr size_t SN_PARALLEL_THRESHOLD = 12288;
constexpr size_t QN_PARALLEL_THRESHOLD = 8192;


// Sort thresholds
#ifdef TUNED_SORT_BOOST_THRESHOLD
constexpr size_t SORT_BOOST_THRESHOLD = TUNED_SORT_BOOST_THRESHOLD;
#else
constexpr size_t SORT_BOOST_THRESHOLD = 512;      // std::sort -> Boost spreadsort (safe default)
#endif
constexpr size_t SORT_TBB_FLOAT_THRESHOLD = 6144; // spreadsort -> TBB (float)
constexpr size_t SORT_TBB_INT_THRESHOLD = 8192;   // spreadsort -> TBB (integer)

} // namespace robscale::qnsn

#endif // ROBSCALE_QNSN_THRESHOLDS_H
