#ifndef ROBSCALE_QNSN_THRESHOLDS_H
#define ROBSCALE_QNSN_THRESHOLDS_H

#include <cstddef>

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
constexpr size_t QN_EXACT_THRESHOLD = 64;

// Sn stack threshold: Max size for stack-allocated working array
constexpr size_t SN_STACK_THRESHOLD = 2048;

// Parallel thresholds (serial -> RcppParallel/TBB)
constexpr size_t SN_PARALLEL_THRESHOLD = 12288;
constexpr size_t QN_PARALLEL_THRESHOLD = 8192;


// Sort thresholds
constexpr size_t SORT_BOOST_THRESHOLD = 256;       // std::sort -> Boost spreadsort
constexpr size_t SORT_TBB_FLOAT_THRESHOLD = 6144; // spreadsort -> TBB (float)
constexpr size_t SORT_TBB_INT_THRESHOLD = 8192;   // spreadsort -> TBB (integer)

} // namespace robscale::qnsn

#endif // ROBSCALE_QNSN_THRESHOLDS_H
