#include "robscale_config.h"
#include "robust_core.h"
#include "pdq_select.h"
#include <Rcpp.h>
// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::depends(BH)]]
#include "qnsn_constants.h"
#include "qnsn_sort_utils.h"
#include "qnsn_runtime_config.h"

#include "worker_compat.h"
#include <algorithm>
#include <memory>
#include <type_traits>
#include <cassert>

namespace robscale::qnsn {

constexpr size_t SN_MAX_STACK = ROBSCALE_SN_STACK_THRESHOLD;
// OPT-S3: Micro-buffer threshold for L1-resident fast path (n <= 128 → 1 KB stack frame).
constexpr size_t SN_MICRO_SIZE = 128;

// --- SHARED INNER LOOP ---

// Walking-window inner loop shared by SnWorker, Tier-1, and Tier-2.
// Computes inner_medians[i] = min_L max(sx[i]-sx[L], sx[L+h]-sx[i]) for each i
// in [begin, end), where L walks monotonically forward (amortised O(n) total).
//
// Parameters:
//   sx      - sorted input array (RESTRICT: does not alias im)
//   n       - total number of elements (used to compute L bounds)
//   im      - output array for inner medians (RESTRICT: does not alias sx)
//   begin   - first index to process (inclusive)
//   end     - last index to process (exclusive)
//   L_init  - initial value of the walking pointer L (caller computes via
//             binary search for sub-range chunks; pass 0 for full-range serial)
template <typename T>
#if defined(__GNUC__) || defined(__clang__)
__attribute__((always_inline))
#endif
static inline void sn_inner_serial(
    const T* ROBSCALE_RESTRICT sx,
    int32_t n,
    T* ROBSCALE_RESTRICT im,
    int32_t begin,
    int32_t end,
    int32_t L_init) noexcept {
  int32_t h = static_cast<int32_t>(static_cast<size_t>(n) / 2);
  int32_t L = L_init;
  for (int32_t i = begin; i < end; ++i) {
    int32_t L_min = (std::max)(0, i - h);
    int32_t L_max = (std::min)(i, n - 1 - h);
    if (L < L_min) L = L_min;
    T candidate = (std::max)(sx[i] - sx[L], sx[L + h] - sx[i]);
    while (L < L_max) {
      T next = (std::max)(sx[i] - sx[L + 1], sx[L + 1 + h] - sx[i]);
      if (candidate < next) break;
      L++;
      candidate = next;
    }
    im[i] = candidate;
  }
}

// --- SN ESTIMATOR WORKER ---

// OPT-S4: ROBSCALE_RESTRICT on pointer members: the compiler can now prove
// sorted_x and results do not alias, enabling load hoisting across the inner loop.
template <typename T>
struct SnWorker : public WorkerBase {
  const T* ROBSCALE_RESTRICT sorted_x;
  size_t n;
  mutable T* ROBSCALE_RESTRICT results; // mutable to allow updating results in const operator()

  SnWorker(const T* sorted_x, size_t n, T* results)
      : sorted_x(sorted_x), n(n), results(results) {}

  void operator()(size_t begin, size_t end) {
    if (ROBSCALE_UNLIKELY(begin >= end)) return;
    int32_t h = static_cast<int32_t>(n / 2);

    int32_t i = static_cast<int32_t>(begin);
    int32_t L_min_init = (std::max)(0, i - h);
    int32_t L_max_init = (std::min)(i, static_cast<int32_t>(n) - 1 - h);

    // Initial binary search for the first element in the chunk to find an optimal starting L
    int32_t low = L_min_init;
    int32_t high = L_max_init;
    int32_t L = L_min_init;
    while (low <= high) {
      int32_t mid = low + (high - low) / 2;
      if (mid == L_max_init) {
        L = mid;
        break;
      }
      double v_mid = (std::max)(static_cast<double>(sorted_x[i] - sorted_x[mid]),
                               static_cast<double>(sorted_x[mid + h] - sorted_x[i]));
      double v_next = (std::max)(static_cast<double>(sorted_x[i] - sorted_x[mid + 1]),
                                static_cast<double>(sorted_x[mid + 1 + h] - sorted_x[i]));
      if (v_mid > v_next) {
        low = mid + 1;
        L = mid + 1;
      } else {
        high = mid - 1;
      }
    }

    // Delegate the walking-window loop to the shared template function.
    sn_inner_serial(sorted_x, static_cast<int32_t>(n), results,
                    static_cast<int32_t>(begin), static_cast<int32_t>(end), L);
  }
};

// OPT-S2: Large-n heap path extracted into a NOINLINE function to isolate the
// heap-allocation frame from the hot small-n stack path inside sn_kernel.
// This mirrors the OPT-M1/OPT-I1 pattern: the compiler can fully inline and
// optimise the small-n branch without heap-allocation machinery polluting the
// instruction cache or preventing stack-frame optimisations.
template <typename T>
ROBSCALE_NOINLINE double sn_kernel_large(const T* sorted_x, size_t n) {
  const auto& config = RuntimeConfig::get();

  // Heap path: only inner_medians (n elements, not 2n — sorted_x provided by caller)
  std::unique_ptr<T[]> inner_medians_buf(new T[n]);
  T* inner_medians = inner_medians_buf.get();

  if (n < config.sn_parallel_threshold) {
    SnWorker<T> worker(sorted_x, n, inner_medians);
    worker(0, n);
  } else {
    SnWorker<T> worker(sorted_x, n, inner_medians);
#ifdef USE_DIRECT_TBB
    tbb::parallel_for(tbb::blocked_range<size_t>(0, n, config.grain_size),
                      [&worker](const tbb::blocked_range<size_t>& r) { worker(r.begin(), r.end()); });
#else
    RcppParallel::parallelFor(0, n, worker, config.grain_size);
#endif
  }

  double raw = robscale::adaptive_lowmedian_select(inner_medians, n);
  return raw * CONST_SN * get_sn_factor(n);
}

// Post-sort kernel: sorted_x is read-only, allocates its own inner_medians.
template <typename T>
double sn_kernel(const T* sorted_x, size_t n) {
  const auto& config = RuntimeConfig::get();

  // OPT-S3: Tier 1 — micro path: n <= SN_MICRO_SIZE (128). inner_medians fits in
  // ~1 KB, fully L1-resident. Avoids the 16 KB SN_MAX_STACK frame for small inputs.
  // OPT-S4: RESTRICT aliases allow the compiler to hoist loads and eliminate
  // aliasing barriers between sorted_x reads and inner_medians writes.
  if (n <= SN_MICRO_SIZE) {
    T inner_medians[SN_MICRO_SIZE];
    // OPT-S4: RESTRICT aliases preserve aliasing-free contract required by sn_inner_serial.
    const T* ROBSCALE_RESTRICT sx = sorted_x;
    T* ROBSCALE_RESTRICT im = inner_medians;
    // Tier-1 serial path: full range, L starts at 0.
    sn_inner_serial(sx, static_cast<int32_t>(n), im, 0, static_cast<int32_t>(n), 0);
    // OPT-S5: n <= 16 fast path — small_sort + direct index avoids FR/pdqselect dispatch.
    if (n <= 16) {
      robscale::small_sort(inner_medians, n);
      double raw = static_cast<double>(inner_medians[(n - 1) / 2]);
      return raw * CONST_SN * get_sn_factor(n);
    }
    double raw = robscale::adaptive_lowmedian_select(inner_medians, n);
    return raw * CONST_SN * get_sn_factor(n);
  }

  // Tier 2 — medium stack path: SN_MICRO_SIZE < n <= sn_stack_threshold (2048). 16 KB frame.
  // OPT-S4: RESTRICT aliases allow the compiler to hoist loads and eliminate
  // aliasing barriers between sorted_x reads and inner_medians writes.
  if (n <= config.sn_stack_threshold) {
    assert(config.sn_stack_threshold <= SN_MAX_STACK);
    T inner_medians[SN_MAX_STACK];
    // OPT-S4: RESTRICT aliases preserve aliasing-free contract required by sn_inner_serial.
    const T* ROBSCALE_RESTRICT sx = sorted_x;
    T* ROBSCALE_RESTRICT im = inner_medians;
    // Tier-2 serial path: full range, L starts at 0.
    sn_inner_serial(sx, static_cast<int32_t>(n), im, 0, static_cast<int32_t>(n), 0);
    double raw = robscale::adaptive_lowmedian_select(inner_medians, n);
    return raw * CONST_SN * get_sn_factor(n);
  }

  // Tier 3 — large heap path (n > sn_stack_threshold).
  return sn_kernel_large(sorted_x, n);
}

// OPT-S1: Large-n heap path extracted into a NOINLINE function to isolate the
// large heap-allocation frame from the hot small-n stack path. This mirrors the
// OPT-M1/OPT-I1 pattern: the compiler can now fully inline and optimise the
// small-n branch without the large-n allocation machinery polluting the
// instruction cache or preventing stack-frame optimisations.
template <typename T>
ROBSCALE_NOINLINE double C_sn_impl_large(const T* x_ptr, size_t n) {
  std::unique_ptr<T[]> sorted_buf(new T[n]);
  T* sorted_x = sorted_buf.get();

  for (size_t i = 0; i < n; ++i) {
    if constexpr (std::is_floating_point_v<T>) {
      if (ROBSCALE_UNLIKELY(!std::isfinite(x_ptr[i]))) {
        if (std::isnan(x_ptr[i]))
          Rcpp::stop("There are NAs in the data yet na.rm is FALSE");
        else
          Rcpp::stop("'x' must not contain non-finite values (Inf, -Inf, NaN)");
      }
    }
    sorted_x[i] = x_ptr[i];
  }
  // n > sn_stack_threshold >= SN_MICRO_SIZE (128) > 16, so n <= 16 is
  // structurally unreachable here — optimized_sort is always correct.
  optimized_sort(sorted_x, sorted_x + n);

  return sn_kernel(sorted_x, n);
}

// OPT-S3: Medium stack path (SN_MICRO_SIZE < n <= sn_stack_threshold) extracted
// into a NOINLINE function so the compiler cannot merge the 16 KB SN_MAX_STACK
// frame with the 1 KB SN_MICRO_SIZE frame in the hot micro path.
template <typename T>
ROBSCALE_NOINLINE double C_sn_impl_medium(const T* x_ptr, size_t n) {
  assert(n <= SN_MAX_STACK);
  T sorted_x[SN_MAX_STACK];
  for (size_t i = 0; i < n; ++i) {
    if constexpr (std::is_floating_point_v<T>) {
      if (ROBSCALE_UNLIKELY(!std::isfinite(x_ptr[i]))) {
        if (std::isnan(x_ptr[i]))
          Rcpp::stop("There are NAs in the data yet na.rm is FALSE");
        else
          Rcpp::stop("'x' must not contain non-finite values (Inf, -Inf, NaN)");
      }
    }
    sorted_x[i] = x_ptr[i];
  }
  // n > SN_MICRO_SIZE (128) > 16 — small_sort branch is structurally unreachable.
  optimized_sort(sorted_x, sorted_x + n);
  return sn_kernel(sorted_x, n);
}

template <typename T>
double C_sn_impl(const T* x_ptr, size_t n) {
  if (ROBSCALE_UNLIKELY(n < 2)) return R_NaReal;
  if (ROBSCALE_UNLIKELY(n > 6060000000ULL)) {
    Rcpp::stop("robscale Error: sample size n > 6.06 * 10^9 overflows 64-bit boundaries.");
  }

  const auto& config = RuntimeConfig::get();

  // OPT-S3: Tier 1 — micro path: n <= SN_MICRO_SIZE (128). sorted_x fits in ~1 KB (L1-resident).
  if (n <= SN_MICRO_SIZE) {
    T sorted_x[SN_MICRO_SIZE];
    for (size_t i = 0; i < n; ++i) {
      if constexpr (std::is_floating_point_v<T>) {
        if (ROBSCALE_UNLIKELY(!std::isfinite(x_ptr[i]))) {
          if (std::isnan(x_ptr[i]))
            Rcpp::stop("There are NAs in the data yet na.rm is FALSE");
          else
            Rcpp::stop("'x' must not contain non-finite values (Inf, -Inf, NaN)");
        }
      }
      sorted_x[i] = x_ptr[i];
    }
    if (n <= 16) {
      robscale::small_sort(sorted_x, n);
    } else {
      optimized_sort(sorted_x, sorted_x + n);
    }
    return sn_kernel(sorted_x, n);
  }

  // Tier 2 — medium stack path: SN_MICRO_SIZE < n <= sn_stack_threshold.
  if (n <= config.sn_stack_threshold) {
    assert(config.sn_stack_threshold <= SN_MAX_STACK);
    return C_sn_impl_medium(x_ptr, n);
  }

  return C_sn_impl_large(x_ptr, n);
}

// Sorted variant: input MUST be sorted ascending. No copy, no sort, no NaN scan.
template <typename T>
double C_sn_impl_sorted(const T* sorted_x, size_t n) {
  if (ROBSCALE_UNLIKELY(n < 2)) return R_NaReal;
  assert(std::is_sorted(sorted_x, sorted_x + n));
  return sn_kernel(sorted_x, n);
}

// OPT-S7: Workspace-accepting overload of sn_kernel for large-n heap path.
// When workspace != nullptr AND n > sn_stack_threshold, the provided buffer
// (caller guarantees >= n elements) is used as inner_medians instead of
// heap-allocating. For n <= sn_stack_threshold the workspace is ignored and
// the existing stack path (Tier 1 or Tier 2) is used unchanged.
template <typename T>
double sn_kernel(const T* sorted_x, size_t n, T* workspace) {
  const auto& config = RuntimeConfig::get();
  if (n <= config.sn_stack_threshold) {
    // Tier 1 / Tier 2: stack path unchanged — workspace ignored.
    return sn_kernel(sorted_x, n);
  }
  // Tier 3 heap path: use caller-supplied workspace instead of heap-allocating.
  // workspace != nullptr is a pre-condition for n > sn_stack_threshold callers.
  T* inner_medians = workspace;
  if (n < config.sn_parallel_threshold) {
    SnWorker<T> worker(sorted_x, n, inner_medians);
    worker(0, n);
  } else {
    SnWorker<T> worker(sorted_x, n, inner_medians);
#ifdef USE_DIRECT_TBB
    tbb::parallel_for(tbb::blocked_range<size_t>(0, n, config.grain_size),
                      [&worker](const tbb::blocked_range<size_t>& r) { worker(r.begin(), r.end()); });
#else
    RcppParallel::parallelFor(0, n, worker, config.grain_size);
#endif
  }
  double raw = robscale::adaptive_lowmedian_select(inner_medians, n);
  return raw * CONST_SN * get_sn_factor(n);
}

// OPT-S7: Workspace-accepting overload of C_sn_impl_sorted.
// Delegates to the workspace-accepting sn_kernel. workspace may be nullptr
// only when n <= sn_stack_threshold (stack path ignores it); for n above
// the threshold the caller must supply a buffer of >= n doubles.
// Pre-condition: n >= 2 (guarded by the n < 2 check below).
template <typename T>
double C_sn_impl_sorted(const T* sorted_x, size_t n, T* workspace) {
  if (ROBSCALE_UNLIKELY(n < 2)) return R_NaReal;
  assert(std::is_sorted(sorted_x, sorted_x + n));
  return sn_kernel(sorted_x, n, workspace);
}

// Explicit template instantiations
template double C_sn_impl<double>(const double*, size_t);
template double C_sn_impl_sorted<double>(const double*, size_t);
template double C_sn_impl_sorted<double>(const double*, size_t, double*);

} // namespace robscale::qnsn

// [[Rcpp::export]]
double C_sn_fast(Rcpp::NumericVector x) { 
  return robscale::qnsn::C_sn_impl(x.begin(), static_cast<size_t>(x.size())); 
}

// [[Rcpp::export]]
double C_sn_int_fast(Rcpp::IntegerVector x) { 
  return robscale::qnsn::C_sn_impl(x.begin(), static_cast<size_t>(x.size())); 
}

// [[Rcpp::export]]
double C_get_sn_factor(int n) {
  return robscale::qnsn::get_sn_factor(static_cast<size_t>(n));
}

