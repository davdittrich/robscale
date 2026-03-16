#include "robscale_config.h"
#include "robust_core.h"
#include "pdq_select.h"
#include <Rcpp.h>
// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::depends(BH)]]
#include "qnsn_constants.h"
#include "qnsn_sort_utils.h"
#include "qnsn_runtime_config.h"

#include <RcppParallel.h>
#ifdef USE_DIRECT_TBB
// Include TBB headers
#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>
#include <tbb/blocked_range.h>
#endif
#include <algorithm>
#include <memory>
#include <type_traits>
#include <cassert>

#ifdef USE_DIRECT_TBB
struct WorkerBase {};
using SplitType = tbb::split;
#else
using WorkerBase = RcppParallel::Worker;
using SplitType = RcppParallel::Split;
#endif

namespace robscale::qnsn {

constexpr size_t SN_MAX_STACK = ROBSCALE_SN_STACK_THRESHOLD;

// --- SN ESTIMATOR WORKER ---

template <typename T>
struct SnWorker : public WorkerBase {
  const T* sorted_x;
  size_t n;
  mutable T* results; // mutable to allow updating results in const operator()

  SnWorker(const T* sorted_x, size_t n, T* results)
      : sorted_x(sorted_x), n(n), results(results) {}

#ifdef USE_DIRECT_TBB
  void operator()(const tbb::blocked_range<size_t>& r) const {
    const_cast<SnWorker*>(this)->operator()(r.begin(), r.end());
  }
#endif

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

    for (; i < static_cast<int32_t>(end); ++i) {
      int32_t L_min = (std::max)(0, i - h);
      int32_t L_max = (std::min)(i, static_cast<int32_t>(n) - 1 - h);
      if (L < L_min) L = L_min;

      T candidate = (std::max)(sorted_x[i] - sorted_x[L], sorted_x[L + h] - sorted_x[i]);
      while (L < L_max) {
        T next = (std::max)(sorted_x[i] - sorted_x[L + 1], sorted_x[L + 1 + h] - sorted_x[i]);
        if (candidate < next) break;
        L++;
        candidate = next;
      }
      results[i] = candidate;
    }
  }
};

// Post-sort kernel: sorted_x is read-only, allocates its own inner_medians.
template <typename T>
double sn_kernel(const T* sorted_x, size_t n) {
  const auto& config = RuntimeConfig::get();

  if (n <= config.sn_stack_threshold) {
    assert(config.sn_stack_threshold <= SN_MAX_STACK);
    T inner_medians[SN_MAX_STACK];
    int32_t h = static_cast<int32_t>(n / 2);
    int32_t L = 0;
    for (int32_t i = 0; i < static_cast<int32_t>(n); ++i) {
      int32_t L_min = (std::max)(0, i - h);
      int32_t L_max = (std::min)(i, static_cast<int32_t>(n) - 1 - h);
      if (L < L_min) L = L_min;
      T candidate = (std::max)(sorted_x[i] - sorted_x[L], sorted_x[L + h] - sorted_x[i]);
      while (L < L_max) {
        T next = (std::max)(sorted_x[i] - sorted_x[L + 1], sorted_x[L + 1 + h] - sorted_x[i]);
        if (candidate < next) break;
        L++;
        candidate = next;
      }
      inner_medians[i] = candidate;
    }
    double raw = robscale::adaptive_lowmedian_select(inner_medians, n);
    return raw * CONST_SN * get_sn_factor(n);
  }

  // Heap path: only inner_medians (n elements, not 2n — sorted_x provided by caller)
  std::unique_ptr<T[]> inner_medians_buf(new T[n]);
  T* inner_medians = inner_medians_buf.get();

  if (n < config.sn_parallel_threshold) {
    SnWorker<T> worker(sorted_x, n, inner_medians);
    worker(0, n);
  } else {
    SnWorker<T> worker(sorted_x, n, inner_medians);
#ifdef USE_DIRECT_TBB
    tbb::parallel_for(tbb::blocked_range<size_t>(0, n, config.grain_size), worker);
#else
    RcppParallel::parallelFor(0, n, worker, config.grain_size);
#endif
  }

  double raw = robscale::adaptive_lowmedian_select(inner_medians, n);
  return raw * CONST_SN * get_sn_factor(n);
}

template <typename T>
double C_sn_impl(const T* x_ptr, size_t n) {
  if (ROBSCALE_UNLIKELY(n < 2)) return R_NaReal;
  if (ROBSCALE_UNLIKELY(n > 6060000000ULL)) {
    Rcpp::stop("robscale Error: sample size n > 6.06 * 10^9 overflows 64-bit boundaries.");
  }

  const auto& config = RuntimeConfig::get();

  if (n <= config.sn_stack_threshold) {
    assert(config.sn_stack_threshold <= SN_MAX_STACK);
    T sorted_x[SN_MAX_STACK];
    for (size_t i = 0; i < n; ++i) {
      if constexpr (std::is_floating_point_v<T>) {
        if (ROBSCALE_UNLIKELY(!std::isfinite(x_ptr[i]))) return R_NaReal;
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

  // Heap path: copy + sort, then delegate to kernel
  std::unique_ptr<T[]> sorted_buf(new T[n]);
  T* sorted_x = sorted_buf.get();

  for (size_t i = 0; i < n; ++i) {
    if constexpr (std::is_floating_point_v<T>) {
      if (ROBSCALE_UNLIKELY(!std::isfinite(x_ptr[i]))) return R_NaReal;
    }
    sorted_x[i] = x_ptr[i];
  }
  optimized_sort(sorted_x, sorted_x + n);

  return sn_kernel(sorted_x, n);
}

// Sorted variant: input MUST be sorted ascending. No copy, no sort, no NaN scan.
template <typename T>
double C_sn_impl_sorted(const T* sorted_x, size_t n) {
  if (ROBSCALE_UNLIKELY(n < 2)) return R_NaReal;
  assert(std::is_sorted(sorted_x, sorted_x + n));
  return sn_kernel(sorted_x, n);
}

// Explicit template instantiations
template double C_sn_impl<double>(const double*, size_t);
template double C_sn_impl<float>(const float*, size_t);
template double C_sn_impl_sorted<double>(const double*, size_t);

} // namespace robscale::qnsn

// [[Rcpp::export]]
double C_sn_fast(Rcpp::NumericVector x) { 
  return robscale::qnsn::C_sn_impl(x.begin(), static_cast<size_t>(x.size())); 
}

// [[Rcpp::export]]
double C_sn_float(Rcpp::NumericVector x) {
  std::vector<float> x_float(x.begin(), x.end());
  return robscale::qnsn::C_sn_impl(x_float.data(), x_float.size());
}

// [[Rcpp::export]]
double C_sn_int_fast(Rcpp::IntegerVector x) { 
  return robscale::qnsn::C_sn_impl(x.begin(), static_cast<size_t>(x.size())); 
}

// [[Rcpp::export]]
double C_get_sn_factor(int n) {
  return robscale::qnsn::get_sn_factor(static_cast<size_t>(n));
}
