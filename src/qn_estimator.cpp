#pragma GCC optimize("O3")
#include <Rcpp.h>
// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::depends(BH)]]
#include "robust_core.h"
#include "qnsn_constants.h"
#include "qnsn_sort_utils.h"
#include "qnsn_dispatcher.h"
#include "qnsn_runtime_config.h"

#include <RcppParallel.h>
#ifdef USE_DIRECT_TBB
#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>
#include <tbb/blocked_range.h>
#endif
#include <algorithm>
#include <memory>
#include <new>
#include <type_traits>
#include <cassert>
#include <cstdlib>

#ifndef NA_REAL
#define NA_REAL R_NaReal
#endif

#ifdef USE_DIRECT_TBB
struct WorkerBase {};
using SplitType = tbb::split;
#else
using WorkerBase = RcppParallel::Worker;
using SplitType = RcppParallel::Split;
#endif

namespace robscale::qnsn {

// --- UTILITIES ---

// Low-median via O(n) selection
template <typename T>
inline double lowmedian_ptr(T* arr, size_t n) {
  if (ROBSCALE_UNLIKELY(n == 0)) return 0.0;
  size_t h = (n - 1) / 2;
  robscale::floyd_rivest_select(arr, arr + h, arr + n);
  return static_cast<double>(arr[h]);
}

// --- QN ESTIMATOR HELPERS ---

template <typename T>
inline T whimed_cpp(T* a, int32_t* iw, size_t n, int64_t target) {
  if (ROBSCALE_UNLIKELY(n == 0)) return T(0);
  if (n == 1) return a[0];

  size_t l = 0, r = n - 1;
  int64_t t = target;

  while (l < r) {
    T pivot = a[l + (r - l) / 2];
    size_t i = l, j = l;
    while (j <= r) {
      if (a[j] < pivot) {
        std::swap(a[i], a[j]);
        std::swap(iw[i], iw[j]);
        i++;
      }
      j++;
    }
    int64_t wleft = 0;
    for (size_t idx = l; idx < i; ++idx) wleft += iw[idx];

    if (wleft > t) {
      r = (i > l) ? i - 1 : l;
    } else {
      size_t i_eq = i, j_eq = i;
      while (j_eq <= r) {
        if (a[j_eq] == pivot) {
          std::swap(a[i_eq], a[j_eq]);
          std::swap(iw[i_eq], iw[j_eq]);
          i_eq++;
        }
        j_eq++;
      }
      int64_t weq = 0;
      for (size_t idx = i; idx < i_eq; ++idx) weq += iw[idx];

      if (wleft + weq > t) return pivot;
      else {
        t -= (wleft + weq);
        l = i_eq;
      }
    }
  }
  return a[l];
}

template <typename T>
struct QnCountWorker : public WorkerBase {
  const T* x;
  size_t n;
  double trial;
  uint64_t sumP = 0;
  uint64_t sumQ = 0;

  QnCountWorker(const T* x, size_t n, double trial)
      : x(x), n(n), trial(trial) {}
  QnCountWorker(const QnCountWorker& other, SplitType)
      : x(other.x), n(other.n), trial(other.trial) {}

#ifdef USE_DIRECT_TBB
  void operator()(const tbb::blocked_range<size_t>& r) {
    this->operator()(r.begin(), r.end());
  }
#endif

  void operator()(size_t begin, size_t end) {
    if (ROBSCALE_UNLIKELY(begin >= end)) return;
    size_t i = begin;
    if (i == 0) i = 1;
    if (ROBSCALE_UNLIKELY(i >= end)) return;

    size_t jp = static_cast<size_t>(std::upper_bound(x, x + i, static_cast<double>(x[i]) - trial) - x);
    size_t jq = static_cast<size_t>(std::lower_bound(x, x + i, static_cast<double>(x[i]) - trial) - x);

    for (; i < end; ++i) {
      double target = static_cast<double>(x[i]) - trial;
      while (ROBSCALE_UNLIKELY(jp < i && static_cast<double>(x[jp]) <= target)) jp++;
      while (ROBSCALE_UNLIKELY(jq < jp && static_cast<double>(x[jq]) < target)) jq++;
      sumP += (i - jp);
      sumQ += (i - jq);
    }
  }

  void join(const QnCountWorker& other) {
    sumP += other.sumP;
    sumQ += other.sumQ;
  }
};

template <typename T>
struct QnRefineWorker : public WorkerBase {
  const T* x;
  size_t n;
  double trial;
  bool is_sumP;
  int32_t* bounds;

  QnRefineWorker(const T* x, size_t n, double trial, bool is_sumP, int32_t* bounds)
      : x(x), n(n), trial(trial), is_sumP(is_sumP), bounds(bounds) {}

#ifdef USE_DIRECT_TBB
  void operator()(const tbb::blocked_range<size_t>& r) const {
    const_cast<QnRefineWorker*>(this)->operator()(r.begin(), r.end());
  }
#endif

  void operator()(size_t begin, size_t end) {
    if (begin >= end) return;
    size_t i = begin;
    if (i == 0) i = 1;
    if (i >= end) return;

    size_t j = is_sumP ? (std::upper_bound(x, x + i, static_cast<double>(x[i]) - trial) - x)
                       : (std::lower_bound(x, x + i, static_cast<double>(x[i]) - trial) - x);

    for (; i < end; ++i) {
      double target = static_cast<double>(x[i]) - trial;
      if (is_sumP) {
        while (ROBSCALE_UNLIKELY(j < i && static_cast<double>(x[j]) <= target)) j++;
        int32_t jj_bound = static_cast<int32_t>(i - j);
        if (jj_bound < bounds[i]) bounds[i] = jj_bound;
      } else {
        while (ROBSCALE_UNLIKELY(j < i && static_cast<double>(x[j]) < target)) j++;
        int32_t jj_bound = static_cast<int32_t>(i - j + 1);
        if (jj_bound > bounds[i]) bounds[i] = jj_bound;
      }
    }
  }
};

template <typename T>
double C_qn_impl(const T* x_ptr, size_t n) {
  if (ROBSCALE_UNLIKELY(n < 2)) return NA_REAL;
  if (ROBSCALE_UNLIKELY(n > 6060000000ULL)) {
    Rcpp::stop("robscale Error: sample size n > 6.06 * 10^9 natively overflows 64-bit boundaries.");
  }

  auto &config = RuntimeConfig::get();
  
  if (n <= config.qn_exact_threshold) {
    size_t num_pairs = n * (n - 1) / 2;
    size_t h_qn = n / 2 + 1;
    size_t k_target = h_qn * (h_qn - 1) / 2;

    std::unique_ptr<T[]> sorted_uninit(new T[n]);
    for (size_t i = 0; i < n; i++) {
      if constexpr (std::is_floating_point_v<T>) {
        if (ROBSCALE_UNLIKELY(!std::isfinite(x_ptr[i]))) return NA_REAL;
      }
      sorted_uninit[i] = x_ptr[i];
    }
    std::sort(sorted_uninit.get(), sorted_uninit.get() + n);

    std::unique_ptr<double[]> diffs(new double[num_pairs]);
    Dispatcher::qn_brute_force(sorted_uninit.get(), n, diffs.get(), config);

    robscale::floyd_rivest_select(diffs.get(), diffs.get() + k_target - 1, diffs.get() + num_pairs);
    double raw = diffs[k_target - 1];
    return raw * CONST_QN * get_qn_factor(n);
  }

  // Properly-typed allocations to avoid alignment UB
  std::unique_ptr<T[]> sorted_x_buf;
  std::unique_ptr<float[]> work_buf;
  std::unique_ptr<int32_t[]> bounds_buf;
  try {
    sorted_x_buf.reset(new T[n]);
    work_buf.reset(new float[n]);
    bounds_buf.reset(new int32_t[3 * n]);
  } catch (const std::bad_alloc& e) {
    Rcpp::stop("robscale Out of Memory: failed to allocate Qn arena.");
  }
  T* sorted_x = sorted_x_buf.get();
  float* work = work_buf.get();
  int32_t* iweight = bounds_buf.get();
  int32_t* left = bounds_buf.get() + n;
  int32_t* right = bounds_buf.get() + 2 * n;

  for (size_t i = 0; i < n; i++) {
    if constexpr (std::is_floating_point_v<T>) {
      if (ROBSCALE_UNLIKELY(!std::isfinite(x_ptr[i]))) return NA_REAL;
    }
    sorted_x[i] = x_ptr[i];
  }
  optimized_sort(sorted_x, sorted_x + n);

  size_t h = n / 2 + 1;
  uint64_t k_target = static_cast<uint64_t>(h) * (h - 1) / 2;

  std::fill_n(left, n, 1);
  for (size_t i = 0; i < n; ++i) right[i] = static_cast<int32_t>(i);

  uint64_t nL = 0;
  uint64_t nR = static_cast<uint64_t>(n) * (n - 1) / 2;

  while (nR - nL > n) {
    size_t m = 0;
    for (size_t i = 1; i < n; ++i) {
      if (left[i] <= right[i]) {
        int32_t w = right[i] - left[i] + 1;
        int32_t jj = left[i] + w / 2;
        work[m] = static_cast<float>(static_cast<double>(sorted_x[i]) - static_cast<double>(sorted_x[i - jj]));
        iweight[m] = w;
        m += 1;
      }
    }

    double trial = whimed_cpp(work, iweight, m, static_cast<int64_t>((nR - nL) / 2));

    QnCountWorker<T> countWorker(sorted_x, n, trial);
    if (n > config.qn_parallel_threshold) {
#ifdef USE_DIRECT_TBB
      tbb::parallel_reduce(tbb::blocked_range<size_t>(1, n, 1024), countWorker);
#else
      RcppParallel::parallelReduce(1, n, countWorker, 1024);
#endif
    } else {
      countWorker(1, n);
    }

    if (k_target <= countWorker.sumP) {
      QnRefineWorker<T> refineWorker(sorted_x, n, trial, true, right);
      if (n > config.qn_parallel_threshold) {
#ifdef USE_DIRECT_TBB
        tbb::parallel_for(tbb::blocked_range<size_t>(1, n, 1024), refineWorker);
#else
        RcppParallel::parallelFor(1, n, refineWorker, 1024);
#endif
      } else refineWorker(1, n);
      nR = countWorker.sumP;
    } else if (k_target > countWorker.sumQ) {
      QnRefineWorker<T> refineWorker(sorted_x, n, trial, false, left);
      if (n > config.qn_parallel_threshold) {
#ifdef USE_DIRECT_TBB
        tbb::parallel_for(tbb::blocked_range<size_t>(1, n, 1024), refineWorker);
#else
        RcppParallel::parallelFor(1, n, refineWorker, 1024);
#endif
      } else refineWorker(1, n);
      nL = countWorker.sumQ;
    } else {
      return trial * CONST_QN * get_qn_factor(n);
    }
  }

  size_t num_final = static_cast<size_t>(nR - nL);
  std::unique_ptr<double[]> final_diffs(new double[num_final]);
  size_t fd_idx = 0;
  for (size_t i = 1; i < n; ++i) {
    for (int32_t jj = left[i]; jj <= right[i]; ++jj) {
      final_diffs[fd_idx++] = static_cast<double>(sorted_x[i]) - static_cast<double>(sorted_x[i - jj]);
    }
  }
  
  double* k_ptr = final_diffs.get() + (k_target - nL - 1);
  robscale::floyd_rivest_select(final_diffs.get(), k_ptr, final_diffs.get() + num_final);
  double raw = *k_ptr;
  return raw * CONST_QN * get_qn_factor(n);
}

} // namespace robscale::qnsn

// --- R EXPORTS ---



// [[Rcpp::export]]
double C_qn_fast(Rcpp::NumericVector x) { 
  return robscale::qnsn::C_qn_impl(x.begin(), static_cast<size_t>(x.size())); 
}

// [[Rcpp::export]]
double C_qn_int_fast(Rcpp::IntegerVector x) { 
  return robscale::qnsn::C_qn_impl(x.begin(), static_cast<size_t>(x.size())); 
}

// [[Rcpp::export]]
Rcpp::List get_qnsn_config() {
  auto &config = robscale::qnsn::RuntimeConfig::get();
  return Rcpp::List::create(
    Rcpp::Named("simd_level") = (int)config.hw.simd_level,
    Rcpp::Named("qn_exact_threshold") = (int)config.qn_exact_threshold,
    Rcpp::Named("qn_parallel_threshold") = (int)config.qn_parallel_threshold,
    Rcpp::Named("sn_stack_threshold") = (int)config.sn_stack_threshold,
    Rcpp::Named("sn_parallel_threshold") = (int)config.sn_parallel_threshold,
    Rcpp::Named("l2_cache_size") = (int)config.hw.l2_cache_size,
    Rcpp::Named("num_logical_cores") = (int)config.hw.num_logical_cores
  );
}



// [[Rcpp::export]]
double C_get_qn_factor(int n) {
  return robscale::qnsn::get_qn_factor(static_cast<size_t>(n));
}
