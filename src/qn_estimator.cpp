#include "robscale_config.h"
#include "robust_core.h"
#include "pdq_select.h"
#include <Rcpp.h>
// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::depends(BH)]]
#include "qnsn_constants.h"
#include "qnsn_sort_utils.h"
#include "qnsn_kernels.h"

#include "qnsn_dispatcher.h"
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
#include <vector>

#ifdef USE_DIRECT_TBB
struct WorkerBase {};
using SplitType = tbb::split;
#else
using WorkerBase = RcppParallel::Worker;
using SplitType = RcppParallel::Split;
#endif

namespace robscale::qnsn {

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
  uint64_t sumP;
  uint64_t sumQ;

  QnCountWorker(const T* x, size_t n, double trial)
      : x(x), n(n), trial(trial), sumP(0), sumQ(0) {}

  QnCountWorker(const QnCountWorker& other, SplitType)
      : x(other.x), n(other.n), trial(other.trial), sumP(0), sumQ(0) {}

#ifdef USE_DIRECT_TBB
  void operator()(const tbb::blocked_range<size_t>& r) {
    operator()(r.begin(), r.end());
  }
#endif

  void operator()(size_t begin, size_t end) {
    if (begin >= end) return;
    size_t i = begin;
    if (i == 0) i = 1;
    if (i >= end) return;

    size_t jP = std::upper_bound(x, x + i, static_cast<double>(x[i]) - trial) - x;
    size_t jQ = std::lower_bound(x, x + i, static_cast<double>(x[i]) - trial) - x;

    for (; i < end; ++i) {
      double target = static_cast<double>(x[i]) - trial;
      while (jP < i && static_cast<double>(x[jP]) <= target) jP++;
      while (jQ < i && static_cast<double>(x[jQ]) < target) jQ++;
      sumP += (i - jP);
      sumQ += (i - jQ);
    }
  }

  void join(const QnCountWorker& other) {
    sumP += other.sumP;
    sumQ += other.sumQ;
  }
};

// QnCandidateWorker and QnCandidateFulfiller removed — the struct-based
// blocked_range pattern used r.begin()/grain as a block index, which is
// unsafe because TBB does not guarantee range splits at grain multiples.
// Replaced with explicit block-index parallel_for in qn_refinement_kernel.

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

// Post-sort brute-force kernel: sorted_x is read-only, n <= 64.
template <typename T>
double qn_brute_force_kernel(const T* sorted_x, size_t n) {
  size_t num_pairs = n * (n - 1) / 2;
  size_t h_qn = n / 2 + 1;
  size_t k_target = h_qn * (h_qn - 1) / 2;

  constexpr size_t QN_MAX_PAIRS =
      ROBSCALE_QN_EXACT_THRESHOLD * (ROBSCALE_QN_EXACT_THRESHOLD - 1) / 2;
  double diffs_buf[QN_MAX_PAIRS];

  const auto& config = RuntimeConfig::get();
  if (n > 16) {
    Dispatcher::qn_brute_force(sorted_x, n, diffs_buf, config);
  } else {
    qn_brute_force_scalar(sorted_x, n, diffs_buf);
  }

  robscale::floyd_rivest_select(diffs_buf, diffs_buf + k_target - 1, diffs_buf + num_pairs);
  return diffs_buf[k_target - 1] * CONST_QN * get_qn_factor(n);
}

// Original entry point: validate + copy + sort + kernel
template <typename T>
double qn_brute_force_exact(const T* x_ptr, size_t n) {
  T sorted_buf[64];
  for (size_t i = 0; i < n; i++) {
    if constexpr (std::is_floating_point_v<T>) {
      if (ROBSCALE_UNLIKELY(!std::isfinite(x_ptr[i]))) return R_NaReal;
    }
    sorted_buf[i] = x_ptr[i];
  }
  if (n <= 16) {
    robscale::small_sort(sorted_buf, n);
  } else {
    std::sort(sorted_buf, sorted_buf + n);
  }
  return qn_brute_force_kernel(sorted_buf, n);
}

// Post-sort refinement kernel: sorted_x is read-only, n > qn_exact_threshold.
template <typename T>
double qn_refinement_kernel(const T* sorted_x, size_t n) {
  const auto& config = RuntimeConfig::get();

  std::unique_ptr<float[]> work_buf;
  std::unique_ptr<int32_t[]> bounds_buf;
  try {
    work_buf.reset(new float[n]);
    bounds_buf.reset(new int32_t[3 * n]);
  } catch (const std::bad_alloc& e) {
    Rcpp::stop("robscale Out of Memory: failed to allocate Qn arena.");
  }
  float* work = work_buf.get();
  int32_t* iweight = bounds_buf.get();
  int32_t* left = bounds_buf.get() + n;
  int32_t* right = bounds_buf.get() + 2 * n;

  size_t h = n / 2 + 1;
  uint64_t k_target = static_cast<uint64_t>(h) * (h - 1) / 2;

  std::fill_n(left, n, 1);
  for (size_t i = 0; i < n; ++i) right[i] = static_cast<int32_t>(i);

  uint64_t nL = 0;
  uint64_t nR = static_cast<uint64_t>(n) * (n - 1) / 2;

  while (nR - nL > n) {
    // --- Candidate generation: collect weighted medians ---
    size_t m = 0;
    bool used_parallel = false;
#ifdef USE_DIRECT_TBB
    if (n > config.qn_parallel_threshold) {
      // Explicit block-index parallel_for: each block_idx maps to a fixed
      // [begin, end) range, so count phase and fill phase see identical splits.
      size_t g = config.grain_size;
      size_t num_blocks = (n - 1 + g - 1) / g;  // blocks covering [1, n)
      std::vector<size_t> block_offsets(num_blocks, 0);

      // Count: how many candidates per block
      tbb::parallel_for(size_t(0), num_blocks, [&](size_t block_idx) {
        size_t begin = 1 + block_idx * g;
        size_t end = std::min(begin + g, n);
        size_t count = 0;
        for (size_t i = begin; i < end; ++i) {
          if (left[i] <= right[i]) count++;
        }
        block_offsets[block_idx] = count;
      });

      // Prefix sum -> write offsets
      size_t current = 0;
      for (size_t b = 0; b < num_blocks; ++b) {
        size_t c = block_offsets[b];
        block_offsets[b] = current;
        current += c;
      }
      m = current;

      // Fill: each block writes at its deterministic offset
      tbb::parallel_for(size_t(0), num_blocks, [&](size_t block_idx) {
        size_t begin = 1 + block_idx * g;
        size_t end = std::min(begin + g, n);
        size_t o = block_offsets[block_idx];
        for (size_t i = begin; i < end; ++i) {
          if (left[i] <= right[i]) {
            int32_t w_val = right[i] - left[i] + 1;
            int32_t jj = left[i] + w_val / 2;
            work[o] = static_cast<float>(static_cast<double>(sorted_x[i]) - static_cast<double>(sorted_x[i - jj]));
            iweight[o] = w_val;
            o++;
          }
        }
      });
      used_parallel = true;
    }
#endif
    if (!used_parallel) {
      for (size_t i = 1; i < n; ++i) {
        if (left[i] <= right[i]) {
          int32_t w_val = right[i] - left[i] + 1;
          int32_t jj = left[i] + w_val / 2;
          work[m] = static_cast<float>(static_cast<double>(sorted_x[i]) - static_cast<double>(sorted_x[i - jj]));
          iweight[m] = w_val;
          m += 1;
        }
      }
    }

    double trial = whimed_cpp(work, iweight, m, static_cast<int64_t>((nR - nL) / 2));

    // --- Count P/Q sums ---
    QnCountWorker<T> countWorker(sorted_x, n, trial);
    if (n > config.qn_parallel_threshold) {
#ifdef USE_DIRECT_TBB
      tbb::parallel_reduce(tbb::blocked_range<size_t>(1, n, config.get_dynamic_grain_size(n)), countWorker);
#else
      RcppParallel::parallelReduce(1, n, countWorker, config.get_dynamic_grain_size(n));
#endif
    } else {
      countWorker(1, n);
    }

    // --- Refine bounds ---
    if (k_target <= countWorker.sumP) {
      QnRefineWorker<T> refineWorker(sorted_x, n, trial, true, right);
      if (n > config.qn_parallel_threshold) {
        size_t g = config.get_dynamic_grain_size(n);
#ifdef USE_DIRECT_TBB
        tbb::parallel_for(tbb::blocked_range<size_t>(1, n, g), refineWorker);
#else
        RcppParallel::parallelFor(1, n, refineWorker, g);
#endif
      } else refineWorker(1, n);
      nR = countWorker.sumP;
    } else if (k_target > countWorker.sumQ) {
      QnRefineWorker<T> refineWorker(sorted_x, n, trial, false, left);
      if (n > config.qn_parallel_threshold) {
#ifdef USE_DIRECT_TBB
        tbb::parallel_for(tbb::blocked_range<size_t>(1, n, config.get_dynamic_grain_size(n)), refineWorker);
#else
        RcppParallel::parallelFor(1, n, refineWorker, config.get_dynamic_grain_size(n));
#endif
      } else refineWorker(1, n);
      nL = countWorker.sumQ;
    } else {
      return trial * CONST_QN * get_qn_factor(n);
    }
  }

  // --- Final selection on surviving diffs ---
  std::unique_ptr<double[]> diffs(new double[nR - nL]);

  bool filled_parallel = false;
#ifdef USE_DIRECT_TBB
  if (n > config.qn_parallel_threshold) {
    // Explicit block-index parallel_for — same pattern as candidate generation.
    size_t g = config.get_dynamic_grain_size(n);
    size_t num_blocks = (n - 1 + g - 1) / g;  // blocks covering [1, n)
    std::vector<size_t> block_offsets(num_blocks, 0);

    // Count: total diffs per block
    tbb::parallel_for(size_t(0), num_blocks, [&](size_t block_idx) {
      size_t begin = 1 + block_idx * g;
      size_t end = std::min(begin + g, n);
      size_t count = 0;
      for (size_t i = begin; i < end; ++i) {
        if (left[i] <= right[i]) count += (right[i] - left[i] + 1);
      }
      block_offsets[block_idx] = count;
    });

    // Prefix sum -> write offsets
    size_t current = 0;
    for (size_t b = 0; b < num_blocks; ++b) {
      size_t c = block_offsets[b];
      block_offsets[b] = current;
      current += c;
    }

    // Fill: each block writes diffs at its deterministic offset
    tbb::parallel_for(size_t(0), num_blocks, [&](size_t block_idx) {
      size_t begin = 1 + block_idx * g;
      size_t end = std::min(begin + g, n);
      size_t o = block_offsets[block_idx];
      for (size_t i = begin; i < end; ++i) {
        if (left[i] <= right[i]) {
          for (int32_t j = left[i]; j <= right[i]; ++j) {
            diffs[o++] = static_cast<double>(sorted_x[i]) - static_cast<double>(sorted_x[i - j]);
          }
        }
      }
    });
    filled_parallel = true;
  }
#endif
  if (!filled_parallel) {
    size_t offset = 0;
    for (size_t i = 1; i < n; ++i) {
      for (int32_t j = left[i]; j <= right[i]; ++j) {
        diffs[offset++] = static_cast<double>(sorted_x[i]) - static_cast<double>(sorted_x[i - j]);
      }
    }
  }

  {
    const size_t window = static_cast<size_t>(nR - nL);
    double* sel_ptr = diffs.get();
    const size_t k_off = static_cast<size_t>(k_target - nL - 1);
    if (window <= config.pdq_qn_final_threshold)
      robscale::floyd_rivest_select(sel_ptr, sel_ptr + k_off, sel_ptr + window);
    else
      miniselect::pdqselect(sel_ptr, sel_ptr + k_off, sel_ptr + window);
  }
  double final_raw = diffs[k_target - nL - 1];
  return final_raw * CONST_QN * get_qn_factor(n);
}

template <typename T>
double C_qn_impl(const T* x_ptr, size_t n) {
  if (ROBSCALE_UNLIKELY(n < 2)) return R_NaReal;
  if (ROBSCALE_UNLIKELY(n > 6060000000ULL)) {
    Rcpp::stop("robscale Error: sample size n > 6.06 * 10^9 natively overflows 64-bit boundaries.");
  }

  const auto& config = RuntimeConfig::get();

  if (n <= config.qn_exact_threshold) {
    return qn_brute_force_exact(x_ptr, n);
  }

  // Copy + sort, then delegate to refinement kernel
  std::unique_ptr<T[]> sorted_x_buf;
  try {
    sorted_x_buf.reset(new T[n]);
  } catch (const std::bad_alloc& e) {
    Rcpp::stop("robscale Out of Memory: failed to allocate Qn sorted buffer.");
  }
  T* sorted_x = sorted_x_buf.get();

  for (size_t i = 0; i < n; i++) {
    if constexpr (std::is_floating_point_v<T>) {
      if (ROBSCALE_UNLIKELY(!std::isfinite(x_ptr[i]))) return R_NaReal;
    }
    sorted_x[i] = x_ptr[i];
  }
  optimized_sort(sorted_x, sorted_x + n);

  return qn_refinement_kernel(sorted_x, n);
}

// Sorted variant: input MUST be sorted ascending. No copy, no sort, no NaN scan.
template <typename T>
double C_qn_impl_sorted(const T* sorted_x, size_t n) {
  if (ROBSCALE_UNLIKELY(n < 2)) return R_NaReal;
  assert(std::is_sorted(sorted_x, sorted_x + n));

  const auto& config = RuntimeConfig::get();
  if (n <= config.qn_exact_threshold) {
    return qn_brute_force_kernel(sorted_x, n);
  }
  return qn_refinement_kernel(sorted_x, n);
}

// Explicit template instantiations
template double C_qn_impl_sorted<double>(const double*, size_t);
template double C_qn_impl<double>(const double*, size_t);
template double C_qn_impl<int>(const int*, size_t);

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
  const auto& config = robscale::qnsn::RuntimeConfig::get();
  return Rcpp::List::create(
    Rcpp::Named("simd_level") = (int)config.hw.simd_level,
    Rcpp::Named("qn_exact_threshold") = (int)config.qn_exact_threshold,
    Rcpp::Named("qn_parallel_threshold") = (int)config.qn_parallel_threshold,
    Rcpp::Named("sn_stack_threshold") = (int)config.sn_stack_threshold,
    Rcpp::Named("sn_parallel_threshold") = (int)config.sn_parallel_threshold,
    Rcpp::Named("sort_boost_threshold") = (int)config.sort_boost_threshold,
    Rcpp::Named("sort_tbb_threshold") = (int)config.sort_tbb_threshold,
    Rcpp::Named("pdq_median_threshold") = (int)config.pdq_median_threshold,
    Rcpp::Named("pdq_robscale_threshold") = (int)config.pdq_robscale_threshold,
    Rcpp::Named("pdq_lowmedian_threshold") = (int)config.pdq_lowmedian_threshold,
    Rcpp::Named("pdq_qn_final_threshold") = (int)config.pdq_qn_final_threshold,
    Rcpp::Named("grain_size") = (int)config.grain_size,
    Rcpp::Named("l2_cache_size") = (int)config.hw.l2_cache_size,
    Rcpp::Named("l2_per_core") = (int)config.hw.l2_per_core,
    Rcpp::Named("num_logical_cores") = (int)config.hw.num_logical_cores,
    Rcpp::Named("num_physical_cores") = (int)config.hw.num_physical_cores,
    Rcpp::Named("has_tuned_sort_thresholds") = false,
#ifdef ROBSCALE_HAS_AVX2_DISPATCH
    Rcpp::Named("has_avx2") =
        (config.hw.simd_level >= robscale::qnsn::SIMDLevel::AVX2),
#else
    Rcpp::Named("has_avx2") = false,
#endif
#if defined(_OPENMP) || defined(ROBSCALE_HAS_OMP_SIMD)
    Rcpp::Named("has_omp_simd") = true
#else
    Rcpp::Named("has_omp_simd") = false
#endif
  );
}

// [[Rcpp::export]]
double C_get_qn_factor(int n) {
  return robscale::qnsn::get_qn_factor(static_cast<size_t>(n));
}
