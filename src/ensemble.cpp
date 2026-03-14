#include "robscale_config.h"
#include "estimators_internal.h"
#include <Rcpp.h>
#include <cstring>
#include <cmath>
#include <memory>
#include <algorithm>

#include <RcppParallel.h>
#ifdef USE_DIRECT_TBB
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#endif

// XorShift32 PRNG for deterministic bootstrap resampling
struct XorShift32 {
  uint32_t state;
  explicit XorShift32(uint32_t seed) : state(seed ? seed : 1) {}
  uint32_t next() {
    state ^= state << 13;
    state ^= state >> 17;
    state ^= state << 5;
    return state;
  }
};

static constexpr int N_ESTIMATORS = 7;
static constexpr int64_t ENSEMBLE_PARALLEL_THRESHOLD = 10000;

// Run all 7 estimators on a single bootstrap replicate.
// resample, work1..work4 are per-task scratch buffers (each n doubles).
static void ensemble_one_replicate(
    const double* xp, int n, int r,
    double* boot_row,
    double* resample, double* work1, double* work2,
    double* work3, double* work4) {

  XorShift32 rng(static_cast<uint32_t>(r + 12345));
  for (int i = 0; i < n; ++i) {
    resample[i] = xp[rng.next() % n];
  }

  // 0: sd_c4 (read-only)
  boot_row[0] = robscale::internal::sd_c4(resample, n);

  // 1: gmd (needs sorted copy)
  std::memcpy(work4, resample, n * sizeof(double));
  boot_row[1] = robscale::internal::gmd(work4, n);

  // 2: mad (needs buf + dev)
  boot_row[2] = robscale::internal::mad_from_data(resample, work1, work2, n);

  // 3: iqr (needs two buffers)
  boot_row[3] = robscale::internal::iqr(resample, work1, work2, n);

  // 4: sn (read-only, internally allocates)
  boot_row[4] = robscale::internal::sn(resample, n);

  // 5: qn (read-only, internally allocates)
  boot_row[5] = robscale::internal::qn(resample, n);

  // 6: robScale (needs buf + dev)
  boot_row[6] = robscale::internal::rob_scale(resample, work1, work2, n);
}

// [[Rcpp::export]]
double cpp_scale_ensemble(Rcpp::NumericVector x, int n_boot) {
  int n = x.size();
  if (n < 2) return NA_REAL;
  if (n_boot < 2) n_boot = 2;

  const double* xp = x.begin();

  // Allocate boot_results: n_boot * N_ESTIMATORS (shared, non-overlapping writes)
  std::unique_ptr<double[]> boot_mem(new double[static_cast<size_t>(n_boot) * N_ESTIMATORS]);
  double* boot_results = boot_mem.get();

  int64_t work_size = static_cast<int64_t>(n) * n_boot;

#ifdef USE_DIRECT_TBB
  if (work_size >= ENSEMBLE_PARALLEL_THRESHOLD) {
    // --- TBB parallel path ---
    // Each task allocates its own workspace: resample + 4 work buffers = 5n doubles.
    // Captures only raw pointers (xp, boot_results) — NO Rcpp objects.
    tbb::parallel_for(
      tbb::blocked_range<int>(0, n_boot),
      [xp, n, boot_results](const tbb::blocked_range<int>& range) {
        // Per-task workspace allocation (5n doubles)
        std::unique_ptr<double[]> ws(new double[5 * static_cast<size_t>(n)]);
        double* resample = ws.get();
        double* work1 = resample + n;
        double* work2 = work1 + n;
        double* work3 = work2 + n;
        double* work4 = work3 + n;

        for (int r = range.begin(); r < range.end(); ++r) {
          double* row = boot_results + static_cast<size_t>(r) * N_ESTIMATORS;
          ensemble_one_replicate(xp, n, r, row, resample, work1, work2, work3, work4);
        }
      }
    );
  } else
#endif
  {
    // --- Sequential path ---
    std::unique_ptr<double[]> ws(new double[5 * static_cast<size_t>(n)]);
    double* resample = ws.get();
    double* work1 = resample + n;
    double* work2 = work1 + n;
    double* work3 = work2 + n;
    double* work4 = work3 + n;

    for (int r = 0; r < n_boot; ++r) {
      double* row = boot_results + static_cast<size_t>(r) * N_ESTIMATORS;
      ensemble_one_replicate(xp, n, r, row, resample, work1, work2, work3, work4);
    }
  }

  // Compute inverse-variance weights
  double vars[N_ESTIMATORS];
  double means[N_ESTIMATORS];
  for (int j = 0; j < N_ESTIMATORS; ++j) {
    double sum = 0.0;
    int count = 0;
    for (int r = 0; r < n_boot; ++r) {
      double v = boot_results[static_cast<size_t>(r) * N_ESTIMATORS + j];
      if (std::isfinite(v) && v > 0.0) {
        sum += v;
        ++count;
      }
    }
    if (count < 2) {
      vars[j] = 1e30; // effectively zero weight
      means[j] = 0.0;
      continue;
    }
    means[j] = sum / count;

    double sq_sum = 0.0;
    for (int r = 0; r < n_boot; ++r) {
      double v = boot_results[static_cast<size_t>(r) * N_ESTIMATORS + j];
      if (std::isfinite(v) && v > 0.0) {
        double d = v - means[j];
        sq_sum += d * d;
      }
    }
    vars[j] = sq_sum / (count - 1.0);
  }

  double inv_vars[N_ESTIMATORS];
  double weight_sum = 0.0;
  for (int j = 0; j < N_ESTIMATORS; ++j) {
    inv_vars[j] = 1.0 / std::max(1e-30, vars[j]);
    weight_sum += inv_vars[j];
  }

  // Final estimates on original data (fresh workspace)
  std::unique_ptr<double[]> final_ws(new double[5 * static_cast<size_t>(n)]);
  double* work1 = final_ws.get();
  double* work2 = work1 + n;
  double* work3 = work2 + n;
  double* work4 = work3 + n;

  double estimates[N_ESTIMATORS];
  estimates[0] = robscale::internal::sd_c4(xp, n);

  std::memcpy(work4, xp, n * sizeof(double));
  estimates[1] = robscale::internal::gmd(work4, n);

  estimates[2] = robscale::internal::mad_from_data(xp, work1, work2, n);
  estimates[3] = robscale::internal::iqr(xp, work1, work2, n);
  estimates[4] = robscale::internal::sn(xp, n);
  estimates[5] = robscale::internal::qn(xp, n);
  estimates[6] = robscale::internal::rob_scale(xp, work1, work2, n);

  // Weighted combination
  double result = 0.0;
  for (int j = 0; j < N_ESTIMATORS; ++j) {
    if (std::isfinite(estimates[j]) && estimates[j] > 0.0) {
      result += (inv_vars[j] / weight_sum) * estimates[j];
    }
  }

  return result;
}
