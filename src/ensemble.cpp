#include "robscale_config.h"
#include "estimators_internal.h"
#include <Rcpp.h>
#include <cstring>
#include <cmath>
#include <memory>
#include <algorithm>
#include <vector>
#include <chrono>

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
// Sort resample once; use sorted-aware estimators to avoid redundant sorts.
// resample, work1, work2 are per-task scratch buffers (each n doubles).
static void ensemble_one_replicate(
    const double* xp, int n, int r,
    double* boot_row,
    double* resample, double* work1, double* work2) {

  XorShift32 rng(static_cast<uint32_t>(r + 12345));
  for (int i = 0; i < n; ++i) {
    resample[i] = xp[rng.next() % n];
  }

  // Sort resample once — all estimators below exploit sorted order
  if (n <= 16) {
    robscale::small_sort(resample, n);
  } else {
    std::sort(resample, resample + n);
  }

  // 0: sd_c4 (read-only, order-agnostic)
  boot_row[0] = robscale::internal::sd_c4(resample, n);

  // 1: gmd — inline sorted kernel, O(n) weighted sum
  {
    double sum = 0.0;
#if defined(_OPENMP) || defined(ROBSCALE_HAS_OMP_SIMD)
    #pragma omp simd reduction(+:sum)
#endif
    for (int i = 0; i < n; ++i)
      sum += (2.0 * (i + 1) - n - 1.0) * resample[i];
    boot_row[1] = (n < 2) ? 0.0
      : robscale::GMD_CONSISTENCY * 2.0 * sum / (static_cast<double>(n) * (n - 1));
  }

  // 2: mad — O(1) median from sorted data, then deviations
  {
    double med = robscale::median_sorted(resample, static_cast<size_t>(n));
    for (int i = 0; i < n; ++i) work2[i] = std::abs(resample[i] - med);
    boot_row[2] = (n < 2) ? 0.0
      : robscale::MAD_CONSISTENCY * robscale::median_select(work2, static_cast<size_t>(n));
  }

  // 3: iqr — direct index reads with interpolation
  {
    double h1 = (n - 1.0) * 0.25;
    int lo1 = static_cast<int>(h1);
    double frac1 = h1 - lo1;
    double q1 = resample[lo1];
    if (frac1 > 0.0 && lo1 + 1 < n)
      q1 += frac1 * (resample[lo1 + 1] - q1);

    double h3 = (n - 1.0) * 0.75;
    int lo3 = static_cast<int>(h3);
    double frac3 = h3 - lo3;
    double q3 = resample[lo3];
    if (frac3 > 0.0 && lo3 + 1 < n)
      q3 += frac3 * (resample[lo3 + 1] - q3);

    boot_row[3] = (n < 2) ? 0.0 : (q3 - q1) * robscale::IQR_CONSISTENCY;
  }

  // 4: sn — sorted variant, skip redundant copy+sort
  boot_row[4] = robscale::internal::sn_sorted(resample, n);

  // 5: qn — sorted variant, skip redundant copy+sort
  boot_row[5] = robscale::internal::qn_sorted(resample, n);

  // 6: robScale — O(1) median from sorted data, then Newton-Raphson
  {
    if (n < 4) {
      boot_row[6] = 0.0;
    } else {
      double t = robscale::median_sorted(resample, static_cast<size_t>(n));
      for (int i = 0; i < n; ++i) work2[i] = std::abs(resample[i] - t);
      double s_init = robscale::MAD_CONSISTENCY * robscale::median_select(work2, static_cast<size_t>(n));
      if (s_init <= 1e-4) {
        boot_row[6] = robscale::adm_core(resample, n, t, robscale::ADM_CONSISTENCY);
      } else {
        boot_row[6] = rob_scale_compute(resample, static_cast<size_t>(n),
                                         t, s_init, 80, 1.4901161e-8, work1);
      }
    }
  }
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
    // Each task allocates its own workspace: resample + 2 work buffers = 3n doubles.
    // Captures only raw pointers (xp, boot_results) — NO Rcpp objects.
    tbb::parallel_for(
      tbb::blocked_range<int>(0, n_boot),
      [xp, n, boot_results](const tbb::blocked_range<int>& range) {
        std::unique_ptr<double[]> ws(new double[3 * static_cast<size_t>(n)]);
        double* resample = ws.get();
        double* work1 = resample + n;
        double* work2 = work1 + n;

        for (int r = range.begin(); r < range.end(); ++r) {
          double* row = boot_results + static_cast<size_t>(r) * N_ESTIMATORS;
          ensemble_one_replicate(xp, n, r, row, resample, work1, work2);
        }
      }
    );
  } else
#endif
  {
    // --- Sequential path ---
    std::unique_ptr<double[]> ws(new double[3 * static_cast<size_t>(n)]);
    double* resample = ws.get();
    double* work1 = resample + n;
    double* work2 = work1 + n;

    for (int r = 0; r < n_boot; ++r) {
      double* row = boot_results + static_cast<size_t>(r) * N_ESTIMATORS;
      ensemble_one_replicate(xp, n, r, row, resample, work1, work2);
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
  std::unique_ptr<double[]> final_ws(new double[3 * static_cast<size_t>(n)]);
  double* work1 = final_ws.get();
  double* work2 = work1 + n;
  double* work_gmd = work2 + n;

  double estimates[N_ESTIMATORS];
  estimates[0] = robscale::internal::sd_c4(xp, n);

  std::memcpy(work_gmd, xp, n * sizeof(double));
  estimates[1] = robscale::internal::gmd(work_gmd, n);

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

// Helper: compute all 7 estimators on arbitrary (non-const-safe) data.
// work1, work2 must each be >= n doubles.  x is NOT modified.
static void compute_all_estimators(const double* x, int n, double* results,
                                   double* work1, double* work2) {
  results[0] = robscale::internal::sd_c4(x, n);

  std::memcpy(work1, x, n * sizeof(double));
  results[1] = robscale::internal::gmd(work1, n);

  results[2] = robscale::internal::mad_from_data(x, work1, work2, n);
  results[3] = robscale::internal::iqr(x, work1, work2, n);
  results[4] = robscale::internal::sn(x, n);
  results[5] = robscale::internal::qn(x, n);
  results[6] = robscale::internal::rob_scale(x, work1, work2, n);
}

// [[Rcpp::export]]
Rcpp::List cpp_scale_ensemble_ci(Rcpp::NumericVector x, int n_boot,
                                 double level, int method_code) {
  int n = x.size();
  if (n < 2) {
    return Rcpp::List::create(
      Rcpp::Named("estimate")    = NA_REAL,
      Rcpp::Named("ci_lower")    = NA_REAL,
      Rcpp::Named("ci_upper")    = NA_REAL,
      Rcpp::Named("level")       = level,
      Rcpp::Named("method_code") = method_code,
      Rcpp::Named("estimates")   = Rcpp::NumericVector(N_ESTIMATORS, NA_REAL),
      Rcpp::Named("weights")     = Rcpp::NumericVector(N_ESTIMATORS, NA_REAL),
      Rcpp::Named("boot_lowers") = Rcpp::NumericVector(N_ESTIMATORS, NA_REAL),
      Rcpp::Named("boot_uppers") = Rcpp::NumericVector(N_ESTIMATORS, NA_REAL),
      Rcpp::Named("z0")          = Rcpp::NumericVector(N_ESTIMATORS + 1, NA_REAL),
      Rcpp::Named("acc")         = Rcpp::NumericVector(N_ESTIMATORS + 1, NA_REAL),
      Rcpp::Named("boot_sds")    = Rcpp::NumericVector(N_ESTIMATORS + 1, NA_REAL)
    );
  }

  const double* xp = x.begin();

  // Reduce n_boot for parametric tier
  int actual_nboot = (method_code == 2) ? std::min(n_boot, 50) : n_boot;
  if (actual_nboot < 2) actual_nboot = 2;

  // --- Bootstrap loop (identical to cpp_scale_ensemble) ---
  std::unique_ptr<double[]> boot_mem(
    new double[static_cast<size_t>(actual_nboot) * N_ESTIMATORS]);
  double* boot_results = boot_mem.get();

  int64_t work_size = static_cast<int64_t>(n) * actual_nboot;

#ifdef USE_DIRECT_TBB
  if (work_size >= ENSEMBLE_PARALLEL_THRESHOLD) {
    tbb::parallel_for(
      tbb::blocked_range<int>(0, actual_nboot),
      [xp, n, boot_results](const tbb::blocked_range<int>& range) {
        std::unique_ptr<double[]> ws(new double[3 * static_cast<size_t>(n)]);
        double* resample = ws.get();
        double* work1 = resample + n;
        double* work2 = work1 + n;
        for (int r = range.begin(); r < range.end(); ++r) {
          double* row = boot_results + static_cast<size_t>(r) * N_ESTIMATORS;
          ensemble_one_replicate(xp, n, r, row, resample, work1, work2);
        }
      }
    );
  } else
#endif
  {
    std::unique_ptr<double[]> ws(new double[3 * static_cast<size_t>(n)]);
    double* resample = ws.get();
    double* work1 = resample + n;
    double* work2 = work1 + n;
    for (int r = 0; r < actual_nboot; ++r) {
      double* row = boot_results + static_cast<size_t>(r) * N_ESTIMATORS;
      ensemble_one_replicate(xp, n, r, row, resample, work1, work2);
    }
  }

  // --- Means and variances ---
  double vars[N_ESTIMATORS], means[N_ESTIMATORS];
  for (int j = 0; j < N_ESTIMATORS; ++j) {
    double sum = 0.0;
    int count = 0;
    for (int r = 0; r < actual_nboot; ++r) {
      double v = boot_results[static_cast<size_t>(r) * N_ESTIMATORS + j];
      if (std::isfinite(v) && v > 0.0) { sum += v; ++count; }
    }
    if (count < 2) { vars[j] = 1e30; means[j] = 0.0; continue; }
    means[j] = sum / count;
    double sq_sum = 0.0;
    for (int r = 0; r < actual_nboot; ++r) {
      double v = boot_results[static_cast<size_t>(r) * N_ESTIMATORS + j];
      if (std::isfinite(v) && v > 0.0) {
        double d = v - means[j];
        sq_sum += d * d;
      }
    }
    vars[j] = sq_sum / (count - 1.0);
  }

  // --- Inverse-variance weights ---
  double inv_vars[N_ESTIMATORS], weight_sum = 0.0;
  for (int j = 0; j < N_ESTIMATORS; ++j) {
    inv_vars[j] = 1.0 / std::max(1e-30, vars[j]);
    weight_sum += inv_vars[j];
  }
  double weights[N_ESTIMATORS];
  for (int j = 0; j < N_ESTIMATORS; ++j)
    weights[j] = inv_vars[j] / weight_sum;

  // --- Point estimates on original data ---
  double estimates[N_ESTIMATORS];
  {
    std::unique_ptr<double[]> final_ws(new double[3 * static_cast<size_t>(n)]);
    double* w1 = final_ws.get();
    double* w2 = w1 + n;
    double* w3 = w2 + n;
    estimates[0] = robscale::internal::sd_c4(xp, n);
    std::memcpy(w3, xp, n * sizeof(double));
    estimates[1] = robscale::internal::gmd(w3, n);
    estimates[2] = robscale::internal::mad_from_data(xp, w1, w2, n);
    estimates[3] = robscale::internal::iqr(xp, w1, w2, n);
    estimates[4] = robscale::internal::sn(xp, n);
    estimates[5] = robscale::internal::qn(xp, n);
    estimates[6] = robscale::internal::rob_scale(xp, w1, w2, n);
  }

  // Ensemble point estimate
  double ensemble_est = 0.0;
  for (int j = 0; j < N_ESTIMATORS; ++j) {
    if (std::isfinite(estimates[j]) && estimates[j] > 0.0)
      ensemble_est += weights[j] * estimates[j];
  }

  // --- Outputs ---
  double boot_lowers[N_ESTIMATORS], boot_uppers[N_ESTIMATORS];
  double z0_arr[N_ESTIMATORS + 1], acc_arr[N_ESTIMATORS + 1];
  double boot_sd_arr[N_ESTIMATORS + 1];
  for (int j = 0; j <= N_ESTIMATORS; ++j) {
    z0_arr[j] = NA_REAL;
    acc_arr[j] = NA_REAL;
    boot_sd_arr[j] = NA_REAL;
  }
  double ensemble_lower = NA_REAL, ensemble_upper = NA_REAL;

  // Compute ensemble bootstrap values (needed by all tiers)
  std::vector<double> ensemble_boot(actual_nboot);
  for (int r = 0; r < actual_nboot; ++r) {
    double val = 0.0;
    for (int j = 0; j < N_ESTIMATORS; ++j) {
      double v = boot_results[static_cast<size_t>(r) * N_ESTIMATORS + j];
      if (std::isfinite(v) && v > 0.0) val += weights[j] * v;
    }
    ensemble_boot[r] = val;
  }

  if (method_code == 0) {
    // ========================= BCa =========================
    // Step 1: Bias correction z0
    for (int j = 0; j < N_ESTIMATORS; ++j) {
      int below = 0;
      for (int r = 0; r < actual_nboot; ++r) {
        if (boot_results[static_cast<size_t>(r) * N_ESTIMATORS + j] < estimates[j])
          ++below;
      }
      double prop = static_cast<double>(below) / actual_nboot;
      prop = std::max(0.5 / actual_nboot,
                      std::min(1.0 - 0.5 / actual_nboot, prop));
      z0_arr[j] = R::qnorm(prop, 0.0, 1.0, 1, 0);
    }
    // z0 for ensemble
    {
      int below = 0;
      for (int r = 0; r < actual_nboot; ++r)
        if (ensemble_boot[r] < ensemble_est) ++below;
      double prop = static_cast<double>(below) / actual_nboot;
      prop = std::max(0.5 / actual_nboot,
                      std::min(1.0 - 0.5 / actual_nboot, prop));
      z0_arr[N_ESTIMATORS] = R::qnorm(prop, 0.0, 1.0, 1, 0);
    }

    // Step 2: Jackknife acceleration
    std::vector<double> jack_flat(static_cast<size_t>(n) * N_ESTIMATORS);
    {
      std::unique_ptr<double[]> loo(new double[n - 1]);
      std::unique_ptr<double[]> jw1(new double[n]);
      std::unique_ptr<double[]> jw2(new double[n]);
      for (int i = 0; i < n; ++i) {
        int k = 0;
        for (int ii = 0; ii < n; ++ii)
          if (ii != i) loo[k++] = xp[ii];
        compute_all_estimators(loo.get(), n - 1,
                               &jack_flat[static_cast<size_t>(i) * N_ESTIMATORS],
                               jw1.get(), jw2.get());
      }
    }

    // Per-estimator acceleration
    for (int j = 0; j < N_ESTIMATORS; ++j) {
      double jack_mean = 0.0;
      for (int i = 0; i < n; ++i)
        jack_mean += jack_flat[static_cast<size_t>(i) * N_ESTIMATORS + j];
      jack_mean /= n;
      double sum2 = 0.0, sum3 = 0.0;
      for (int i = 0; i < n; ++i) {
        double L = jack_mean -
                   jack_flat[static_cast<size_t>(i) * N_ESTIMATORS + j];
        double L2 = L * L;
        sum2 += L2;
        sum3 += L2 * L;
      }
      acc_arr[j] = (sum2 > 0.0) ? sum3 / (6.0 * std::pow(sum2, 1.5)) : 0.0;
    }
    // Ensemble acceleration
    {
      std::vector<double> ens_jack(n);
      for (int i = 0; i < n; ++i) {
        double val = 0.0;
        for (int j = 0; j < N_ESTIMATORS; ++j) {
          double v = jack_flat[static_cast<size_t>(i) * N_ESTIMATORS + j];
          if (std::isfinite(v) && v > 0.0) val += weights[j] * v;
        }
        ens_jack[i] = val;
      }
      double jack_mean = 0.0;
      for (int i = 0; i < n; ++i) jack_mean += ens_jack[i];
      jack_mean /= n;
      double sum2 = 0.0, sum3 = 0.0;
      for (int i = 0; i < n; ++i) {
        double L = jack_mean - ens_jack[i];
        double L2 = L * L;
        sum2 += L2;
        sum3 += L2 * L;
      }
      acc_arr[N_ESTIMATORS] =
        (sum2 > 0.0) ? sum3 / (6.0 * std::pow(sum2, 1.5)) : 0.0;
    }

    // Step 3: Adjusted percentiles
    double alpha = 1.0 - level;
    double z_lo = R::qnorm(alpha / 2.0, 0.0, 1.0, 1, 0);
    double z_hi = R::qnorm(1.0 - alpha / 2.0, 0.0, 1.0, 1, 0);

    std::vector<double> col(actual_nboot);
    for (int j = 0; j < N_ESTIMATORS; ++j) {
      for (int r = 0; r < actual_nboot; ++r)
        col[r] = boot_results[static_cast<size_t>(r) * N_ESTIMATORS + j];
      std::sort(col.begin(), col.end());

      double numer_lo = z0_arr[j] + z_lo;
      double numer_hi = z0_arr[j] + z_hi;
      double a1 = R::pnorm(
        z0_arr[j] + numer_lo / (1.0 - acc_arr[j] * numer_lo),
        0.0, 1.0, 1, 0);
      double a2 = R::pnorm(
        z0_arr[j] + numer_hi / (1.0 - acc_arr[j] * numer_hi),
        0.0, 1.0, 1, 0);
      a1 = std::max(0.5 / actual_nboot,
                     std::min(1.0 - 0.5 / actual_nboot, a1));
      a2 = std::max(0.5 / actual_nboot,
                     std::min(1.0 - 0.5 / actual_nboot, a2));
      int idx1 = std::max(0, std::min(actual_nboot - 1,
                   static_cast<int>(std::floor(a1 * actual_nboot))));
      int idx2 = std::max(0, std::min(actual_nboot - 1,
                   static_cast<int>(std::floor(a2 * actual_nboot))));
      boot_lowers[j] = col[idx1];
      boot_uppers[j] = col[idx2];
    }

    // Ensemble BCa CI
    {
      std::sort(ensemble_boot.begin(), ensemble_boot.end());
      double numer_lo = z0_arr[N_ESTIMATORS] + z_lo;
      double numer_hi = z0_arr[N_ESTIMATORS] + z_hi;
      double a1 = R::pnorm(
        z0_arr[N_ESTIMATORS] +
          numer_lo / (1.0 - acc_arr[N_ESTIMATORS] * numer_lo),
        0.0, 1.0, 1, 0);
      double a2 = R::pnorm(
        z0_arr[N_ESTIMATORS] +
          numer_hi / (1.0 - acc_arr[N_ESTIMATORS] * numer_hi),
        0.0, 1.0, 1, 0);
      a1 = std::max(0.5 / actual_nboot,
                     std::min(1.0 - 0.5 / actual_nboot, a1));
      a2 = std::max(0.5 / actual_nboot,
                     std::min(1.0 - 0.5 / actual_nboot, a2));
      int idx1 = std::max(0, std::min(actual_nboot - 1,
                   static_cast<int>(std::floor(a1 * actual_nboot))));
      int idx2 = std::max(0, std::min(actual_nboot - 1,
                   static_cast<int>(std::floor(a2 * actual_nboot))));
      ensemble_lower = ensemble_boot[idx1];
      ensemble_upper = ensemble_boot[idx2];
    }

  } else if (method_code == 1) {
    // ====================== Percentile ======================
    double alpha = 1.0 - level;
    int lo_idx = std::max(0,
      static_cast<int>(std::floor(alpha / 2.0 * actual_nboot)));
    int hi_idx = std::min(actual_nboot - 1,
      static_cast<int>(std::floor((1.0 - alpha / 2.0) * actual_nboot)));

    std::vector<double> col(actual_nboot);
    for (int j = 0; j < N_ESTIMATORS; ++j) {
      for (int r = 0; r < actual_nboot; ++r)
        col[r] = boot_results[static_cast<size_t>(r) * N_ESTIMATORS + j];
      std::sort(col.begin(), col.end());
      boot_lowers[j] = col[lo_idx];
      boot_uppers[j] = col[hi_idx];
    }

    std::sort(ensemble_boot.begin(), ensemble_boot.end());
    ensemble_lower = ensemble_boot[lo_idx];
    ensemble_upper = ensemble_boot[hi_idx];

  } else {
    // ====================== Parametric ======================
    double z = R::qnorm(1.0 - (1.0 - level) / 2.0, 0.0, 1.0, 1, 0);

    for (int j = 0; j < N_ESTIMATORS; ++j) {
      double sd_j = std::sqrt(std::max(0.0, vars[j]));
      boot_sd_arr[j] = sd_j;
      boot_lowers[j] = estimates[j] - z * sd_j;
      boot_uppers[j] = estimates[j] + z * sd_j;
    }

    double ens_mean = 0.0;
    for (int r = 0; r < actual_nboot; ++r) ens_mean += ensemble_boot[r];
    ens_mean /= actual_nboot;
    double ens_sq = 0.0;
    for (int r = 0; r < actual_nboot; ++r) {
      double d = ensemble_boot[r] - ens_mean;
      ens_sq += d * d;
    }
    double ens_sd = std::sqrt(ens_sq / (actual_nboot - 1.0));
    boot_sd_arr[N_ESTIMATORS] = ens_sd;
    ensemble_lower = ensemble_est - z * ens_sd;
    ensemble_upper = ensemble_est + z * ens_sd;
  }

  // --- Build return list ---
  Rcpp::NumericVector est_vec(N_ESTIMATORS), wt_vec(N_ESTIMATORS);
  Rcpp::NumericVector bl_vec(N_ESTIMATORS), bu_vec(N_ESTIMATORS);
  Rcpp::NumericVector z0_vec(N_ESTIMATORS + 1), acc_vec(N_ESTIMATORS + 1);
  Rcpp::NumericVector sd_vec(N_ESTIMATORS + 1);

  for (int j = 0; j < N_ESTIMATORS; ++j) {
    est_vec[j] = estimates[j];
    wt_vec[j]  = weights[j];
    bl_vec[j]  = boot_lowers[j];
    bu_vec[j]  = boot_uppers[j];
  }
  for (int j = 0; j <= N_ESTIMATORS; ++j) {
    z0_vec[j]  = z0_arr[j];
    acc_vec[j] = acc_arr[j];
    sd_vec[j]  = boot_sd_arr[j];
  }

  return Rcpp::List::create(
    Rcpp::Named("estimate")    = ensemble_est,
    Rcpp::Named("ci_lower")    = ensemble_lower,
    Rcpp::Named("ci_upper")    = ensemble_upper,
    Rcpp::Named("level")       = level,
    Rcpp::Named("method_code") = method_code,
    Rcpp::Named("estimates")   = est_vec,
    Rcpp::Named("weights")     = wt_vec,
    Rcpp::Named("boot_lowers") = bl_vec,
    Rcpp::Named("boot_uppers") = bu_vec,
    Rcpp::Named("z0")          = z0_vec,
    Rcpp::Named("acc")         = acc_vec,
    Rcpp::Named("boot_sds")    = sd_vec
  );
}

// ============================================================================
// Benchmark harness — temporary, guarded by ROBSCALE_BENCH_INTERNALS
// ============================================================================
#ifdef ROBSCALE_BENCH_INTERNALS

// [[Rcpp::export]]
double C_sn_sorted_test(Rcpp::NumericVector x) {
  return robscale::qnsn::C_sn_impl_sorted<double>(x.begin(), static_cast<size_t>(x.size()));
}

// [[Rcpp::export]]
double C_qn_sorted_test(Rcpp::NumericVector x) {
  return robscale::qnsn::C_qn_impl_sorted<double>(x.begin(), static_cast<size_t>(x.size()));
}

// [[Rcpp::export]]
Rcpp::DataFrame bench_ensemble_internals(Rcpp::NumericVector x, int n_boot, int reps) {
  using Clock = std::chrono::high_resolution_clock;
  int n = x.size();
  const double* xp = x.begin();

  Rcpp::IntegerVector rep_col(reps);
  Rcpp::NumericVector total_us_col(reps);

  size_t ws_size = 3 * static_cast<size_t>(n);
  std::unique_ptr<double[]> ws(new double[ws_size]);
  std::unique_ptr<double[]> boot_mem(new double[static_cast<size_t>(n_boot) * N_ESTIMATORS]);

  double* resample = ws.get();
  double* work1 = resample + n;
  double* work2 = work1 + n;
  double* boot_results = boot_mem.get();

  for (int rep = 0; rep < reps; ++rep) {
    auto t0 = Clock::now();
    for (int r = 0; r < n_boot; ++r) {
      double* row = boot_results + static_cast<size_t>(r) * N_ESTIMATORS;
      ensemble_one_replicate(xp, n, r, row, resample, work1, work2);
    }
    auto t1 = Clock::now();

    rep_col[rep] = rep + 1;
    total_us_col[rep] = std::chrono::duration<double, std::micro>(t1 - t0).count();
  }

  return Rcpp::DataFrame::create(
    Rcpp::Named("rep") = rep_col,
    Rcpp::Named("total_us") = total_us_col,
    Rcpp::Named("stringsAsFactors") = false
  );
}

// [[Rcpp::export]]
double bench_sort_cost_us(Rcpp::NumericVector x, int reps) {
  using Clock = std::chrono::high_resolution_clock;
  int n = x.size();
  const double* xp = x.begin();

  std::unique_ptr<double[]> buf(new double[n]);
  std::vector<double> timings(reps);

  for (int rep = 0; rep < reps; ++rep) {
    std::memcpy(buf.get(), xp, n * sizeof(double));
    auto t0 = Clock::now();
    if (n <= 16) {
      robscale::small_sort(buf.get(), n);
    } else {
      std::sort(buf.get(), buf.get() + n);
    }
    auto t1 = Clock::now();
    timings[rep] = std::chrono::duration<double, std::micro>(t1 - t0).count();
  }

  std::sort(timings.begin(), timings.end());
  return timings[reps / 2]; // median
}

// [[Rcpp::export]]
Rcpp::NumericVector bench_estimator_breakdown(Rcpp::NumericVector x, int reps) {
  using Clock = std::chrono::high_resolution_clock;
  int n = x.size();
  const double* xp = x.begin();

  // 7 estimators, reps timings each
  std::vector<std::vector<double>> timings(N_ESTIMATORS, std::vector<double>(reps));

  // Workspace: 5n
  std::unique_ptr<double[]> ws(new double[5 * static_cast<size_t>(n)]);
  double* resample = ws.get();
  double* work1 = resample + n;
  double* work2 = work1 + n;
  double* work3 = work2 + n;
  double* work4 = work3 + n;

  double boot_row[N_ESTIMATORS];

  for (int rep = 0; rep < reps; ++rep) {
    // Generate a fresh resample for each rep
    XorShift32 rng(static_cast<uint32_t>(rep + 99999));
    for (int i = 0; i < n; ++i) {
      resample[i] = xp[rng.next() % n];
    }

    // 0: sd_c4
    {
      auto t0 = Clock::now();
      boot_row[0] = robscale::internal::sd_c4(resample, n);
      auto t1 = Clock::now();
      timings[0][rep] = std::chrono::duration<double, std::micro>(t1 - t0).count();
    }

    // 1: gmd (needs sorted copy)
    {
      std::memcpy(work4, resample, n * sizeof(double));
      auto t0 = Clock::now();
      boot_row[1] = robscale::internal::gmd(work4, n);
      auto t1 = Clock::now();
      timings[1][rep] = std::chrono::duration<double, std::micro>(t1 - t0).count();
    }

    // 2: mad
    {
      auto t0 = Clock::now();
      boot_row[2] = robscale::internal::mad_from_data(resample, work1, work2, n);
      auto t1 = Clock::now();
      timings[2][rep] = std::chrono::duration<double, std::micro>(t1 - t0).count();
    }

    // 3: iqr
    {
      auto t0 = Clock::now();
      boot_row[3] = robscale::internal::iqr(resample, work1, work2, n);
      auto t1 = Clock::now();
      timings[3][rep] = std::chrono::duration<double, std::micro>(t1 - t0).count();
    }

    // 4: sn
    {
      auto t0 = Clock::now();
      boot_row[4] = robscale::internal::sn(resample, n);
      auto t1 = Clock::now();
      timings[4][rep] = std::chrono::duration<double, std::micro>(t1 - t0).count();
    }

    // 5: qn
    {
      auto t0 = Clock::now();
      boot_row[5] = robscale::internal::qn(resample, n);
      auto t1 = Clock::now();
      timings[5][rep] = std::chrono::duration<double, std::micro>(t1 - t0).count();
    }

    // 6: robScale
    {
      auto t0 = Clock::now();
      boot_row[6] = robscale::internal::rob_scale(resample, work1, work2, n);
      auto t1 = Clock::now();
      timings[6][rep] = std::chrono::duration<double, std::micro>(t1 - t0).count();
    }
  }

  // Compute medians
  Rcpp::NumericVector result(N_ESTIMATORS);
  Rcpp::CharacterVector names = {"sd_c4", "gmd", "mad", "iqr", "sn", "qn", "rob_scale"};
  for (int j = 0; j < N_ESTIMATORS; ++j) {
    std::sort(timings[j].begin(), timings[j].end());
    result[j] = timings[j][reps / 2];
  }
  result.attr("names") = names;
  return result;
}

#endif // ROBSCALE_BENCH_INTERNALS
