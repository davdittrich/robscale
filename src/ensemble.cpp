#include "robscale_config.h"
#include "vshaped_mad.h"
#include "estimators_internal.h"
#include <Rcpp.h>
#include <cstring>
#include <cmath>
#include <memory>
#include <algorithm>
#include <vector>

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

// Remap estimator_id (API order: 0=gmd, 1=sd_c4, 2=mad, ...) to the
// compute_all_estimators results[] order (0=sd_c4, 1=gmd, 2=mad, ...).
// Only indices 0 and 1 are swapped; 2-6 are identity.
static const int kEstIdToAllIdx[N_ESTIMATORS] = {1, 0, 2, 3, 4, 5, 6};

#ifdef USE_DIRECT_TBB
static constexpr int64_t ENSEMBLE_PARALLEL_THRESHOLD = 10000;
#endif

// Compute one bootstrap replicate and write N_ESTIMATORS estimates into the
// transposed (7 × nboot) boot_results matrix.
//
// Parameters:
//   xp        — pointer to the original data (length n); never modified.
//   n         — number of observations.
//   r         — replicate index (0-based); seeds the XorShift32 PRNG as
//               (r + 12345) to give each replicate a distinct, reproducible
//               sequence.
//   base      — start of the 7×nboot matrix; estimator j for replicate r
//               is written to base[j * nboot + r].
//   nboot     — total number of bootstrap replicates (used as column stride).
//   resample  — scratch buffer of length n (overwritten with the resample).
//   work1     — scratch buffer of length n for estimators that need it.
//   work2     — scratch buffer of length n for estimators that need it.
//
// WU-LAYOUT-1: transposed layout (7×nboot instead of nboot×7) makes the
// mean/variance reduction pass sequential (j outer, r inner) and removes
// the stride-7 access pattern in EnsembleCore::run().
//
// The resample is drawn with replacement using a XorShift32 PRNG (three XOR
// shifts: <<13, >>17, <<5), then sorted once in ascending order.  All seven
// estimators (sd_c4, gmd, mad, iqr, sn, qn, adm) exploit the sorted order to
// avoid redundant sorting passes.
static void ensemble_one_replicate(
    const double* xp, int n, int r,
    double* base, int nboot,
    double* resample, double* work1, double* work2,
    float* qn_work, int32_t* qn_iweight,
    int32_t* qn_left, int32_t* qn_right) {
  // Helper: write estimator j for replicate r.
  // base[j * nboot + r] in the 7×nboot column-major matrix.
  auto write_est = [base, nboot, r](int j, double v) {
    base[static_cast<size_t>(j) * nboot + r] = v;
  };

  XorShift32 rng(static_cast<uint32_t>(r + 12345));
  const uint32_t un32 = static_cast<uint32_t>(n);
  for (int i = 0; i < n; ++i) {
    // Lemire's bias-free bounded random: (rng * n) >> 32
    resample[i] = xp[static_cast<uint32_t>((static_cast<uint64_t>(rng.next()) * un32) >> 32)];
  }

  // Sort resample once — all estimators below exploit sorted order.
  // OPT-G4: boost::float_sort (serial radix) for n > sort_boost_threshold.
  //   tbb::parallel_sort is PROHIBITED: nested TBB inside tbb::parallel_for
  //   causes oversubscription (same constraint as rob_scale_compute here).
  if (n <= 16) {
    robscale::small_sort(resample, n);
  } else {
    if (static_cast<size_t>(n) <= ROBSCALE_SORT_BOOST_THRESHOLD)
      std::sort(resample, resample + n);
    else
      boost::sort::spreadsort::float_sort(resample, resample + n);
  }

  // 0: sd_c4 (read-only, order-agnostic)
  write_est(0, robscale::internal::sd_c4(resample, n));

  // 1: gmd — WU-GMD-1: shared kernel in robust_core.h (AVX2 FMA or scalar).
  // OPT-G5: scale precomputed once outside the accumulation loop.
  {
    const double gmd_scale = (n < 2) ? 0.0
      : robscale::GMD_CONSISTENCY * 2.0
        / (static_cast<double>(n) * (n - 1));
    write_est(1, robscale::gmd_weighted_sum(resample, n, gmd_scale));
  }

  // 2: mad — V-shaped O(log n) deviation median (OPT-M7)
  // work2 is reused as scratch; estimator 6 overwrites it independently below.
  {
    double med = robscale::median_sorted(resample, static_cast<size_t>(n));
    write_est(2, (n < 2) ? 0.0
      : robscale::MAD_CONSISTENCY * robscale::vshaped_mad(resample, n, med, work2));
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

    write_est(3, (n < 2) ? 0.0 : (q3 - q1) * robscale::IQR_CONSISTENCY);
  }

  // 4: sn — sorted variant, skip redundant copy+sort
  // OPT-S7: pass work1 as workspace to avoid heap allocation for n > sn_stack_threshold.
  // work1 is free at this point: rob_scale_compute (estimator 6) writes it below.
  write_est(4, robscale::internal::sn_sorted(resample, n, work1));

  // 5: qn — sorted variant, skip redundant copy+sort
  // Qn workspace uses correctly-typed arrays passed by caller (no aliasing UB).
  // Guard: only for n <= QN_WS_THRESHOLD (2048); beyond that pass nullptr.
  {
    using WS = robscale::qnsn::QnWorkspace;
    static constexpr size_t QN_WS_THRESHOLD = 2048;
    const size_t un = static_cast<size_t>(n);
    WS qn_ws{};
    WS* qn_ws_ptr = nullptr;
    if (un <= QN_WS_THRESHOLD && qn_work != nullptr) {
      qn_ws.work    = qn_work;
      qn_ws.iweight = qn_iweight;
      qn_ws.left    = qn_left;
      qn_ws.right   = qn_right;
      qn_ws_ptr = &qn_ws;
    }
    write_est(5, robscale::internal::qn_sorted(resample, n, qn_ws_ptr));
  }

  // 6: robScale — OPT-9: rob_scale_sorted combines O(1) median + O(log n) MAD
  //    + Newton-Raphson in one call; resample is already sorted.
  write_est(6, robscale::internal::rob_scale_sorted(
      resample, static_cast<size_t>(n), work1));
}

// Helper: compute all 7 estimators on arbitrary (non-const-safe) data.
// work1, work2 must each be >= n doubles.  x is NOT modified.
static void compute_all_estimators(const double* x, int n, double* results,
                                   double* work1, double* work2) {
  results[0] = robscale::internal::sd_c4(x, n);

  std::memcpy(work1, x, n * sizeof(double));
  results[1] = robscale::internal::gmd(work1, n);

  results[2] = robscale::internal::mad_from_data(x, work1, n);
  results[3] = robscale::internal::iqr(x, work1, n);
  results[4] = robscale::internal::sn(x, n);
  results[5] = robscale::internal::qn(x, n);
  results[6] = robscale::internal::rob_scale(x, work1, n);
}

// Shared bootstrap + statistics pipeline used by both exported functions.
// After run(), callers access boot_results, weights, estimates, ensemble_est.
struct EnsembleCore {
  std::unique_ptr<double[]> boot_mem;
  double* boot_results = nullptr;
  int nboot = 0;

  double vars[N_ESTIMATORS];
  double means[N_ESTIMATORS];
  double weights[N_ESTIMATORS];  // normalized inverse-variance weights
  double estimates[N_ESTIMATORS];
  double ensemble_est = 0.0;

  void run(const double* xp, int n, int n_boot) {
    nboot = n_boot;
    // WU-LAYOUT-1: allocate 7×nboot (estimator-major) instead of nboot×7.
    // Mean/variance pass reads boot_results[j*nboot+r] (sequential per j),
    // eliminating the stride-7 access pattern of the old nboot×7 layout.
    boot_mem.reset(new double[static_cast<size_t>(N_ESTIMATORS) * nboot]);
    boot_results = boot_mem.get();

    // --- Bootstrap loop ---
    double* br = boot_results;  // raw pointer safe to capture in lambda

#ifdef USE_DIRECT_TBB
    if (static_cast<int64_t>(n) * nboot >= ENSEMBLE_PARALLEL_THRESHOLD) {
      tbb::parallel_for(
        tbb::blocked_range<int>(0, nboot),
        [xp, n, br, n_boot](const tbb::blocked_range<int>& range) {
          const size_t un = static_cast<size_t>(n);
          std::unique_ptr<double[]> ws(new double[3 * un]);
          double* resample = ws.get();
          double* work1 = resample + n;
          double* work2 = work1 + n;
          // Per-thread typed Qn workspace (no aliasing UB)
          static constexpr size_t QN_WS_THRESHOLD = 2048;
          std::unique_ptr<float[]>   qn_w(un <= QN_WS_THRESHOLD ? new float[un] : nullptr);
          std::unique_ptr<int32_t[]> qn_iw(un <= QN_WS_THRESHOLD ? new int32_t[un] : nullptr);
          std::unique_ptr<int32_t[]> qn_l(un <= QN_WS_THRESHOLD ? new int32_t[un] : nullptr);
          std::unique_ptr<int32_t[]> qn_r(un <= QN_WS_THRESHOLD ? new int32_t[un] : nullptr);
          for (int r = range.begin(); r < range.end(); ++r)
            ensemble_one_replicate(xp, n, r, br, n_boot, resample, work1, work2,
                                   qn_w.get(), qn_iw.get(), qn_l.get(), qn_r.get());
        }
      );
    } else
#endif
    {
      const size_t un = static_cast<size_t>(n);
      std::unique_ptr<double[]> ws(new double[3 * un]);
      double* resample = ws.get();
      double* work1 = resample + n;
      double* work2 = work1 + n;
      // Typed Qn workspace (no aliasing UB)
      static constexpr size_t QN_WS_THRESHOLD = 2048;
      std::unique_ptr<float[]>   qn_w(un <= QN_WS_THRESHOLD ? new float[un] : nullptr);
      std::unique_ptr<int32_t[]> qn_iw(un <= QN_WS_THRESHOLD ? new int32_t[un] : nullptr);
      std::unique_ptr<int32_t[]> qn_l(un <= QN_WS_THRESHOLD ? new int32_t[un] : nullptr);
      std::unique_ptr<int32_t[]> qn_r(un <= QN_WS_THRESHOLD ? new int32_t[un] : nullptr);
      for (int r = 0; r < nboot; ++r)
        ensemble_one_replicate(xp, n, r, br, nboot, resample, work1, work2,
                               qn_w.get(), qn_iw.get(), qn_l.get(), qn_r.get());
    }

    // --- Mean and variance per estimator ---
    // WU-LAYOUT-1: boot_results[j*nboot+r] — sequential reads for fixed j.
    for (int j = 0; j < N_ESTIMATORS; ++j) {
      const double* col = boot_results + static_cast<size_t>(j) * nboot;
      double sum = 0.0;
      int count = 0;
      for (int r = 0; r < nboot; ++r) {
        double v = col[r];
        if (std::isfinite(v) && v > 0.0) { sum += v; ++count; }
      }
      if (count < 2) { vars[j] = 1e30; means[j] = 0.0; continue; }
      means[j] = sum / count;
      double sq_sum = 0.0;
      for (int r = 0; r < nboot; ++r) {
        double v = col[r];
        if (std::isfinite(v) && v > 0.0) {
          double d = v - means[j];
          sq_sum += d * d;
        }
      }
      vars[j] = sq_sum / (count - 1.0);
    }

    // --- Normalized inverse-variance weights ---
    double weight_sum = 0.0;
    for (int j = 0; j < N_ESTIMATORS; ++j)
      weight_sum += 1.0 / std::max(1e-30, vars[j]);
    for (int j = 0; j < N_ESTIMATORS; ++j)
      weights[j] = (1.0 / std::max(1e-30, vars[j])) / weight_sum;

    // --- Point estimates on original data ---
    {
      std::unique_ptr<double[]> ws(new double[2 * static_cast<size_t>(n)]);
      double* w1 = ws.get();
      double* w2 = w1 + n;
      compute_all_estimators(xp, n, estimates, w1, w2);
    }

    // --- Weighted ensemble estimate ---
    ensemble_est = 0.0;
    for (int j = 0; j < N_ESTIMATORS; ++j) {
      if (std::isfinite(estimates[j]) && estimates[j] > 0.0)
        ensemble_est += weights[j] * estimates[j];
    }
  }
};

// Bootstrap CI for a single named estimator.
//
// Runs n_boot resamples, computing only the requested estimator per replicate
// (avoids the 7x overhead of cpp_scale_ensemble_ci for individual methods).
//
// Parameters:
//   x            — data vector (validated: no NAs, n >= 2)
//   est          — pre-computed point estimate on original data
//   estimator_id — 0=gmd, 1=sd, 2=mad, 3=iqr, 4=sn, 5=qn, 6=robScale
//   n_boot       — bootstrap replicates
//   level        — confidence level (e.g. 0.95)
//   method_code  — BCa=0, percentile=1, parametric=2
//                  BCa falls back to percentile when z0 or adjusted quantiles
//                  are non-finite.
//
// Returns List with ci_lower and ci_upper.
// [[Rcpp::export]]
Rcpp::List cpp_single_estimator_ci_bounds(
    Rcpp::NumericVector x, double est, int estimator_id,
    int n_boot, double level, int method_code) {

  using Rcpp::Named;
  int n = x.size();

  if (n < 2 || !std::isfinite(est)) {
    return Rcpp::List::create(
      Named("ci_lower") = NA_REAL,
      Named("ci_upper") = NA_REAL
    );
  }

  const double* xp = x.begin();
  if (n_boot < 2) n_boot = 2;

  // --- Bootstrap loop ---
  std::vector<double> boot_vals(static_cast<size_t>(n_boot));
  const size_t un_ci = static_cast<size_t>(n);
  static constexpr size_t QN_CI_THRESHOLD = 2048;
  // Typed Qn workspace — allocated once before the loop (no aliasing UB).
  std::unique_ptr<float[]>   qn_w_ci(un_ci <= QN_CI_THRESHOLD ? new float[un_ci] : nullptr);
  std::unique_ptr<int32_t[]> qn_iw_ci(un_ci <= QN_CI_THRESHOLD ? new int32_t[un_ci] : nullptr);
  std::unique_ptr<int32_t[]> qn_l_ci(un_ci <= QN_CI_THRESHOLD ? new int32_t[un_ci] : nullptr);
  std::unique_ptr<int32_t[]> qn_r_ci(un_ci <= QN_CI_THRESHOLD ? new int32_t[un_ci] : nullptr);
  {
    std::unique_ptr<double[]> ws(new double[3 * un_ci]);
    double* resample = ws.get();
    double* work1    = resample + n;
    double* work2    = work1   + n;

    for (int r = 0; r < n_boot; ++r) {
      XorShift32 rng(static_cast<uint32_t>(r + 12345));
      const uint32_t un32_ci = static_cast<uint32_t>(n);
      for (int i = 0; i < n; ++i)
        resample[i] = xp[static_cast<uint32_t>((static_cast<uint64_t>(rng.next()) * un32_ci) >> 32)];

      if (n <= 16) {
        robscale::small_sort(resample, n);
      } else if (static_cast<size_t>(n) <= ROBSCALE_SORT_BOOST_THRESHOLD) {
        std::sort(resample, resample + n);
      } else {
        boost::sort::spreadsort::float_sort(resample, resample + n);
      }

      double val = 0.0;
      switch (estimator_id) {
        case 0: {  // gmd
          // WU-GMD-1: shared kernel in robust_core.h (AVX2 FMA or scalar).
          const double gmd_scale = robscale::GMD_CONSISTENCY * 2.0
            / (static_cast<double>(n) * (n - 1));
          val = robscale::gmd_weighted_sum(resample, n, gmd_scale);
          break;
        }
        case 1:  // sd_c4
          val = robscale::internal::sd_c4(resample, n);
          break;
        case 2: {  // mad_scaled
          double med = robscale::median_sorted(resample, static_cast<size_t>(n));
          val = (n < 2) ? 0.0
              : robscale::MAD_CONSISTENCY
                * robscale::vshaped_mad(resample, n, med, work2);
          break;
        }
        case 3: {  // iqr_scaled
          double h1 = (n - 1.0) * 0.25;
          int    lo1 = static_cast<int>(h1);
          double q1  = resample[lo1];
          double f1  = h1 - lo1;
          if (f1 > 0.0 && lo1 + 1 < n) q1 += f1 * (resample[lo1 + 1] - q1);
          double h3 = (n - 1.0) * 0.75;
          int    lo3 = static_cast<int>(h3);
          double q3  = resample[lo3];
          double f3  = h3 - lo3;
          if (f3 > 0.0 && lo3 + 1 < n) q3 += f3 * (resample[lo3 + 1] - q3);
          val = (q3 - q1) * robscale::IQR_CONSISTENCY;
          break;
        }
        case 4:  // sn
          val = robscale::internal::sn_sorted(resample, n, work1);
          break;
        case 5: {  // qn
          using WS = robscale::qnsn::QnWorkspace;
          WS  qn_ws{};
          WS* qn_ws_ptr = nullptr;
          if (un_ci <= QN_CI_THRESHOLD) {
            qn_ws.work    = qn_w_ci.get();
            qn_ws.iweight = qn_iw_ci.get();
            qn_ws.left    = qn_l_ci.get();
            qn_ws.right   = qn_r_ci.get();
            qn_ws_ptr = &qn_ws;
          }
          val = robscale::internal::qn_sorted(resample, n, qn_ws_ptr);
          break;
        }
        default:  // robScale (6)
          val = robscale::internal::rob_scale_sorted(
            resample, static_cast<size_t>(n), work1);
          break;
      }
      boot_vals[r] = val;
    }
  }

  // --- CI computation ---
  double ci_lower = NA_REAL, ci_upper = NA_REAL;

  if (method_code == 0) {
    // ========================= BCa =========================
    // Step 1: bias-correction z0
    int below = 0;
    for (int r = 0; r < n_boot; ++r)
      if (boot_vals[r] < est) ++below;
    double prop = static_cast<double>(below) / n_boot;
    prop = std::max(0.5 / n_boot, std::min(1.0 - 0.5 / n_boot, prop));
    double z0 = R::qnorm(prop, 0.0, 1.0, 1, 0);

    bool bca_ok = std::isfinite(z0);
    if (bca_ok) {
      // Step 2: jackknife acceleration
      std::unique_ptr<double[]> loo(new double[n - 1]);
      std::unique_ptr<double[]> jw1(new double[n]);
      std::unique_ptr<double[]> jw2(new double[n]);
      std::vector<double> jv(static_cast<size_t>(n));
      for (int i = 0; i < n; ++i) {
        int k = 0;
        for (int ii = 0; ii < n; ++ii)
          if (ii != i) loo[k++] = xp[ii];
        double res[N_ESTIMATORS];
        compute_all_estimators(loo.get(), n - 1, res, jw1.get(), jw2.get());
        jv[i] = res[kEstIdToAllIdx[estimator_id]];
      }
      double jack_mean = 0.0;
      for (int i = 0; i < n; ++i) jack_mean += jv[i];
      jack_mean /= n;
      double sum2 = 0.0, sum3 = 0.0;
      for (int i = 0; i < n; ++i) {
        double L  = jack_mean - jv[i];
        double L2 = L * L;
        sum2 += L2;
        sum3 += L2 * L;
      }
      double acc = (sum2 > 0.0) ? sum3 / (6.0 * std::pow(sum2, 1.5)) : 0.0;

      // Step 3: adjusted quantile indices
      double alpha  = 1.0 - level;
      double z_lo   = R::qnorm(alpha / 2.0,       0.0, 1.0, 1, 0);
      double z_hi   = R::qnorm(1.0 - alpha / 2.0, 0.0, 1.0, 1, 0);
      double num_lo = z0 + z_lo;
      double num_hi = z0 + z_hi;
      double a1 = R::pnorm(z0 + num_lo / (1.0 - acc * num_lo), 0.0, 1.0, 1, 0);
      double a2 = R::pnorm(z0 + num_hi / (1.0 - acc * num_hi), 0.0, 1.0, 1, 0);

      if (std::isfinite(a1) && std::isfinite(a2)) {
        a1 = std::max(0.5 / n_boot, std::min(1.0 - 0.5 / n_boot, a1));
        a2 = std::max(0.5 / n_boot, std::min(1.0 - 0.5 / n_boot, a2));
        std::sort(boot_vals.begin(), boot_vals.end());
        int idx1 = std::max(0, std::min(n_boot - 1,
                     static_cast<int>(std::floor(a1 * n_boot))));
        int idx2 = std::max(0, std::min(n_boot - 1,
                     static_cast<int>(std::floor(a2 * n_boot))));
        ci_lower = boot_vals[idx1];
        ci_upper = boot_vals[idx2];
      } else {
        bca_ok = false;
      }
    }
    if (!bca_ok) method_code = 1;  // fall back to percentile
  }

  if (method_code == 1) {
    // ====================== Percentile ======================
    std::sort(boot_vals.begin(), boot_vals.end());
    double alpha  = 1.0 - level;
    int lo_idx = std::max(0,
      static_cast<int>(std::floor(alpha / 2.0 * n_boot)));
    int hi_idx = std::min(n_boot - 1,
      static_cast<int>(std::floor((1.0 - alpha / 2.0) * n_boot)));
    ci_lower = boot_vals[lo_idx];
    ci_upper = boot_vals[hi_idx];
  } else if (method_code == 2) {
    // ====================== Parametric ======================
    double mean_b = 0.0;
    for (int r = 0; r < n_boot; ++r) mean_b += boot_vals[r];
    mean_b /= n_boot;
    double sq = 0.0;
    for (int r = 0; r < n_boot; ++r) {
      double d = boot_vals[r] - mean_b;
      sq += d * d;
    }
    double sd_b = (n_boot > 1) ? std::sqrt(sq / (n_boot - 1.0)) : 0.0;
    double z = R::qnorm(1.0 - (1.0 - level) / 2.0, 0.0, 1.0, 1, 0);
    ci_lower = est - z * sd_b;
    ci_upper = est + z * sd_b;
  }

  return Rcpp::List::create(
    Named("ci_lower") = ci_lower,
    Named("ci_upper") = ci_upper
  );
}

// [[Rcpp::export]]
double cpp_scale_ensemble(Rcpp::NumericVector x, int n_boot) {
  int n = x.size();
  if (n < 2) return NA_REAL;
  if (n_boot < 2) n_boot = 2;

  EnsembleCore core;
  core.run(x.begin(), n, n_boot);
  return core.ensemble_est;
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

  // Reduce n_boot for parametric tier
  int actual_nboot = (method_code == 2) ? std::min(n_boot, 50) : n_boot;
  if (actual_nboot < 2) actual_nboot = 2;

  EnsembleCore core;
  core.run(x.begin(), n, actual_nboot);

  const double* boot_results = core.boot_results;
  const double* weights      = core.weights;
  const double* estimates    = core.estimates;
  const double* vars         = core.vars;
  double        ensemble_est = core.ensemble_est;

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

  // Compute ensemble bootstrap values (needed by all tiers).
  // WU-LAYOUT-1: boot_results[j*actual_nboot+r] — stride-actual_nboot per j.
  // Inner loop has 7 iterations so the scattered reads stay in L1/L2 cache.
  std::vector<double> ensemble_boot(actual_nboot);
  for (int r = 0; r < actual_nboot; ++r) {
    double val = 0.0;
    for (int j = 0; j < N_ESTIMATORS; ++j) {
      double v = boot_results[static_cast<size_t>(j) * actual_nboot + r];
      if (std::isfinite(v) && v > 0.0) val += weights[j] * v;
    }
    ensemble_boot[r] = val;
  }

  if (method_code == 0) {
    // ========================= BCa =========================
    // Step 1: Bias correction z0
    // WU-LAYOUT-1: col = boot_results + j*actual_nboot — sequential reads.
    for (int j = 0; j < N_ESTIMATORS; ++j) {
      const double* col_j = boot_results + static_cast<size_t>(j) * actual_nboot;
      int below = 0;
      for (int r = 0; r < actual_nboot; ++r) {
        if (col_j[r] < estimates[j]) ++below;
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
          if (ii != i) loo[k++] = x[ii];
        compute_all_estimators(loo.get(), n - 1,
                               &jack_flat[static_cast<size_t>(i) * N_ESTIMATORS],
                               jw1.get(), jw2.get());
      }
    }

    // Per-estimator acceleration
    // kEstIdToAllIdx remaps API estimator_id order (0=gmd,1=sd) to
    // compute_all_estimators order (0=sd_c4,1=gmd) in the jack_flat array.
    for (int j = 0; j < N_ESTIMATORS; ++j) {
      const int jj = kEstIdToAllIdx[j];
      double jack_mean = 0.0;
      for (int i = 0; i < n; ++i)
        jack_mean += jack_flat[static_cast<size_t>(i) * N_ESTIMATORS + jj];
      jack_mean /= n;
      double sum2 = 0.0, sum3 = 0.0;
      for (int i = 0; i < n; ++i) {
        double L = jack_mean -
                   jack_flat[static_cast<size_t>(i) * N_ESTIMATORS + jj];
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

    // WU-LAYOUT-1: copy column j (sequential) then sort.
    std::vector<double> col(actual_nboot);
    for (int j = 0; j < N_ESTIMATORS; ++j) {
      const double* col_j = boot_results + static_cast<size_t>(j) * actual_nboot;
      std::copy(col_j, col_j + actual_nboot, col.begin());
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

    // WU-LAYOUT-1: copy sequential column j then sort.
    std::vector<double> col(actual_nboot);
    for (int j = 0; j < N_ESTIMATORS; ++j) {
      const double* col_j = boot_results + static_cast<size_t>(j) * actual_nboot;
      std::copy(col_j, col_j + actual_nboot, col.begin());
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

// Sorted-input variants — thin wrappers, always available for testing.
// [[Rcpp::export]]
double C_sn_sorted(Rcpp::NumericVector x) {
  return robscale::qnsn::C_sn_impl_sorted<double>(x.begin(), static_cast<size_t>(x.size()));
}

// [[Rcpp::export]]
double C_qn_sorted(Rcpp::NumericVector x) {
  return robscale::qnsn::C_qn_impl_sorted<double>(x.begin(), static_cast<size_t>(x.size()));
}
