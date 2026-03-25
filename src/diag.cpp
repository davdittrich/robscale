// diag.cpp — Internal diagnostic and benchmark helpers.
// None of these functions are part of the user-facing API.
// They are compiled into the package to support:
//   M-scale iteration diagnostics: convergence and rho_eval counting for rob_scale
//   Median crossover benchmarks: median_net vs FR-select threshold calibration
//
// All exports are intentionally ugly names (rob_scale_diag_impl,
// bench_median_net_impl, bench_fr_select_impl) to signal internal status.

#include "robust_core.h"
#include "sort_net.h"
#include "pdq_select.h"
#include "robscale_config.h"
#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <algorithm>

// ---------------------------------------------------------------------------
// M-scale iteration diagnostics
//
// Mirrors rob_scale_core/rob_scale_compute exactly but counts:
//   - outer_iters: number of while-loop passes in aitken_iterate
//   - aitken_fires: number of times Aitken extrapolation was accepted
//   - rho_evals: total rho_sum calls (= 2 * outer_iters, always)
//
// Called only via robscale:::rob_scale_diag_impl(x).
// ---------------------------------------------------------------------------

struct DiagStats {
  int outer_iters  = 0;
  int aitken_fires = 0;
  int rho_evals    = 0;
};

// Instrumented version of aitken_iterate.  Identical logic; adds stat
// bookkeeping.  Kept here (not in rob_scale.cpp) to avoid touching
// production code before the RED tests are written.
template <typename RhoSum>
static double aitken_iterate_diag(RhoSum&& rho_sum, double s,
                                   int maxit, double tol, double inv_n,
                                   DiagStats& stats) {
  int k = 0;
  while (k < maxit) {
    ++stats.outer_iters;

    const double s0 = s;
    const double v0 = std::sqrt(2.0 * rho_sum(s0) * inv_n);
    ++stats.rho_evals;
    ++k;
    const double s1 = s0 * v0;
    if (std::abs(v0 - 1.0) <= tol) { s = s1; break; }
    if (k >= maxit)                 { s = s1; break; }

    const double v1 = std::sqrt(2.0 * rho_sum(s1) * inv_n);
    ++stats.rho_evals;
    ++k;
    const double s2 = s1 * v1;
    if (std::abs(v1 - 1.0) <= tol) { s = s2; break; }

    const double d1    = s1 - s0;
    const double d2    = s2 - s1;
    const double denom = d2 - d1;
    if (d1 * d2 > 0.0 && std::abs(d2) < std::abs(d1) &&
        std::abs(denom) > 1e-30 * s0 && k < maxit) {
      const double candidate = s2 - d2 * d2 / denom;
      if (candidate > 0.0) {
        ++stats.aitken_fires;
        s = candidate;
        continue;
      }
    }
    s = (d1 * d2 < 0.0) ? std::sqrt(s1 * s2) : s2;
  }
  return s;
}

// [[Rcpp::export]]
Rcpp::List rob_scale_diag_impl(Rcpp::NumericVector x_r,
                                int maxit = 80,
                                double tol = 1.4901161193847656e-8) {
  size_t n = static_cast<size_t>(x_r.size());
  if (n < 2) {
    return Rcpp::List::create(
      Rcpp::Named("scale")        = 0.0,
      Rcpp::Named("outer_iters")  = 0,
      Rcpp::Named("aitken_fires") = 0,
      Rcpp::Named("rho_evals")    = 0,
      Rcpp::Named("converged")    = true
    );
  }

  // --- Initialisation: mirror rob_scale_core exactly ---
  std::vector<double> w(x_r.begin(), x_r.end());
  std::vector<double> dev(n);

  // Median: for n > ROBSCALE_SORT_NETWORK_THRESHOLD use FR-based select,
  // else use median_net.  Must match the current threshold.
  const bool is_small = (n <= ROBSCALE_SORT_NETWORK_THRESHOLD);
  double t;
  if (is_small) {
    t = robscale::median_net(w.data(), n);
  } else {
    t = robscale::adaptive_robscale_median_select(w.data(), n);
  }

  for (size_t i = 0; i < n; ++i) dev[i] = std::abs(w[i] - t);

  double s_init;
  if (is_small) {
    s_init = robscale::MAD_CONSISTENCY * robscale::median_net(dev.data(), n);
  } else {
    s_init = robscale::MAD_CONSISTENCY *
             robscale::adaptive_robscale_median_select(dev.data(), n);
  }

  if (s_init <= 0.0) {
    return Rcpp::List::create(
      Rcpp::Named("scale")        = 0.0,
      Rcpp::Named("outer_iters")  = 0,
      Rcpp::Named("aitken_fires") = 0,
      Rcpp::Named("rho_evals")    = 0,
      Rcpp::Named("converged")    = true
    );
  }

  // --- Instrumented iteration ---
  const double inv_n = 1.0 / static_cast<double>(n);
  DiagStats stats;

  // Use scalar rho_sum (matches rob_scale_compute scalar path for portability)
  const double* data = x_r.begin();
  auto rho_sum = [&](double sc) -> double {
    const double hisc = 0.5 * robscale::INV_RHO_SCALE_CONST / sc;
    std::vector<double> tmp(n);
    for (size_t i = 0; i < n; ++i) tmp[i] = (data[i] - t) * hisc;
    robscale::bulk_tanh(tmp.data(), static_cast<int>(n));
    double sr = 0.0;
    for (size_t i = 0; i < n; ++i) sr += tmp[i] * tmp[i];
    return sr;
  };

  double s_final = aitken_iterate_diag(rho_sum, s_init, maxit, tol, inv_n, stats);

  // Determine convergence: re-evaluate v at s_final
  const double hisc_final = 0.5 * robscale::INV_RHO_SCALE_CONST / s_final;
  std::vector<double> tmp(n);
  for (size_t i = 0; i < n; ++i) tmp[i] = (data[i] - t) * hisc_final;
  robscale::bulk_tanh(tmp.data(), static_cast<int>(n));
  double sr_final = 0.0;
  for (size_t i = 0; i < n; ++i) sr_final += tmp[i] * tmp[i];
  const double v_final = std::sqrt(2.0 * sr_final * inv_n);
  bool converged = std::abs(v_final - 1.0) <= tol;

  return Rcpp::List::create(
    Rcpp::Named("scale")        = s_final,
    Rcpp::Named("outer_iters")  = stats.outer_iters,
    Rcpp::Named("aitken_fires") = stats.aitken_fires,
    Rcpp::Named("rho_evals")    = stats.rho_evals,
    Rcpp::Named("converged")    = converged
  );
}

// ---------------------------------------------------------------------------
// Median crossover benchmark helpers
//
// Each function copies x into a local buffer, runs the target selection
// algorithm, and returns the median.  The copy ensures:
//   (a) the input vector is not permuted (selection algorithms are in-place)
//   (b) cache state is comparable across repeated calls
//
// Called via:
//   robscale:::bench_median_net_impl(x)
//   robscale:::bench_fr_select_impl(x)
// ---------------------------------------------------------------------------

// [[Rcpp::export]]
double bench_median_net_impl(Rcpp::NumericVector x) {
  std::vector<double> buf(x.begin(), x.end());
  return robscale::median_net(buf.data(), buf.size());
}

// [[Rcpp::export]]
double bench_fr_select_impl(Rcpp::NumericVector x) {
  // Uses median_select from robust_core.h which:
  //   n <= ROBSCALE_SORT_NETWORK_THRESHOLD: calls median_net (same as above)
  //   n >  ROBSCALE_SORT_NETWORK_THRESHOLD: calls floyd_rivest_select
  //     -> which for n < 600 uses std::nth_element
  // To force the FR/nth_element path for n <= threshold we call
  // floyd_rivest_select directly rather than going through median_select.
  std::vector<double> buf(x.begin(), x.end());
  size_t n = buf.size();
  if (n == 0) return 0.0;
  size_t h = (n - 1) / 2;
  robscale::floyd_rivest_select(buf.data(), buf.data() + h, buf.data() + n);
  if (n & 1) return buf[h];
  double v1 = buf[h];
  double v2 = *std::min_element(buf.data() + h + 1, buf.data() + n);
  return (v1 + v2) * 0.5;
}
