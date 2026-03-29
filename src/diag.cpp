// diag.cpp — Internal diagnostic helper for M-scale iteration analysis.
// rob_scale_diag_impl mirrors rob_scale_core (Newton-Raphson) with
// instrumented iteration counts. Used by test-robScale-oscillation.R.

#include "robust_core.h"
#include "pdq_select.h"
#include "robscale_config.h"
#include <Rcpp.h>
#include <vector>
#include <cmath>

// ---------------------------------------------------------------------------
// M-scale iteration diagnostics
//
// Mirrors rob_scale_core/nr_scale_compute (Newton-Raphson) and counts:
//   - outer_iters: number of NR iteration steps
//   - aitken_fires: always 0 (NR has no Aitken; field retained for API compat)
//   - rho_evals: number of NR steps (= outer_iters)
//
// Called only via robscale:::rob_scale_diag_impl(x).
// ---------------------------------------------------------------------------

struct DiagStats {
  int outer_iters  = 0;
  int aitken_fires = 0;
  int rho_evals    = 0;
};

// Instrumented version of nr_scale_compute (Newton-Raphson).
// Mirrors the production scalar NR path exactly, with stat bookkeeping.
static double nr_scale_diag(const double* data, size_t n, double t,
                             double s, int maxit, double tol,
                             DiagStats& stats) {
  const double inv_n = 1.0 / static_cast<double>(n);
  const double hisc_pre = 0.5 * robscale::INV_RHO_SCALE_CONST;

  for (int iter = 0; iter < maxit; ++iter) {
    ++stats.outer_iters;
    ++stats.rho_evals;

    const double hisc = hisc_pre / s;
    double sum_tanh2 = 0.0;
    double sum_u_tanh_sech2 = 0.0;

    for (size_t i = 0; i < n; ++i) {
      double u = (data[i] - t) * hisc;
      double p = std::tanh(u);
      double p2 = p * p;
      sum_tanh2 += p2;
      sum_u_tanh_sech2 += u * p * (1.0 - p2);
    }

    // Match production: numer = sum_tanh2 * inv_n - 0.5
    const double numer = sum_tanh2 * inv_n - 0.5;
    const double denom = 2.0 * inv_n * sum_u_tanh_sech2;

    // GP guard: denominator degenerate → multiplicative fallback
    if (std::abs(denom) <= 1e-14 * s) {
      s *= std::sqrt(2.0 * sum_tanh2 * inv_n);
      continue;
    }

    const double delta_s = s * numer / denom;

    // Neg-s guard: proposed update non-positive → halve
    if (s + delta_s <= 0.0) { s /= 2.0; continue; }

    s += delta_s;
    if (std::abs(delta_s) / s <= tol) break;
  }
  return s;
}

// [[Rcpp::export(rng = false)]]
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

  // Median: n < 8 direct median_net (no SIMD kernel); n >= 8 through
  // adaptive dispatch (SIMD fires at 8/16/32, median_net for other small n).
  double t;
  if (n < 8) {
    t = robscale::median_net(w.data(), n);
  } else {
    t = robscale::adaptive_robscale_median_select(w.data(), n);
  }

  for (size_t i = 0; i < n; ++i) dev[i] = std::abs(w[i] - t);

  double s_init;
  if (n < 8) {
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

  // --- Instrumented NR iteration ---
  DiagStats stats;
  const double* data = x_r.begin();

  double s_final = nr_scale_diag(data, n, t, s_init, maxit, tol, stats);

  // Determine convergence: check final NR step size
  const double hisc_final = 0.5 * robscale::INV_RHO_SCALE_CONST / s_final;
  double sum_tanh2_f = 0.0, sum_ut_f = 0.0;
  for (size_t i = 0; i < n; ++i) {
    double u = (data[i] - t) * hisc_final;
    double p = std::tanh(u);
    double p2 = p * p;
    sum_tanh2_f += p2;
    sum_ut_f += u * p * (1.0 - p2);
  }
  double inv_n = 1.0 / static_cast<double>(n);
  double denom_f = 2.0 * inv_n * sum_ut_f;
  bool converged = (std::abs(denom_f) <= 1e-14 * s_final) ||
                   (std::abs(s_final * (1.0 - 2.0 * inv_n * sum_tanh2_f) / denom_f) / s_final <= tol);

  return Rcpp::List::create(
    Rcpp::Named("scale")        = s_final,
    Rcpp::Named("outer_iters")  = stats.outer_iters,
    Rcpp::Named("aitken_fires") = stats.aitken_fires,
    Rcpp::Named("rho_evals")    = stats.rho_evals,
    Rcpp::Named("converged")    = converged
  );
}
