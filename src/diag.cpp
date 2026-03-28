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

// [[Rcpp::export(rng = false)]]
double bench_median_net_impl(Rcpp::NumericVector x) {
  std::vector<double> buf(x.begin(), x.end());
  return robscale::median_net(buf.data(), buf.size());
}

// [[Rcpp::export(rng = false)]]
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
