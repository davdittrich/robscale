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
#include "simd_median.h"
#include "pdq_select.h"
#include "robscale_config.h"
#include "qnsn_runtime_config.h"
#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <algorithm>
#include <chrono>
#include <cstring>

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

// ---------------------------------------------------------------------------
// SIMD bitonic median benchmark helpers
// ---------------------------------------------------------------------------

// Forces AVX2 SIMD path (falls back to scalar median_net on non-x86).
// [[Rcpp::export(rng = false)]]
double bench_simd_median_avx2_impl(Rcpp::NumericVector x) {
  std::vector<double> buf(x.begin(), x.end());
#ifdef ROBSCALE_HAS_AVX2_DISPATCH
  return robscale::simd::simd_median_dispatch_avx2(buf.data(), buf.size());
#else
  return robscale::median_net(buf.data(), buf.size());
#endif
}

// SIMD selection network (hybrid: SIMD early stages + scalar tail).
// [[Rcpp::export(rng = false)]]
double bench_simd_median_sel_impl(Rcpp::NumericVector x) {
  std::vector<double> buf(x.begin(), x.end());
#ifdef ROBSCALE_HAS_AVX2_DISPATCH
  return robscale::simd::simd_median_sel_dispatch_avx2(buf.data(), buf.size());
#else
  return robscale::median_net(buf.data(), buf.size());
#endif
}

// Auto-dispatch: uses best available SIMD tier, falls back to scalar.
// [[Rcpp::export(rng = false)]]
double bench_simd_median_impl(Rcpp::NumericVector x) {
  std::vector<double> buf(x.begin(), x.end());
  size_t n = buf.size();
#ifdef ROBSCALE_HAS_AVX2_DISPATCH
  const auto& cfg = robscale::qnsn::RuntimeConfig::get();
  if (cfg.hw.simd_level >= robscale::qnsn::SIMDLevel::AVX2)
    return robscale::simd::simd_median_dispatch_avx2(buf.data(), n);
#endif
  return robscale::median_net(buf.data(), n);
}

// ---------------------------------------------------------------------------
// Full-sort benchmark helpers: SIMD bitonic vs scalar sort_net vs std::sort
// ---------------------------------------------------------------------------

// SIMD bitonic full sort (AVX2). Pads to next power-of-2 with +inf.
// Returns buf[0] (minimum) as a correctness check value.
// [[Rcpp::export(rng = false)]]
double bench_simd_sort_impl(Rcpp::NumericVector x) {
  std::vector<double> buf(x.begin(), x.end());
  size_t n = buf.size();
#ifdef ROBSCALE_HAS_AVX2_DISPATCH
  robscale::simd::simd_sort_dispatch_avx2(buf.data(), n);
#else
  robscale::small_sort(buf.data(), n);
#endif
  return buf[0];
}

// SIMD hybrid sorting network (SIMD early stages + scalar tail).
// [[Rcpp::export(rng = false)]]
double bench_simd_sort_net_impl(Rcpp::NumericVector x) {
  std::vector<double> buf(x.begin(), x.end());
#ifdef ROBSCALE_HAS_AVX2_DISPATCH
  robscale::simd::simd_sort_net_dispatch_avx2(buf.data(), buf.size());
#else
  robscale::small_sort(buf.data(), buf.size());
#endif
  return buf[0];
}

// Scalar sorting network (sort_net_N via small_sort).
// [[Rcpp::export(rng = false)]]
double bench_sort_net_impl(Rcpp::NumericVector x) {
  std::vector<double> buf(x.begin(), x.end());
  robscale::small_sort(buf.data(), buf.size());
  return buf[0];
}

// std::sort baseline.
// [[Rcpp::export(rng = false)]]
double bench_std_sort_impl(Rcpp::NumericVector x) {
  std::vector<double> buf(x.begin(), x.end());
  std::sort(buf.data(), buf.data() + buf.size());
  return buf[0];
}

// ---------------------------------------------------------------------------
// Tight-loop sort benchmark: times only the sort, zero R overhead.
//
// Pre-fills `reps` copies of the input array, then times the sort loop
// using std::chrono.  Returns per-sort nanoseconds.
//
// method: 0 = scalar sort_net (small_sort)
//         1 = std::sort
//         2 = SIMD hybrid sort_net
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// Tight-loop median benchmark: times only the median, zero R overhead.
//
// method: 0 = scalar median_net
//         1 = Floyd-Rivest / nth_element
//         2 = SIMD hybrid selection network
// ---------------------------------------------------------------------------

// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector bench_median_timed_impl(Rcpp::NumericVector x,
                                            int reps, int method) {
  const size_t n = static_cast<size_t>(x.size());
  const size_t stride = n;

  std::vector<double> pool(stride * reps);
  for (int r = 0; r < reps; ++r)
    std::memcpy(pool.data() + r * stride, x.begin(), n * sizeof(double));

  volatile double sink = 0.0;

  auto t0 = std::chrono::high_resolution_clock::now();
  switch (method) {
    case 0:  // scalar median_net
      for (int r = 0; r < reps; ++r) {
        double* buf = pool.data() + r * stride;
        sink += robscale::median_net(buf, n);
      }
      break;
    case 1:  // Floyd-Rivest / nth_element
      for (int r = 0; r < reps; ++r) {
        double* buf = pool.data() + r * stride;
        if (n == 0) { sink += 0.0; break; }
        size_t h = (n - 1) / 2;
        robscale::floyd_rivest_select(buf, buf + h, buf + n);
        if (n & 1) { sink += buf[h]; }
        else {
          double v1 = buf[h];
          double v2 = *std::min_element(buf + h + 1, buf + n);
          sink += (v1 + v2) * 0.5;
        }
      }
      break;
    case 2:  // SIMD hybrid selection network
#ifdef ROBSCALE_HAS_AVX2_DISPATCH
      for (int r = 0; r < reps; ++r) {
        double* buf = pool.data() + r * stride;
        sink += robscale::simd::simd_median_sel_dispatch_avx2(buf, n);
      }
#else
      for (int r = 0; r < reps; ++r) {
        double* buf = pool.data() + r * stride;
        sink += robscale::median_net(buf, n);
      }
#endif
      break;
  }
  auto t1 = std::chrono::high_resolution_clock::now();

  double total_ns = std::chrono::duration<double, std::nano>(t1 - t0).count();
  (void)sink;
  return Rcpp::NumericVector::create(total_ns / reps);
}

// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector bench_sort_timed_impl(Rcpp::NumericVector x,
                                          int reps, int method) {
  const size_t n = static_cast<size_t>(x.size());
  const size_t stride = n;

  // Pre-fill reps copies (each sort is destructive, needs fresh data)
  std::vector<double> pool(stride * reps);
  for (int r = 0; r < reps; ++r)
    std::memcpy(pool.data() + r * stride, x.begin(), n * sizeof(double));

  volatile double sink = 0.0;  // prevent dead-code elimination

  auto t0 = std::chrono::high_resolution_clock::now();
  switch (method) {
    case 0:
      for (int r = 0; r < reps; ++r) {
        double* buf = pool.data() + r * stride;
        robscale::small_sort(buf, n);
        sink += buf[0];
      }
      break;
    case 1:
      for (int r = 0; r < reps; ++r) {
        double* buf = pool.data() + r * stride;
        std::sort(buf, buf + n);
        sink += buf[0];
      }
      break;
    case 2:
#ifdef ROBSCALE_HAS_AVX2_DISPATCH
      for (int r = 0; r < reps; ++r) {
        double* buf = pool.data() + r * stride;
        robscale::simd::simd_sort_net_dispatch_avx2(buf, n);
        sink += buf[0];
      }
#else
      for (int r = 0; r < reps; ++r) {
        double* buf = pool.data() + r * stride;
        robscale::small_sort(buf, n);
        sink += buf[0];
      }
#endif
      break;
  }
  auto t1 = std::chrono::high_resolution_clock::now();

  double total_ns = std::chrono::duration<double, std::nano>(t1 - t0).count();
  double per_sort_ns = total_ns / reps;
  (void)sink;
  return Rcpp::NumericVector::create(per_sort_ns);
}
