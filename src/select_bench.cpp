// select_bench.cpp — Direct FR vs pdqselect timing for threshold calibration.
//
// These functions isolate the selection step from the surrounding estimator
// logic so callers can find the true FR/pdqselect crossover n for each
// call site (median, low-median, arbitrary-k) without interference from
// tanh iteration, inner_medians computation, or refinement loops.
//
// Each function copies the input to a local heap buffer so that bench::mark()
// can call it repeatedly with a clean, unsorted array every time.

#include "robscale_config.h"
#include "selection.h"            // robscale::floyd_rivest_select
#include "miniselect/pdqselect.h"
#include <Rcpp.h>
#include <memory>
#include <cstring>
#include <cmath>

// ---------------------------------------------------------------------------
// Median (n/2 order statistic, low-median convention) — robScale / IQR / MAD
// ---------------------------------------------------------------------------

// [[Rcpp::export]]
double sel_fr_median(Rcpp::NumericVector x) {
  size_t n = static_cast<size_t>(x.size());
  if (n == 0) return 0.0;
  std::unique_ptr<double[]> buf(new double[n]);
  std::memcpy(buf.get(), x.begin(), n * sizeof(double));
  size_t h = (n - 1) / 2;
  robscale::floyd_rivest_select(buf.get(), buf.get() + h, buf.get() + n);
  return buf[h];
}

// [[Rcpp::export]]
double sel_pdq_median(Rcpp::NumericVector x) {
  size_t n = static_cast<size_t>(x.size());
  if (n == 0) return 0.0;
  std::unique_ptr<double[]> buf(new double[n]);
  std::memcpy(buf.get(), x.begin(), n * sizeof(double));
  size_t h = (n - 1) / 2;
  miniselect::pdqselect(buf.get(), buf.get() + h, buf.get() + n);
  return buf[h];
}

// ---------------------------------------------------------------------------
// Low-median (no even-n averaging) — Sn inner_medians
// ---------------------------------------------------------------------------

// [[Rcpp::export]]
double sel_fr_lowmedian(Rcpp::NumericVector x) {
  size_t n = static_cast<size_t>(x.size());
  if (n == 0) return 0.0;
  std::unique_ptr<double[]> buf(new double[n]);
  std::memcpy(buf.get(), x.begin(), n * sizeof(double));
  size_t h = (n - 1) / 2;
  robscale::floyd_rivest_select(buf.get(), buf.get() + h, buf.get() + n);
  return buf[h];
}

// [[Rcpp::export]]
double sel_pdq_lowmedian(Rcpp::NumericVector x) {
  size_t n = static_cast<size_t>(x.size());
  if (n == 0) return 0.0;
  std::unique_ptr<double[]> buf(new double[n]);
  std::memcpy(buf.get(), x.begin(), n * sizeof(double));
  size_t h = (n - 1) / 2;
  miniselect::pdqselect(buf.get(), buf.get() + h, buf.get() + n);
  return buf[h];
}

// ---------------------------------------------------------------------------
// Arbitrary-k selection — Qn final diff-window
// k is 0-based (k=0 selects the minimum)
// ---------------------------------------------------------------------------

// [[Rcpp::export]]
double sel_fr_kth(Rcpp::NumericVector x, int k) {
  size_t n = static_cast<size_t>(x.size());
  if (n == 0 || k < 0 || static_cast<size_t>(k) >= n) return R_NaReal;
  std::unique_ptr<double[]> buf(new double[n]);
  std::memcpy(buf.get(), x.begin(), n * sizeof(double));
  size_t kk = static_cast<size_t>(k);
  robscale::floyd_rivest_select(buf.get(), buf.get() + kk, buf.get() + n);
  return buf[kk];
}

// [[Rcpp::export]]
double sel_pdq_kth(Rcpp::NumericVector x, int k) {
  size_t n = static_cast<size_t>(x.size());
  if (n == 0 || k < 0 || static_cast<size_t>(k) >= n) return R_NaReal;
  std::unique_ptr<double[]> buf(new double[n]);
  std::memcpy(buf.get(), x.begin(), n * sizeof(double));
  size_t kk = static_cast<size_t>(k);
  miniselect::pdqselect(buf.get(), buf.get() + kk, buf.get() + n);
  return buf[kk];
}
