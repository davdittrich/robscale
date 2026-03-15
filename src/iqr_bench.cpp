// iqr_bench.cpp — Benchmark-only IQR selection-algorithm variants
//
// Six implementations of Type 7 scaled IQR, each exposed via Rcpp::export.
// Only linked when the package is installed; not part of the public API.

#include "robscale_config.h"
#include "selection.h"             // robscale::floyd_rivest_select
#include "miniselect/pdqselect.h"  // miniselect::pdqselect{,_branchless}
#include <Rcpp.h>
#include <algorithm>               // std::nth_element, std::min_element
#include <cstring>
#include <memory>

// ---------- helpers ----------------------------------------------------------

// Type 7 quantile constant
static constexpr double IQR_K = 0.741301109252801;

// Arena allocation: micro / stack / heap tiering (matches iqr.cpp)
struct Arena {
  double  micro[ROBSCALE_MICRO_BUFFER_SIZE];
  double  stack[4096];              // generous stack tier
  std::unique_ptr<double[]> heap;
  double* buf;

  explicit Arena(int need) {
    if (need <= ROBSCALE_MICRO_BUFFER_SIZE)
      buf = micro;
    else if (need <= 4096)
      buf = stack;
    else {
      heap.reset(new double[need]);
      buf = heap.get();
    }
  }
};

// Type 7 quantile interpolation: given a partitioned buffer where buf[lo] is
// the lo-th order statistic, compute the interpolated quantile value.
static inline double interp_q7(double* buf, int n, int lo, double frac) {
  double q = buf[lo];
  if (frac > 0.0 && lo + 1 < n) {
    // next order statistic = min of buf[lo+1 .. n-1]
    double nv = buf[lo + 1];
    for (int i = lo + 2; i < n; ++i)
      if (buf[i] < nv) nv = buf[i];
    q += frac * (nv - q);
  }
  return q;
}

// Compute Q1/Q3 target indices for Type 7 quantile
struct Q13Idx {
  int lo1, lo3;
  double frac1, frac3;
};

static inline Q13Idx q13_indices(int n) {
  double h1 = (n - 1.0) * 0.25;
  double h3 = (n - 1.0) * 0.75;
  return {
    static_cast<int>(h1),
    static_cast<int>(h3),
    h1 - static_cast<int>(h1),
    h3 - static_cast<int>(h3)
  };
}

// ---------- Variant 1: current (dual Floyd-Rivest, 2 memcpys) ---------------

// [[Rcpp::export]]
double iqr_current(Rcpp::NumericVector x) {
  int n = x.size();
  if (n < 2) return 0.0;

  Arena a(n * 2);
  double* buf1 = a.buf;
  double* buf2 = a.buf + n;
  auto idx = q13_indices(n);

  std::memcpy(buf1, x.begin(), n * sizeof(double));
  robscale::floyd_rivest_select(buf1, buf1 + idx.lo1, buf1 + n);
  double q1 = interp_q7(buf1, n, idx.lo1, idx.frac1);

  std::memcpy(buf2, x.begin(), n * sizeof(double));
  robscale::floyd_rivest_select(buf2, buf2 + idx.lo3, buf2 + n);
  double q3 = interp_q7(buf2, n, idx.lo3, idx.frac3);

  return (q3 - q1) * IQR_K;
}

// ---------- Variant 2: FR incremental (single copy, Q3 on partition) --------

// [[Rcpp::export]]
double iqr_fr_incremental(Rcpp::NumericVector x) {
  int n = x.size();
  if (n < 2) return 0.0;

  Arena a(n);
  double* buf = a.buf;
  auto idx = q13_indices(n);

  std::memcpy(buf, x.begin(), n * sizeof(double));

  // Q1: full array
  robscale::floyd_rivest_select(buf, buf + idx.lo1, buf + n);
  double q1 = interp_q7(buf, n, idx.lo1, idx.frac1);

  // Q3: only on buf[lo1+1 .. n-1] — everything <= Q1 is irrelevant
  int start = idx.lo1 + 1;
  robscale::floyd_rivest_select(buf + start, buf + idx.lo3, buf + n);
  double q3 = interp_q7(buf, n, idx.lo3, idx.frac3);

  return (q3 - q1) * IQR_K;
}

// ---------- Variant 3: pdqselect dual (2 memcpys) ---------------------------

// [[Rcpp::export]]
double iqr_pdq_dual(Rcpp::NumericVector x) {
  int n = x.size();
  if (n < 2) return 0.0;

  Arena a(n * 2);
  double* buf1 = a.buf;
  double* buf2 = a.buf + n;
  auto idx = q13_indices(n);

  std::memcpy(buf1, x.begin(), n * sizeof(double));
  miniselect::pdqselect(buf1, buf1 + idx.lo1, buf1 + n);
  double q1 = interp_q7(buf1, n, idx.lo1, idx.frac1);

  std::memcpy(buf2, x.begin(), n * sizeof(double));
  miniselect::pdqselect(buf2, buf2 + idx.lo3, buf2 + n);
  double q3 = interp_q7(buf2, n, idx.lo3, idx.frac3);

  return (q3 - q1) * IQR_K;
}

// ---------- Variant 4: pdqselect incremental --------------------------------

// [[Rcpp::export]]
double iqr_pdq_incremental(Rcpp::NumericVector x) {
  int n = x.size();
  if (n < 2) return 0.0;

  Arena a(n);
  double* buf = a.buf;
  auto idx = q13_indices(n);

  std::memcpy(buf, x.begin(), n * sizeof(double));

  miniselect::pdqselect(buf, buf + idx.lo1, buf + n);
  double q1 = interp_q7(buf, n, idx.lo1, idx.frac1);

  int start = idx.lo1 + 1;
  miniselect::pdqselect(buf + start, buf + idx.lo3, buf + n);
  double q3 = interp_q7(buf, n, idx.lo3, idx.frac3);

  return (q3 - q1) * IQR_K;
}

// ---------- Variant 5: pdqselect_branchless incremental ---------------------

// [[Rcpp::export]]
double iqr_pdq_branchless_inc(Rcpp::NumericVector x) {
  int n = x.size();
  if (n < 2) return 0.0;

  Arena a(n);
  double* buf = a.buf;
  auto idx = q13_indices(n);

  std::memcpy(buf, x.begin(), n * sizeof(double));

  miniselect::pdqselect_branchless(buf, buf + idx.lo1, buf + n);
  double q1 = interp_q7(buf, n, idx.lo1, idx.frac1);

  int start = idx.lo1 + 1;
  miniselect::pdqselect_branchless(buf + start, buf + idx.lo3, buf + n);
  double q3 = interp_q7(buf, n, idx.lo3, idx.frac3);

  return (q3 - q1) * IQR_K;
}

// ---------- Variant 6: std::nth_element incremental -------------------------

// [[Rcpp::export]]
double iqr_nth_incremental(Rcpp::NumericVector x) {
  int n = x.size();
  if (n < 2) return 0.0;

  Arena a(n);
  double* buf = a.buf;
  auto idx = q13_indices(n);

  std::memcpy(buf, x.begin(), n * sizeof(double));

  std::nth_element(buf, buf + idx.lo1, buf + n);
  double q1 = interp_q7(buf, n, idx.lo1, idx.frac1);

  int start = idx.lo1 + 1;
  std::nth_element(buf + start, buf + idx.lo3, buf + n);
  double q3 = interp_q7(buf, n, idx.lo3, idx.frac3);

  return (q3 - q1) * IQR_K;
}
