// mad_bench.cpp — Benchmark-only MAD selection-algorithm variants
//
// Five implementations of constant * median(|x - median(x)|), each exposed
// via Rcpp::export.  Only linked when the package is installed; not part of
// the public API.

#include "robscale_config.h"
#include "robust_core.h"                    // robscale::median_select
#include "selection.h"                      // robscale::floyd_rivest_select
#include "miniselect/pdqselect.h"           // miniselect::pdqselect{,_branchless}
#include "miniselect/floyd_rivest_select.h" // miniselect::floyd_rivest_select
#include <Rcpp.h>
#include <algorithm>                        // std::nth_element
#include <cstring>
#include <memory>

// ---------- arena allocation (matches mad.cpp tiering) ----------------------

struct MadArena {
  double  micro[ROBSCALE_MICRO_BUFFER_SIZE];
  static constexpr int STACK_SIZE = 2048;
  double  stack[STACK_SIZE * 2];
  std::unique_ptr<double[]> heap;
  double* w;    // copy of input  (n doubles)
  double* dev;  // deviation buf  (n doubles)

  explicit MadArena(int n) {
    double* arena;
    if (n <= 64) {
      arena = micro;
    } else if (n <= STACK_SIZE) {
      arena = stack;
    } else {
      heap.reset(new double[n * 2]);
      arena = heap.get();
    }
    w   = arena;
    dev = arena + n;
  }
};

// ---------- median helpers per algorithm ------------------------------------

// Even-n median: after selecting the (n-1)/2-th element, scan for min of the
// upper partition to get the (n/2)-th element, then average.  Odd-n: just the
// middle element.

// Helper: scan for min in buf[start..n-1]
static inline double scan_min(const double* buf, int start, int n) {
  double v = buf[start];
  for (int i = start + 1; i < n; ++i)
    if (buf[i] < v) v = buf[i];
  return v;
}

// --- current production path (sorting nets + our FR) ---
static inline double median_current(double* buf, int n) {
  return robscale::median_select(buf, static_cast<size_t>(n));
}

// --- pdqselect ---
static inline double median_pdq(double* buf, int n) {
  int h = (n - 1) / 2;
  miniselect::pdqselect(buf, buf + h, buf + n);
  if (n & 1) return buf[h];
  return (buf[h] + scan_min(buf, h + 1, n)) * 0.5;
}

// --- pdqselect_branchless ---
static inline double median_pdq_branchless(double* buf, int n) {
  int h = (n - 1) / 2;
  miniselect::pdqselect_branchless(buf, buf + h, buf + n);
  if (n & 1) return buf[h];
  return (buf[h] + scan_min(buf, h + 1, n)) * 0.5;
}

// --- std::nth_element ---
static inline double median_nth(double* buf, int n) {
  int h = (n - 1) / 2;
  std::nth_element(buf, buf + h, buf + n);
  if (n & 1) return buf[h];
  return (buf[h] + scan_min(buf, h + 1, n)) * 0.5;
}

// --- miniselect FR ---
static inline double median_miniselect_fr(double* buf, int n) {
  int h = (n - 1) / 2;
  miniselect::floyd_rivest_select(buf, buf + h, buf + n);
  if (n & 1) return buf[h];
  return (buf[h] + scan_min(buf, h + 1, n)) * 0.5;
}

// ---------- MAD template (DRY) ----------------------------------------------

template <double (*MedianFn)(double*, int)>
static double mad_variant(const double* xp, int n) {
  if (n < 1) return NA_REAL;
  if (n == 1) return 0.0;

  MadArena a(n);
  std::memcpy(a.w, xp, n * sizeof(double));
  double med = MedianFn(a.w, n);

  for (int i = 0; i < n; ++i) a.dev[i] = std::abs(xp[i] - med);
  double mad_raw = MedianFn(a.dev, n);
  return robscale::MAD_CONSISTENCY * mad_raw;
}

// ---------- Rcpp-exported variants ------------------------------------------

// [[Rcpp::export]]
double mad_bench_current(Rcpp::NumericVector x) {
  return mad_variant<median_current>(x.begin(), x.size());
}

// [[Rcpp::export]]
double mad_bench_pdq(Rcpp::NumericVector x) {
  return mad_variant<median_pdq>(x.begin(), x.size());
}

// [[Rcpp::export]]
double mad_bench_pdq_branchless(Rcpp::NumericVector x) {
  return mad_variant<median_pdq_branchless>(x.begin(), x.size());
}

// [[Rcpp::export]]
double mad_bench_nth(Rcpp::NumericVector x) {
  return mad_variant<median_nth>(x.begin(), x.size());
}

// [[Rcpp::export]]
double mad_bench_miniselect_fr(Rcpp::NumericVector x) {
  return mad_variant<median_miniselect_fr>(x.begin(), x.size());
}
