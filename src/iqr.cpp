#include "robscale_config.h"
#include "robust_core.h"
#include "selection.h"
#include <Rcpp.h>
#include <cstring>
#include <memory>

// Type 7 quantile via selection (R default)
static ROBSCALE_INLINE double quantile7_select(double* buf, int n, double p) {
  double h = (n - 1.0) * p;
  int lo = static_cast<int>(h);
  double frac = h - lo;

  robscale::floyd_rivest_select(buf, buf + lo, buf + n);
  double q = buf[lo];

  if (frac > 0.0 && lo + 1 < n) {
    // Find minimum of elements after lo (the next order statistic)
    double next_val = buf[lo + 1];
    for (int i = lo + 2; i < n; ++i) {
      if (buf[i] < next_val) next_val = buf[i];
    }
    q += frac * (next_val - q);
  }
  return q;
}

// [[Rcpp::export]]
double iqr_impl(Rcpp::NumericVector x, double constant) {
  int n = x.size();
  if (n < 2) return 0.0;

  // Arena allocation: need two copies since selection is destructive
  double buf_micro[ROBSCALE_MICRO_BUFFER_SIZE];
  constexpr int STACK_SIZE = 2048;
  double buf_stack[STACK_SIZE];
  std::unique_ptr<double[]> heap;
  double* buf;

  if (n <= ROBSCALE_MICRO_BUFFER_SIZE / 2) {
    buf = buf_micro;
  } else if (n <= STACK_SIZE / 2) {
    buf = buf_stack;
  } else {
    heap.reset(new double[n * 2]);
    buf = heap.get();
  }

  double* buf2 = buf + n;

  // Q1 on first copy
  std::memcpy(buf, x.begin(), n * sizeof(double));
  double q1 = quantile7_select(buf, n, 0.25);

  // Q3 on second copy (selection is destructive)
  std::memcpy(buf2, x.begin(), n * sizeof(double));
  double q3 = quantile7_select(buf2, n, 0.75);

  return (q3 - q1) * constant;
}
