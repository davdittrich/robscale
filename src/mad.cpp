#include "robscale_config.h"
#include "pdq_select.h"
#include <Rcpp.h>
#include <cstring>
#include <memory>

// MAD with auto-median: fused single-buffer approach.
// After median selection, w is a permutation of x. Compute deviations
// in-place: {|w[i] - med|} == {|x[i] - med|} as multisets, and MAD
// is order-invariant. This halves memory from 2n to n doubles.
// [[Rcpp::export]]
double mad_impl_auto(Rcpp::NumericVector x, double constant) {
  int n = x.size();
  if (n < 1) return NA_REAL;
  if (n == 1) return 0.0;

  const double* xp = x.begin();

  // Single buffer of size n (was 2n)
  double buf_micro[64];
  constexpr int STACK_SIZE = 2048;
  double buf_stack[STACK_SIZE];
  std::unique_ptr<double[]> heap;
  double* w;

  if (n <= 64) {
    w = buf_micro;
  } else if (n <= STACK_SIZE) {
    w = buf_stack;
  } else {
    heap.reset(new double[n]);
    w = heap.get();
  }

  // Step 1: copy and select median (destroys ordering of w)
  std::memcpy(w, xp, n * sizeof(double));
  double med = robscale::adaptive_median_select(w, n);

  // Step 2: compute absolute deviations in-place
  for (int i = 0; i < n; ++i) w[i] = std::abs(w[i] - med);

  // Step 3: select median of deviations
  double mad_raw = robscale::adaptive_median_select(w, static_cast<size_t>(n));
  return constant * mad_raw;
}

// MAD with user-supplied center
// [[Rcpp::export]]
double mad_impl_center(Rcpp::NumericVector x, double center, double constant) {
  int n = x.size();
  if (n < 1) return NA_REAL;
  if (n == 1) return 0.0;

  const double* xp = x.begin();

  double buf_micro[64];
  constexpr int STACK_SIZE = 2048;
  double buf_stack[STACK_SIZE];
  std::unique_ptr<double[]> heap;
  double* dev;

  if (n <= 64) {
    dev = buf_micro;
  } else if (n <= STACK_SIZE) {
    dev = buf_stack;
  } else {
    heap.reset(new double[n]);
    dev = heap.get();
  }

  for (int i = 0; i < n; ++i) dev[i] = std::abs(xp[i] - center);
  double mad_raw = robscale::adaptive_median_select(dev, static_cast<size_t>(n));
  return constant * mad_raw;
}
