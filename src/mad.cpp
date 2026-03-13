#include "robscale_config.h"
#include "robust_core.h"
#include <Rcpp.h>
#include <cstring>
#include <memory>

// MAD with auto-median: uses existing mad_select() infrastructure
// [[Rcpp::export]]
double mad_impl_auto(Rcpp::NumericVector x, double constant) {
  int n = x.size();
  if (n < 1) return NA_REAL;
  if (n == 1) return 0.0;

  const double* xp = x.begin();

  // Arena: need space for median copy + deviation buffer
  double buf_micro[ROBSCALE_MICRO_BUFFER_SIZE];
  constexpr int STACK_SIZE = 2048;
  double buf_stack[STACK_SIZE * 2];
  std::unique_ptr<double[]> heap;
  double* arena;

  if (n <= 64) {
    arena = buf_micro;
  } else if (n <= STACK_SIZE) {
    arena = buf_stack;
  } else {
    heap.reset(new double[n * 2]);
    arena = heap.get();
  }

  double* w = arena;
  double* dev = arena + n;

  // Compute median via selection
  std::memcpy(w, xp, n * sizeof(double));
  double med = robscale::median_select(w, n);

  // Compute MAD using existing kernel (uses constant from parameter)
  for (int i = 0; i < n; ++i) dev[i] = std::abs(xp[i] - med);
  double mad_raw = robscale::median_select(dev, static_cast<size_t>(n));
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
  double mad_raw = robscale::median_select(dev, static_cast<size_t>(n));
  return constant * mad_raw;
}
