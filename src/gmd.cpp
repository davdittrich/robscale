#include "robscale_config.h"
#include "robust_core.h"
#include "qnsn_sort_utils.h"
#include <Rcpp.h>
#include <cstring>
#include <memory>

// GMD sorted kernel: assumes x is already sorted ascending
static ROBSCALE_INLINE double gmd_sorted(const double* x, int n, double constant) {
  double sum = 0.0;
  for (int i = 0; i < n; ++i)
    sum += (2.0 * (i + 1) - n - 1.0) * x[i];
  return constant * 2.0 * sum / (static_cast<double>(n) * (n - 1));
}

// [[Rcpp::export]]
double gmd_impl(Rcpp::NumericVector x, double constant) {
  int n = x.size();
  if (n < 2) return 0.0;

  // Arena allocation
  double buf_micro[ROBSCALE_MICRO_BUFFER_SIZE];
  constexpr int STACK_SIZE = 2048;
  double buf_stack[STACK_SIZE];
  std::unique_ptr<double[]> heap;
  double* buf;

  if (n <= ROBSCALE_MICRO_BUFFER_SIZE) {
    buf = buf_micro;
  } else if (n <= STACK_SIZE) {
    buf = buf_stack;
  } else {
    heap.reset(new double[n]);
    buf = heap.get();
  }

  std::memcpy(buf, x.begin(), n * sizeof(double));

  // Tiered sort: sorting networks for n<=16, then optimized_sort
  if (n <= 16) {
    robscale::small_sort(buf, n);
  } else {
    robscale::qnsn::optimized_sort(buf, buf + n);
  }

  return gmd_sorted(buf, n, constant);
}
