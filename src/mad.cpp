#include "robscale_config.h"
#include "pdq_select.h"
#include <Rcpp.h>
#include <cstring>
#include <memory>

// OPT-M1: Extract large-n path as ROBSCALE_NOINLINE so buf_stack[2048] is
// never allocated in the entry frame when n <= 64.
// OPT-M2/M3: ROBSCALE_RESTRICT on helpers + #pragma omp simd (via bulk_abs_diff).
// OPT-M4: bulk_abs_diff / bulk_abs_diff_inplace shared SIMD kernels.

static ROBSCALE_NOINLINE
double mad_impl_auto_large(const double* ROBSCALE_RESTRICT xp, int n, double constant) {
  constexpr int STACK_SIZE = 2048;
  double buf_stack[STACK_SIZE];
  std::unique_ptr<double[]> heap;
  double* w = (n <= STACK_SIZE)
    ? buf_stack
    : (heap.reset(new double[n]), heap.get());
  std::memcpy(w, xp, n * sizeof(double));
  double med = robscale::adaptive_median_select(w, static_cast<size_t>(n));
  robscale::bulk_abs_diff_inplace(w, n, med);
  double mad_raw = robscale::adaptive_median_select(w, static_cast<size_t>(n));
  return constant * mad_raw;
}

// MAD with auto-median: fused single-buffer approach.
// After median selection, w is a permutation of x. Compute deviations
// in-place: {|w[i] - med|} == {|x[i] - med|} as multisets, and MAD
// is order-invariant. This halves memory from 2n to n doubles.
// [[Rcpp::export]]
double mad_impl_auto(Rcpp::NumericVector x, double constant) {
  int n = x.size();
  if (n < 1) return NA_REAL;
  if (n == 1) return 0.0;
  if (n <= 64) {
    double buf_micro[64];
    const double* xp = x.begin();
    std::memcpy(buf_micro, xp, n * sizeof(double));
    double med = robscale::adaptive_median_select(buf_micro, static_cast<size_t>(n));
    robscale::bulk_abs_diff_inplace(buf_micro, n, med);
    return constant * robscale::adaptive_median_select(buf_micro, static_cast<size_t>(n));
  }
  return mad_impl_auto_large(x.begin(), n, constant);
}

static ROBSCALE_NOINLINE
double mad_impl_center_large(const double* ROBSCALE_RESTRICT xp, int n,
                              double center, double constant) {
  constexpr int STACK_SIZE = 2048;
  double buf_stack[STACK_SIZE];
  std::unique_ptr<double[]> heap;
  double* dev = (n <= STACK_SIZE)
    ? buf_stack
    : (heap.reset(new double[n]), heap.get());
  robscale::bulk_abs_diff(dev, xp, n, center);
  double mad_raw = robscale::adaptive_median_select(dev, static_cast<size_t>(n));
  return constant * mad_raw;
}

// MAD with user-supplied center
// [[Rcpp::export]]
double mad_impl_center(Rcpp::NumericVector x, double center, double constant) {
  int n = x.size();
  if (n < 1) return NA_REAL;
  if (n == 1) return 0.0;
  if (n <= 64) {
    double buf_micro[64];
    const double* xp = x.begin();
    robscale::bulk_abs_diff(buf_micro, xp, n, center);
    return constant * robscale::adaptive_median_select(buf_micro, static_cast<size_t>(n));
  }
  return mad_impl_center_large(x.begin(), n, center, constant);
}
