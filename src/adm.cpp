#include "robust_core.h"
#include "pdq_select.h"   // for adaptive_median_select — same pattern as mad.cpp:2
#include <Rcpp.h>
#include <memory>

// [[Rcpp::export]]
double adm_impl(Rcpp::NumericVector x, double center, double constant) {
  return robscale::adm_core(x.begin(), (int)x.size(), center, constant);
}

// Large-n helper: allocates its own stack/heap buffer so that adm_impl_auto's
// small-n frame stays lean (buf_micro[128] only).
// ROBSCALE_NOINLINE ensures the compiler does not merge the two stack frames.
// ROBSCALE_RESTRICT: eliminates aliasing so the compiler can optimise the copy.
// Stack buffer 2048: matches mad.cpp; covers the majority of ensemble call sites.
static ROBSCALE_NOINLINE
double adm_large_n(const double* ROBSCALE_RESTRICT xp, int n, double constant) {
  constexpr int STACK_BUF_SIZE = 2048;
  double buf_stack[STACK_BUF_SIZE];
  std::unique_ptr<double[]> heap;
  double* buf;
  if (n <= STACK_BUF_SIZE) {
    buf = buf_stack;
  } else {
    heap.reset(new double[n]);
    buf = heap.get();
  }
  std::memcpy(buf, xp, n * sizeof(double));
  double med = robscale::adaptive_median_select(buf, static_cast<size_t>(n));
  return robscale::adm_core(buf, n, med, constant);
}

// [[Rcpp::export]]
double adm_impl_auto(Rcpp::NumericVector x, double constant) {
  int n = x.size();
  if (ROBSCALE_UNLIKELY(n == 0)) return 0.0;
  if (ROBSCALE_UNLIKELY(n == 1)) return 0.0;
  if (n <= ROBSCALE_MICRO_BUFFER_SIZE) {
    // Only buf_micro lives in this frame — keeps the hot path stack-lean.
    double buf_micro[ROBSCALE_MICRO_BUFFER_SIZE];
    std::memcpy(buf_micro, x.begin(), n * sizeof(double));
    double med = robscale::adaptive_median_select(buf_micro, static_cast<size_t>(n));
    return robscale::adm_core(buf_micro, n, med, constant);
  }
  return adm_large_n(x.begin(), n, constant);
}
