#include "robust_core.h"
#include <Rcpp.h>
#include <memory>
#ifndef STACK_BUF_SIZE
#define STACK_BUF_SIZE 1024
#endif

// [[Rcpp::export]]
double adm_impl(Rcpp::NumericVector x, double center, double constant) {
  return robscale::adm_core(x.begin(), (int)x.size(), center, constant);
}

// Large-n helper: allocates its own stack/heap buffer so that adm_impl_auto's
// small-n frame stays lean (buf_micro[64] only, ~512 bytes vs 8256 bytes before).
// ROBSCALE_NOINLINE ensures the compiler does not merge the two stack frames.
// OPT-A: adm_core is called with the hot copy buffer, not the cold SEXP pointer.
static ROBSCALE_NOINLINE
double adm_large_n(const double* xp, int n, double constant) {
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
  double med = robscale::median_select(buf, n);
  return robscale::adm_core(buf, n, med, constant);
}

// [[Rcpp::export]]
double adm_impl_auto(Rcpp::NumericVector x, double constant) {
  int n = x.size();
  if (n == 0) return 0.0;
  if (n <= 64) {
    // OPT-B: only buf_micro lives in this frame (~512 bytes).
    // OPT-A: pass hot buf_micro to adm_core instead of the cold x.begin() pointer.
    double buf_micro[64];
    std::memcpy(buf_micro, x.begin(), n * sizeof(double));
    double med = robscale::median_select(buf_micro, n);
    return robscale::adm_core(buf_micro, n, med, constant);
  }
  return adm_large_n(x.begin(), n, constant);
}
