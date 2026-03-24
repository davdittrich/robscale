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

// ------------------------------------------------------------
// WU-ADM0: diagnostic exports for H2H benchmarking
// ------------------------------------------------------------

// C_adm_orig: frozen baseline — self-contained, calls no external symbol
// that any later WU modifies.  The adm_core_frozen lambda inlines the
// original single-accumulator loop so this function remains immune to
// changes made to robscale::adm_core in WU-ADM1 onwards.
// [[Rcpp::export]]
double C_adm_orig(Rcpp::NumericVector x) {
  int n = x.size();
  if (n == 0) return 0.0;
  const double* xp = x.begin();
  auto adm_core_frozen = [](const double* px, int pn, double center,
                             double constant) -> double {
    double sum = 0.0;
    for (int i = 0; i < pn; ++i) sum += std::abs(px[i] - center);
    return constant * sum / pn;
  };
  if (n <= 64) {
    double buf[64];
    std::memcpy(buf, xp, n * sizeof(double));
    double med = robscale::median_select(buf, n);
    return adm_core_frozen(buf, n, med, robscale::ADM_CONSISTENCY);
  }
  constexpr int FROZEN_STACK_SIZE = 1024;
  double buf_stack[FROZEN_STACK_SIZE];
  std::unique_ptr<double[]> heap;
  double* buf = (n <= FROZEN_STACK_SIZE)
    ? buf_stack
    : (heap.reset(new double[n]), heap.get());
  std::memcpy(buf, xp, n * sizeof(double));
  double med = robscale::median_select(buf, n);
  return adm_core_frozen(buf, n, med, robscale::ADM_CONSISTENCY);
}

// C_adm_fast: thin wrapper — always calls the current adm_impl_auto so
// it automatically picks up any kernel improvements from later WUs.
// [[Rcpp::export]]
double C_adm_fast(Rcpp::NumericVector x) {
  return adm_impl_auto(x, robscale::ADM_CONSISTENCY);
}
