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

// ------------------------------------------------------------
// WU-ADM3: diagnostic export for adm_core_sorted benchmarking
// ------------------------------------------------------------

// C_adm_core_sorted: benchmarks adm_core_sorted vs adm_core on sorted input.
// Caller MUST pass a sorted vector; center is derived as median_sorted.
// Remove after WU-ADM3 validation is complete (cleanup in WU-ADM4).
// [[Rcpp::export]]
double C_adm_core_sorted(Rcpp::NumericVector x) {
  int n = x.size();
  if (n == 0) return 0.0;
  double med = robscale::median_sorted(x.begin(), static_cast<size_t>(n));
  return robscale::adm_core_sorted(x.begin(), n, med, robscale::ADM_CONSISTENCY);
}
