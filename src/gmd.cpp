#include "robscale_config.h"
#include "robust_core.h"
#include "qnsn_sort_utils.h"
#include <Rcpp.h>
#include <cstring>
#include <memory>

// GMD sorted kernel: assumes x is already sorted ascending.
// OPT-G6: ROBSCALE_RESTRICT allows the compiler to assume no aliasing.
// OPT-G5: caller precomputes scale = constant * 2 / (n*(n-1)) once.
static ROBSCALE_INLINE double gmd_sorted(const double* ROBSCALE_RESTRICT x,
                                          int n, double scale) {
  double sum = 0.0;
#if defined(_OPENMP) || defined(ROBSCALE_HAS_OMP_SIMD)
  #pragma omp simd reduction(+:sum)
#endif
  for (int i = 0; i < n; ++i)
    sum += (2.0 * (i + 1) - n - 1.0) * x[i];
  return scale * sum;
}

// Shared sort+kernel logic.
// OPT-G5: computes scale once before the loop.
// OPT-G6: ROBSCALE_RESTRICT on both pointer parameters.
static double gmd_core(const double* ROBSCALE_RESTRICT xp, int n,
                       double* ROBSCALE_RESTRICT buf, double constant) {
  const double scale = constant * 2.0 / (static_cast<double>(n) * (n - 1));
  std::memcpy(buf, xp, n * sizeof(double));
  if (n <= 16)
    robscale::small_sort(buf, n);
  else
    robscale::qnsn::optimized_sort(buf, buf + n);
  return gmd_sorted(buf, n, scale);
}

// OPT-G1: small-n path with minimal stack frame (1 KB arena only).
// ROBSCALE_NOINLINE isolates this frame so buf_stack[2048] is never
// allocated when n <= ROBSCALE_MICRO_BUFFER_SIZE.
ROBSCALE_NOINLINE
static double gmd_impl_small(const double* xp, int n, double constant) {
  double arena[ROBSCALE_MICRO_BUFFER_SIZE];  // 128 doubles = 1 KB
  return gmd_core(xp, n, arena, constant);
}

// [[Rcpp::export]]
double gmd_impl(Rcpp::NumericVector x, double constant) {
  int n = x.size();
  if (n < 2) return 0.0;
  if (n <= ROBSCALE_MICRO_BUFFER_SIZE)
    return gmd_impl_small(x.begin(), n, constant);
  constexpr int STACK_SIZE = ROBSCALE_SN_STACK_THRESHOLD;
  double buf_stack[STACK_SIZE];
  std::unique_ptr<double[]> heap;
  double* buf = (n <= STACK_SIZE)
    ? buf_stack
    : (heap.reset(new double[n]), heap.get());
  return gmd_core(x.begin(), n, buf, constant);
}

