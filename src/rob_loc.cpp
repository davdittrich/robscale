#include "robscale_config.h"
#include "robust_core.h"
#include <Rcpp.h>
#include <memory>

/**
 * Portably optimized robLoc kernel.
 */
static ROBSCALE_INLINE double rob_loc_compute(const double* ROBSCALE_RESTRICT xp,
                                              size_t n, double t, double s,
                                              int maxit, double tol,
                                              double* ROBSCALE_RESTRICT tmp) {
  double inv_s = 1.0 / s;
  double half_inv_s = 0.5 * inv_s;

  for (int k = 0; k < maxit; ++k) {
    for (size_t i = 0; i < n; ++i)
      tmp[i] = (xp[i] - t) * half_inv_s;

    robscale::bulk_tanh(tmp, (int)n);

    double sum_psi = 0.0, sum_dpsi = 0.0;
    for (size_t i = 0; i < n; ++i) {
      double p = tmp[i];
      sum_psi += p;
      sum_dpsi += 1.0 - p * p;
    }
    double v = 2.0 * s * sum_psi / sum_dpsi;
    t += v;
    if (std::abs(v) <= tol) break;
  }
  return t;
}

/**
 * Shared core: median, MAD, iteration.
 * Called by both the small-n and large-n entry points.
 */
static double rob_loc_core(const double* xp, size_t n,
                           double* buf, double* dev,
                           bool has_scale, double scale_val,
                           int maxit, double tol) {
  std::memcpy(buf, xp, n * sizeof(double));
  double med = robscale::median_select(buf, n);

  int minobs = has_scale ? 3 : 4;
  if (ROBSCALE_UNLIKELY(n < (size_t)minobs)) return med;

  double s = has_scale ? scale_val : robscale::mad_select(xp, (int)n, med, dev);
  if (ROBSCALE_UNLIKELY(s == 0.0)) return med;

  return rob_loc_compute(xp, n, med, s, maxit, tol, buf);
}

/**
 * Small-n entry point (n <= 64): minimal stack frame (~1KB).
 */
ROBSCALE_NOINLINE
static double rob_loc_impl_small(const double* xp, size_t n,
                                 bool has_scale, double scale_val,
                                 int maxit, double tol) {
  double arena[ROBSCALE_MICRO_BUFFER_SIZE]; // 128 doubles = 1KB
  return rob_loc_core(xp, n, arena, arena + n,
                      has_scale, scale_val, maxit, tol);
}

// [[Rcpp::export]]
double rob_loc_impl(Rcpp::NumericVector x, bool has_scale, double scale_val,
                    int maxit, double tol) {
  size_t n = (size_t)x.size();
  if (ROBSCALE_UNLIKELY(n == 0)) return 0.0;

  const double* xp = x.begin();

  // Small-n: dispatch to noinline helper with minimal stack frame
  if (n <= 64)
    return rob_loc_impl_small(xp, n, has_scale, scale_val, maxit, tol);

  // Large-n: stack or heap arena
  constexpr size_t SCALE_STACK_SIZE = 2048;
  double buf_stack[SCALE_STACK_SIZE * 2];
  std::unique_ptr<double[]> heap;
  double* arena;

  if (ROBSCALE_LIKELY(n <= SCALE_STACK_SIZE)) {
    arena = buf_stack;
  } else {
    heap.reset(new double[n * 2]);
    arena = heap.get();
  }

  return rob_loc_core(xp, n, arena, arena + n,
                      has_scale, scale_val, maxit, tol);
}
