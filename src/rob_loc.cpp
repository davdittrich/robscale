#include "robscale_config.h"
#include "robust_core.h"
#include <Rcpp.h>

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

// [[Rcpp::export]]
double rob_loc_impl(Rcpp::NumericVector x, bool has_scale, double scale_val,
                    int maxit, double tol) {
  size_t n = (size_t)x.size();
  if (ROBSCALE_UNLIKELY(n == 0)) return 0.0;

  // Optimized Arena: dual stack-buffer allocation for n <= 2048
  constexpr size_t SCALE_STACK_SIZE = 2048;
  double buf_micro[128]; // 64 * 2
  double buf_stack[SCALE_STACK_SIZE * 2];
  std::unique_ptr<double[]> heap;
  double* arena;

  if (n <= 64) {
    arena = buf_micro;
  } else if (ROBSCALE_LIKELY(n <= SCALE_STACK_SIZE)) {
    arena = buf_stack;
  } else {
    heap.reset(new double[n * 2]);
    arena = heap.get();
  }
  
  double* buf = arena;
  double* dev = arena + n;
  const double* xp = x.begin();

  // Initial median calculation
  std::memcpy(buf, xp, n * sizeof(double));
  double med = robscale::median_select(buf, n);
  
  int minobs = has_scale ? 3 : 4;
  if (ROBSCALE_UNLIKELY(n < (size_t)minobs)) return med;

  // Scale estimation
  double s = has_scale ? scale_val : robscale::mad_select(xp, (int)n, med, dev);
  if (ROBSCALE_UNLIKELY(s == 0.0)) return med;

  return rob_loc_compute(xp, n, med, s, maxit, tol, buf);
}
