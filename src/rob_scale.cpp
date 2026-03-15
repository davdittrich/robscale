#include "robscale_config.h"
#include "robust_core.h"
#include "pdq_select.h"
#include <Rcpp.h>

/**
 * Portably optimized robScale kernel.
 * Non-static: also called by estimators_internal.h for the ensemble.
 */
double rob_scale_compute(const double* ROBSCALE_RESTRICT data,
                         size_t n, double data_offset, double s,
                         int maxit, double tol,
                         double* ROBSCALE_RESTRICT tmp) {
  double inv_n = 1.0 / (double)n;

  for (int k = 0; k < maxit; ++k) {
    double half_inv_sc = 0.5 * robscale::INV_RHO_SCALE_CONST / s;
    for (size_t i = 0; i < n; ++i)
      tmp[i] = (data[i] - data_offset) * half_inv_sc;
    
    robscale::bulk_tanh(tmp, (int)n);

    double sum_rho = 0.0;
    for (size_t i = 0; i < n; ++i) {
      sum_rho += tmp[i] * tmp[i];
    }
    double v = std::sqrt(2.0 * sum_rho * inv_n);
    s *= v;
    if (std::abs(v - 1.0) <= tol) break;
  }
  return s;
}

// [[Rcpp::export]]
double rob_scale_impl(Rcpp::NumericVector x, bool has_loc, double loc_val,
                      double implbound, int maxit, double tol, int fallback) {
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
  
  double* w = arena;     // Used for median/centering
  double* dev = arena + n; // Used for MAD deviations
  const double* xp = x.begin();

  double t, s_init;
  if (has_loc) {
    t = loc_val;
    // Calculate MAD around fixed loc
    for (size_t i = 0; i < n; ++i) dev[i] = std::abs(xp[i] - t);
    s_init = robscale::MAD_CONSISTENCY * robscale::adaptive_robscale_median_select(dev, n);
  } else {
    std::memcpy(w, xp, n * sizeof(double));
    t = robscale::adaptive_robscale_median_select(w, n);
    s_init = robscale::adaptive_mad_select(xp, (int)n, t, dev);
  }

  // Small-sample fallback (matches original revss behavior)
  int minobs = has_loc ? 3 : 4;
  if (ROBSCALE_UNLIKELY(n < (size_t)minobs)) {
    if (s_init <= implbound) {
      if (ROBSCALE_UNLIKELY(fallback == 1)) return R_NaReal;
      if (has_loc) {
        // Re-check MAD around the data's own median (not the supplied loc)
        std::memcpy(w, xp, n * sizeof(double));
        double med_orig = robscale::adaptive_robscale_median_select(w, n);
        double mad_orig = robscale::adaptive_mad_select(xp, (int)n, med_orig, dev);
        return (mad_orig <= implbound)
          ? robscale::adm_core(xp, (int)n, med_orig, robscale::ADM_CONSISTENCY)
          : mad_orig;
      } else {
        return robscale::adm_core(xp, (int)n, t, robscale::ADM_CONSISTENCY);
      }
    }
    return s_init;
  }

  // MAD collapse for n >= minobs
  if (ROBSCALE_UNLIKELY(s_init <= implbound && fallback == 1)) return R_NaReal;
  if (ROBSCALE_UNLIKELY(s_init == 0.0)) {
    return robscale::adm_core(xp, (int)n, t, robscale::ADM_CONSISTENCY);
  }

  const double* data = xp; 
  double data_offset = t;
  
  return rob_scale_compute(data, n, data_offset, s_init, maxit, tol, w);
}
