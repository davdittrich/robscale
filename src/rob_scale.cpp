#pragma GCC optimize("O3")
#include <Rcpp.h>
#include "robust_core.h"

// [[Rcpp::export]]
double rob_scale_impl(Rcpp::NumericVector x, bool has_loc, double loc_val,
                      double implbound, int maxit, double tol, int fallback) {
  int n = x.size();
  const double* xp = x.begin();

  // Arena: [0..n) working copy, [n..2n) abs values / reused as tmp, [2n..3n) extra
  double buf_stack[STACK_BUF_SIZE * 3];
  std::unique_ptr<double[]> heap;
  double* arena;
  if (ROBSCALE_LIKELY(n * 3 <= STACK_BUF_SIZE * 3)) {
    arena = buf_stack;
  } else {
    heap.reset(new double[n * 3]);
    arena = heap.get();
  }
  double* w = arena;
  double* abs_buf = arena + n;
  double* extra = arena + 2 * n;

  double t, s;
  int minobs;

  if (has_loc) {
    for (int i = 0; i < n; ++i) w[i] = xp[i] - loc_val;
    for (int i = 0; i < n; ++i) abs_buf[i] = std::abs(w[i]);
    s = MAD_CONSISTENCY * median_select(abs_buf, n);
    t = 0.0;
    minobs = 3;
  } else {
    std::memcpy(w, xp, n * sizeof(double));
    double med = median_select(w, n);
    s = mad_select(xp, n, med, abs_buf);
    t = med;
    minobs = 4;
  }

  // Small-sample fallback (matches revss behavior)
  if (ROBSCALE_UNLIKELY(n < minobs)) {
    if (s <= implbound) {
      if (ROBSCALE_UNLIKELY(fallback == 1)) return NA_REAL; // "na" fallback
      if (has_loc) {
        std::memcpy(extra, xp, n * sizeof(double));
        double med_orig = median_select(extra, n);
        double mad_orig = mad_select(xp, n, med_orig, abs_buf);
        return (mad_orig <= implbound)
          ? adm_core(xp, n, med_orig, ADM_CONSISTENCY)
          : mad_orig;
      } else {
        return adm_core(xp, n, t, ADM_CONSISTENCY);
      }
    } else {
      return s;
    }
  }

  // MAD collapse for n >= minobs
  if (ROBSCALE_UNLIKELY(s <= implbound && fallback == 1)) return NA_REAL;
  if (ROBSCALE_UNLIKELY(s == 0.0)) return adm_core(xp, n, t, ADM_CONSISTENCY);

  // Multiplicative iteration with bulk_tanh
  const double* data = has_loc ? w : xp;
  double data_offset = has_loc ? 0.0 : t;
  double inv_n = 1.0 / n;
  double* tmp = abs_buf;

  for (int k = 0; k < maxit; ++k) {
    double half_inv_sc = 0.5 * INV_RHO_SCALE_CONST / s;
    for (int i = 0; i < n; ++i)
      tmp[i] = (data[i] - data_offset) * half_inv_sc;
    bulk_tanh(tmp, n);

    double sum_rho = 0.0;
    for (int i = 0; i < n; ++i) {
      sum_rho += tmp[i] * tmp[i];
    }
    double v = std::sqrt(2.0 * sum_rho * inv_n);
    s *= v;
    if (std::abs(v - 1.0) <= tol) break;
  }

  return s;
}
