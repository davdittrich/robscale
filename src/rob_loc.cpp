#include <Rcpp.h>
#include "robust_core.h"

// [[Rcpp::export]]
double rob_loc_impl(Rcpp::NumericVector x, bool has_scale, double scale_val,
                    int maxit, double tol) {
  int n = x.size();
  int minobs = has_scale ? 3 : 4;

  // Arena: [0..n) for x copy (median selection / reused as tmp), [n..2n) for MAD deviations
  double buf_stack[STACK_BUF_SIZE * 2];
  double* arena = (n * 2 <= STACK_BUF_SIZE * 2) ? buf_stack : new double[n * 2];
  double* buf = arena;
  double* dev = arena + n;

  std::memcpy(buf, x.begin(), n * sizeof(double));
  double med = median_select(buf, n);

  if (n < minobs) {
    if (n * 2 > STACK_BUF_SIZE * 2) delete[] arena;
    return med;
  }

  double s = has_scale ? scale_val : mad_select(x.begin(), n, med, dev);

  if (s == 0.0) {
    if (n * 2 > STACK_BUF_SIZE * 2) delete[] arena;
    return med;
  }

  // True Newton-Raphson iteration with bulk_tanh
  double t = med;
  const double* xp = x.begin();
  double inv_s = 1.0 / s;
  double* tmp = buf;

  for (int k = 0; k < maxit; ++k) {
    double half_inv_s = 0.5 * inv_s;
    for (int i = 0; i < n; ++i)
      tmp[i] = (xp[i] - t) * half_inv_s;
    bulk_tanh(tmp, n);

    double sum_psi = 0.0, sum_dpsi = 0.0;
    for (int i = 0; i < n; ++i) {
      double p = tmp[i];
      sum_psi += p;
      sum_dpsi += 1.0 - p * p;
    }
    double v = 2.0 * s * sum_psi / sum_dpsi;
    t += v;
    if (std::abs(v) <= tol) break;
  }

  if (n * 2 > STACK_BUF_SIZE * 2) delete[] arena;
  return t;
}
