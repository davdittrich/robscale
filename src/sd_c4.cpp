#include "robscale_config.h"
#include "robust_core.h"
#include <Rcpp.h>
#include <cmath>

// Input validation helper: stops with a clear error on first non-finite value.
static void validate_finite(const double* xp, int n) {
  int idx = robscale::find_first_nonfinite(xp, n);
  if (ROBSCALE_UNLIKELY(idx >= 0)) {
    if (std::isnan(xp[idx]))
      Rcpp::stop("There are NAs in the data yet na.rm is FALSE");
    else
      Rcpp::stop("'x' must not contain non-finite values (Inf, -Inf, NaN)");
  }
}

// [[Rcpp::export(rng = false)]]
double sd_c4_impl(Rcpp::NumericVector x) {
  int n = x.size();
  if (n < 2) return NA_REAL;

  const double* xp = x.begin();
  validate_finite(xp, n);

  // Welford's online algorithm for numerically stable variance
  double mean = xp[0];
  double m2 = 0.0;
  for (int i = 1; i < n; ++i) {
    double delta = xp[i] - mean;
    mean += delta / (i + 1.0);
    double delta2 = xp[i] - mean;
    m2 += delta * delta2;
  }

  double sd = std::sqrt(m2 / (n - 1.0));
  return sd / robscale::c4_factor(n);
}
