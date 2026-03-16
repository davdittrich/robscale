#include "robscale_config.h"
#include <Rcpp.h>
#include <cmath>

// [[Rcpp::export]]
double sd_c4_impl(Rcpp::NumericVector x) {
  int n = x.size();
  if (n < 2) return NA_REAL;

  const double* xp = x.begin();

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
