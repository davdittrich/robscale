#ifndef ROBSCALE_VALIDATE_FINITE_H
#define ROBSCALE_VALIDATE_FINITE_H

#include "robust_core.h"
#include <Rcpp.h>

// Shared input validation: scans for first non-finite element and throws
// Rcpp::stop with the appropriate message (NA vs Inf). No allocation,
// short-circuits on first bad value. Used by all Rcpp-exported entry points.
inline void validate_finite(const double* xp, int n) {
  int idx = robscale::find_first_nonfinite(xp, n);
  if (ROBSCALE_UNLIKELY(idx >= 0)) {
    if (std::isnan(xp[idx]))
      Rcpp::stop("There are NAs in the data yet na.rm is FALSE");
    else
      Rcpp::stop("'x' must not contain non-finite values (Inf, -Inf, NaN)");
  }
}

#endif // ROBSCALE_VALIDATE_FINITE_H
