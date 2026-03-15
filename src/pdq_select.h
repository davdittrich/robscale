#ifndef ROBSCALE_PDQ_SELECT_H
#define ROBSCALE_PDQ_SELECT_H

#include "robscale_config.h"
#include "robust_core.h"
#include "qnsn_runtime_config.h"
#include "sort_net.h"
#include "miniselect/pdqselect.h"

namespace robscale {

// Median via pdqselect: sorting networks for n <= 16, pdqselect for n > 16.
// Mirrors median_select() semantics but uses pdqselect instead of Floyd-Rivest.
ROBSCALE_INLINE double pdq_median_select(double* x, size_t n) {
  if (ROBSCALE_UNLIKELY(n == 0)) return 0.0;
  if (n <= 16) {
    robscale::small_sort(x, n);
    if (n & 1) return x[n / 2];
    return (x[(n / 2) - 1] + x[n / 2]) * 0.5;
  }
  size_t h = (n - 1) / 2;
  miniselect::pdqselect(x, x + h, x + n);
  if (n & 1) return x[h];

  // Even n: scan for the next order statistic above h
  double v1 = x[h];
  double v2 = x[h + 1];
  for (size_t i = h + 2; i < n; ++i) {
    if (x[i] < v2) v2 = x[i];
  }
  return (v1 + v2) * 0.5;
}

// Adaptive median selection: FR for medium n, pdqselect for large n.
// Threshold is derived from per-core L2 cache size at startup.
ROBSCALE_INLINE double adaptive_median_select(double* x, size_t n) {
  if (n <= robscale::qnsn::RuntimeConfig::get().pdq_median_threshold)
    return median_select(x, n);
  return pdq_median_select(x, n);
}

// Type 7 quantile interpolation after selection.
// Precondition: buf[lo] is the lo-th order statistic (via pdqselect or similar).
// If frac > 0, scans buf[lo+1..n-1] for the next order statistic.
inline double interp_q7(double* buf, int n, int lo, double frac) {
  double q = buf[lo];
  if (frac > 0.0 && lo + 1 < n) {
    double nv = buf[lo + 1];
    for (int i = lo + 2; i < n; ++i)
      if (buf[i] < nv) nv = buf[i];
    q += frac * (nv - q);
  }
  return q;
}

} // namespace robscale

#endif // ROBSCALE_PDQ_SELECT_H
