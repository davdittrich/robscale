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
  if (n <= ROBSCALE_MEDIAN_NET_THRESHOLD) {
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

// Templated adaptive median: FR below threshold, pdqselect above.
// ThreshField selects which RuntimeConfig threshold to compare against.
// Template instantiation makes the member pointer a compile-time constant —
// no runtime overhead vs a plain conditional on the field directly.
template <size_t robscale::qnsn::RuntimeConfig::*ThreshField>
ROBSCALE_INLINE double adaptive_median_select_t(double* x, size_t n) {
  if (n <= robscale::qnsn::RuntimeConfig::get().*ThreshField)
    return median_select(x, n);
  return pdq_median_select(x, n);
}

// Adaptive median for IQR/MAD (pdq_median_threshold).
ROBSCALE_INLINE double adaptive_median_select(double* x, size_t n) {
  return adaptive_median_select_t<&robscale::qnsn::RuntimeConfig::pdq_median_threshold>(x, n);
}

// Adaptive median for robScale (pdq_robscale_threshold).
// robScale working set: 1–2 warm arrays → lighter cache pressure than MAD.
ROBSCALE_INLINE double adaptive_robscale_median_select(double* x, size_t n) {
  return adaptive_median_select_t<&robscale::qnsn::RuntimeConfig::pdq_robscale_threshold>(x, n);
}

// Adaptive MAD for robScale: fill deviations then select median adaptively.
ROBSCALE_INLINE double adaptive_mad_select(const double* x, int n, double med, double* dev) {
  for (int i = 0; i < n; ++i) dev[i] = std::abs(x[i] - med);
  return robscale::MAD_CONSISTENCY *
         adaptive_robscale_median_select(dev, static_cast<size_t>(n));
}

// Adaptive low-median (no even-n averaging) for Sn inner_medians.
// Templated to support both double and float inner_medians arrays.
template <typename T>
ROBSCALE_INLINE double adaptive_lowmedian_select(T* x, size_t n) {
  if (ROBSCALE_UNLIKELY(n == 0)) return 0.0;
  size_t h = (n - 1) / 2;
  if (n <= robscale::qnsn::RuntimeConfig::get().pdq_lowmedian_threshold)
    robscale::floyd_rivest_select(x, x + h, x + n);
  else
    miniselect::pdqselect(x, x + h, x + n);
  return static_cast<double>(x[h]);
}

// Type 7 quantile interpolation after selection.
// Precondition: buf[lo] is the lo-th order statistic (via pdqselect or similar).
// If frac > 0, scans buf[lo+1..n-1] for the next order statistic.
inline double interp_q7(double* ROBSCALE_RESTRICT buf, int n, int lo, double frac) {
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
