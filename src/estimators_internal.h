#ifndef ROBSCALE_ESTIMATORS_INTERNAL_H
#define ROBSCALE_ESTIMATORS_INTERNAL_H

#include "robscale_config.h"
#include "robust_core.h"
#include <cstring>
#include <cmath>

// Forward declarations for Qn/Sn internal implementations
// (defined in qn_estimator.cpp and sn_estimator.cpp)
namespace robscale::qnsn {
  template <typename T> double C_qn_impl(const T* x, size_t n);
  template <typename T> double C_sn_impl(const T* x, size_t n);
}

// rob_scale_compute: promoted from static in rob_scale.cpp
double rob_scale_compute(const double* ROBSCALE_RESTRICT data,
                         size_t n, double data_offset, double s,
                         int maxit, double tol,
                         double* ROBSCALE_RESTRICT tmp);

namespace robscale { namespace internal {

// GMD: sorts buf in-place, O(n log n)
inline double gmd(double* buf, int n) {
  if (n < 2) return 0.0;
  // Small sort via sorting networks
  if (n <= 16) {
    robscale::small_sort(buf, n);
  } else {
    std::sort(buf, buf + n);  // plain std::sort for ensemble's small resamples
  }
  double sum = 0.0;
  for (int i = 0; i < n; ++i)
    sum += (2.0 * (i + 1) - n - 1.0) * buf[i];
  return GMD_CONSISTENCY * 2.0 * sum / (static_cast<double>(n) * (n - 1));
}

// MAD: uses dev buffer for absolute deviations, O(n)
// buf is used as scratch for median selection (destructive)
// dev is used for absolute deviations
inline double mad(double* buf, double* dev, int n) {
  if (n < 2) return 0.0;
  double med = robscale::median_select(buf, n);
  // Need original data for deviations — but buf was destroyed by median_select.
  // Caller must provide dev[] already filled with |x_i - med|
  // ... Actually, let's take a different approach: accept original data pointer
  // This function won't be called directly — see mad_from_data below
  (void)dev; (void)med;
  return 0.0; // placeholder
}

// MAD from original data: computes median, then MAD
// Needs two buffers: buf (for median selection), dev (for deviations)
inline double mad_from_data(const double* x, double* buf, double* dev, int n) {
  if (n < 2) return 0.0;
  std::memcpy(buf, x, n * sizeof(double));
  double med = robscale::median_select(buf, static_cast<size_t>(n));
  for (int i = 0; i < n; ++i) dev[i] = std::abs(x[i] - med);
  return MAD_CONSISTENCY * robscale::median_select(dev, static_cast<size_t>(n));
}

// IQR: dual selection, O(n)
// Needs two buffers (selection is destructive)
inline double iqr(const double* x, double* buf1, double* buf2, int n) {
  if (n < 2) return 0.0;

  // Q1
  std::memcpy(buf1, x, n * sizeof(double));
  double h1 = (n - 1.0) * 0.25;
  int lo1 = static_cast<int>(h1);
  double frac1 = h1 - lo1;
  robscale::floyd_rivest_select(buf1, buf1 + lo1, buf1 + n);
  double q1 = buf1[lo1];
  if (frac1 > 0.0 && lo1 + 1 < n) {
    double next_val = buf1[lo1 + 1];
    for (int i = lo1 + 2; i < n; ++i)
      if (buf1[i] < next_val) next_val = buf1[i];
    q1 += frac1 * (next_val - q1);
  }

  // Q3
  std::memcpy(buf2, x, n * sizeof(double));
  double h3 = (n - 1.0) * 0.75;
  int lo3 = static_cast<int>(h3);
  double frac3 = h3 - lo3;
  robscale::floyd_rivest_select(buf2, buf2 + lo3, buf2 + n);
  double q3 = buf2[lo3];
  if (frac3 > 0.0 && lo3 + 1 < n) {
    double next_val = buf2[lo3 + 1];
    for (int i = lo3 + 2; i < n; ++i)
      if (buf2[i] < next_val) next_val = buf2[i];
    q3 += frac3 * (next_val - q3);
  }

  return (q3 - q1) * IQR_CONSISTENCY;
}

// SD with c4 correction: read-only, O(n)
inline double sd_c4(const double* x, int n) {
  if (n < 2) return 0.0;
  // Welford's online algorithm
  double mean = x[0];
  double m2 = 0.0;
  for (int i = 1; i < n; ++i) {
    double delta = x[i] - mean;
    mean += delta / (i + 1.0);
    double delta2 = x[i] - mean;
    m2 += delta * delta2;
  }
  double sd = std::sqrt(m2 / (n - 1.0));
  return sd / robscale::c4_factor(n);
}

// Sn: delegates to existing optimized implementation
inline double sn(const double* x, int n) {
  if (n < 2) return 0.0;
  return robscale::qnsn::C_sn_impl<double>(x, static_cast<size_t>(n));
}

// Qn: delegates to existing optimized implementation
inline double qn(const double* x, int n) {
  if (n < 2) return 0.0;
  return robscale::qnsn::C_qn_impl<double>(x, static_cast<size_t>(n));
}

// robScale (M-scale): Newton-Raphson iteration
// Needs two buffers: buf (for median/scratch), dev (for MAD deviations + tanh scratch)
inline double rob_scale(const double* x, double* buf, double* dev, int n) {
  if (n < 4) return 0.0; // minimum for robScale without known location

  // Compute median
  std::memcpy(buf, x, n * sizeof(double));
  double t = robscale::median_select(buf, static_cast<size_t>(n));

  // Compute MAD as initial scale
  for (int i = 0; i < n; ++i) dev[i] = std::abs(x[i] - t);
  double s_init = MAD_CONSISTENCY * robscale::median_select(dev, static_cast<size_t>(n));

  // MAD implosion: return ADM fallback
  if (s_init <= 1e-4) {
    return robscale::adm_core(x, n, t, ADM_CONSISTENCY);
  }

  // Newton-Raphson iteration via promoted rob_scale_compute
  return rob_scale_compute(x, static_cast<size_t>(n), t, s_init, 80, 1.4901161e-8, buf);
}

}} // namespace robscale::internal

#endif // ROBSCALE_ESTIMATORS_INTERNAL_H
