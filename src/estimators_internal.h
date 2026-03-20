#ifndef ROBSCALE_ESTIMATORS_INTERNAL_H
#define ROBSCALE_ESTIMATORS_INTERNAL_H

#include "robscale_config.h"
#include "robust_core.h"
#include "pdq_select.h"
#include <cstring>
#include <cmath>

// Forward declarations for Qn/Sn internal implementations
// (defined in qn_estimator.cpp and sn_estimator.cpp)
namespace robscale::qnsn {
  template <typename T> double C_qn_impl(const T* x, size_t n);
  template <typename T> double C_sn_impl(const T* x, size_t n);
  template <typename T> double C_qn_impl_sorted(const T* x, size_t n);
  template <typename T> double C_sn_impl_sorted(const T* x, size_t n);
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
#if defined(_OPENMP) || defined(ROBSCALE_HAS_OMP_SIMD)
  #pragma omp simd reduction(+:sum)
#endif
  for (int i = 0; i < n; ++i)
    sum += (2.0 * (i + 1) - n - 1.0) * buf[i];
  return GMD_CONSISTENCY * 2.0 * sum / (static_cast<double>(n) * (n - 1));
}

// MAD from original data: fused single-buffer approach.
// After median selection, buf is a permutation of x. Deviations computed
// in-place produce the same multiset as |x[i] - med|.
inline double mad_from_data(const double* x, double* buf, int n) {
  if (n < 2) return 0.0;
  std::memcpy(buf, x, n * sizeof(double));
  double med = robscale::adaptive_median_select(buf, static_cast<size_t>(n));
  for (int i = 0; i < n; ++i) buf[i] = std::abs(buf[i] - med);
  return MAD_CONSISTENCY * robscale::adaptive_median_select(buf, static_cast<size_t>(n));
}

// IQR: incremental pdqselect, O(n)
// Uses buf1 only; buf2 is unused (signature kept for call-site compatibility)
inline double iqr(const double* x, double* buf1, double* buf2, int n) {
  (void)buf2;
  if (n < 2) return 0.0;

  std::memcpy(buf1, x, n * sizeof(double));

  double h1 = (n - 1.0) * 0.25;
  int lo1 = static_cast<int>(h1);
  double frac1 = h1 - lo1;

  double h3 = (n - 1.0) * 0.75;
  int lo3 = static_cast<int>(h3);
  double frac3 = h3 - lo3;

  // Q1: select on full array
  miniselect::pdqselect(buf1, buf1 + lo1, buf1 + n);
  double q1 = robscale::interp_q7(buf1, n, lo1, frac1);

  // Q3: select on buf1[lo1+1 .. n-1] (everything <= Q1 is irrelevant)
  int start = lo1 + 1;
  miniselect::pdqselect(buf1 + start, buf1 + lo3, buf1 + n);
  double q3 = robscale::interp_q7(buf1, n, lo3, frac3);

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

// Sorted variants: input MUST be sorted ascending. No copy, no sort.
inline double sn_sorted(const double* sorted_x, int n) {
  if (n < 2) return 0.0;
  return robscale::qnsn::C_sn_impl_sorted<double>(sorted_x, static_cast<size_t>(n));
}

inline double qn_sorted(const double* sorted_x, int n) {
  if (n < 2) return 0.0;
  return robscale::qnsn::C_qn_impl_sorted<double>(sorted_x, static_cast<size_t>(n));
}

// robScale (M-scale): fused single-buffer approach.
// Median → in-place deviations → MAD → iteration, all on one buffer.
inline double rob_scale(const double* x, double* buf, int n) {
  if (n < 4) return 0.0; // minimum for robScale without known location

  // Compute median (destroys ordering of buf)
  std::memcpy(buf, x, n * sizeof(double));
  double t = robscale::median_select(buf, static_cast<size_t>(n));

  // Compute MAD in-place: deviations on the permuted buf
  for (int i = 0; i < n; ++i) buf[i] = std::abs(buf[i] - t);
  double s_init = MAD_CONSISTENCY * robscale::median_select(buf, static_cast<size_t>(n));

  // MAD implosion: return ADM fallback
  if (s_init <= robscale::IMPLOSION_BOUND) {
    return robscale::adm_core(x, n, t, ADM_CONSISTENCY);
  }

  // Newton-Raphson iteration: buf is reused as scratch (written before read)
  return rob_scale_compute(x, static_cast<size_t>(n), t, s_init, 80, 1.4901161e-8, buf);
}

}} // namespace robscale::internal

#endif // ROBSCALE_ESTIMATORS_INTERNAL_H
