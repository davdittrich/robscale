#ifndef ROBSCALE_ESTIMATORS_INTERNAL_H
#define ROBSCALE_ESTIMATORS_INTERNAL_H

#include "robscale_config.h"
#include "robust_core.h"
#include "vshaped_mad.h"
#include "pdq_select.h"
#include "qnsn_sort_utils.h"
#include "qnsn_kernels.h"
#include <cstring>
#include <cmath>

// Forward declarations for Qn/Sn internal implementations
// (defined in qn_estimator.cpp and sn_estimator.cpp)
namespace robscale::qnsn {
  template <typename T> double C_qn_impl(const T* x, size_t n);
  template <typename T> double C_sn_impl(const T* x, size_t n);
  template <typename T> double C_qn_impl_sorted(const T* x, size_t n);
  template <typename T> double C_qn_impl_sorted(const T* x, size_t n, const QnWorkspace* ws);
  template <typename T> double C_sn_impl_sorted(const T* x, size_t n);
  // OPT-S7: workspace-accepting overload (Tier 3 only; workspace ignored for n <= sn_stack_threshold)
  template <typename T> double C_sn_impl_sorted(const T* x, size_t n, T* workspace);
}

// rob_scale_compute: promoted from static in rob_scale.cpp
// use_avx2: pre-cached AVX2 flag (OPT-E); default false keeps old callers working.
ROBSCALE_HIDDEN
double rob_scale_compute(const double* ROBSCALE_RESTRICT data,
                         size_t n, double data_offset, double s,
                         int maxit, double tol,
                         double* ROBSCALE_RESTRICT tmp,
                         bool use_avx2 = false);

namespace robscale { namespace internal {

// GMD: sorts buf in-place, O(n log n).
// OPT-G3: optimized_sort replaces std::sort (boost::float_sort for large n).
//   Call sites (compute_all_estimators lines 216, 330) are serial — not inside
//   TBB — so tbb::parallel_sort dispatch at large n is safe here.
// OPT-G5: scale precomputed outside the accumulation loop.
inline double gmd(double* buf, int n) {
  if (n < 2) return 0.0;
  if (n <= 16)
    robscale::small_sort(buf, n);
  else
    robscale::qnsn::optimized_sort(buf, buf + n);
  const double scale = GMD_CONSISTENCY * 2.0
    / (static_cast<double>(n) * (n - 1));
  double sum = 0.0;
#if defined(_OPENMP) || defined(ROBSCALE_HAS_OMP_SIMD)
  #pragma omp simd reduction(+:sum)
#endif
  for (int i = 0; i < n; ++i)
    sum += (2.0 * (i + 1) - n - 1.0) * buf[i];
  return scale * sum;
}

// MAD from original data: fused single-buffer approach.
// After median selection, buf is a permutation of x. Deviations computed
// in-place produce the same multiset as |x[i] - med|.
inline double mad_from_data(const double* x, double* buf, int n) {
  if (n < 2) return 0.0;
  std::memcpy(buf, x, n * sizeof(double));
  double med = robscale::adaptive_median_select(buf, static_cast<size_t>(n));
  robscale::bulk_abs_diff_inplace(buf, n, med);
  return MAD_CONSISTENCY * robscale::adaptive_median_select(buf, static_cast<size_t>(n));
}

// IQR: incremental pdqselect, O(n)
// Uses buf1 only; buf2 is unused (signature kept for call-site compatibility)
// OPT-I6: n<=16 handled by small_sort fast path, keeping pdqselect code compact.
// OPT-I3: symmetric Q1 — pdqselect to lo1+1, max-scan [0..lo1] O(0.25n).
inline double iqr(const double* ROBSCALE_RESTRICT x, double* ROBSCALE_RESTRICT buf1, double* buf2, int n) {
  (void)buf2;
  if (n < 2) return 0.0;

  // OPT-I6: n<=16 sort-once-then-index fast path.
  // buf1 is used for the copy — avoids declaring a separate local buffer.
  if (n <= 16) {
    std::memcpy(buf1, x, static_cast<size_t>(n) * sizeof(double));
    double h1 = (n - 1.0) * 0.25;
    int lo1 = static_cast<int>(h1);
    double frac1 = h1 - lo1;
    double h3 = (n - 1.0) * 0.75;
    int lo3 = static_cast<int>(h3);
    double frac3 = h3 - lo3;
    robscale::small_sort(buf1, static_cast<size_t>(n));
    double q1 = buf1[lo1];
    if (frac1 > 0.0) q1 += frac1 * (buf1[lo1 + 1] - q1);
    double q3 = buf1[lo3];
    if (frac3 > 0.0 && lo3 + 1 < n) q3 += frac3 * (buf1[lo3 + 1] - q3);
    return (q3 - q1) * IQR_CONSISTENCY;
  }

  std::memcpy(buf1, x, n * sizeof(double));

  double h1 = (n - 1.0) * 0.25;
  int lo1 = static_cast<int>(h1);
  double frac1 = h1 - lo1;

  double h3 = (n - 1.0) * 0.75;
  int lo3 = static_cast<int>(h3);
  double frac3 = h3 - lo3;

  // Q1
  double q1;
  int q3_start;
  if (frac1 > 0.0) {
    // OPT-I3: select to lo1+1; max-scan [0..lo1] O(0.25n) vs O(0.75n)
    miniselect::pdqselect(buf1, buf1 + lo1 + 1, buf1 + n);
    double q1_next = buf1[lo1 + 1];
    double q1_val = buf1[0];
    for (int i = 1; i <= lo1; ++i)
      if (buf1[i] > q1_val) q1_val = buf1[i];
    q1 = q1_val + frac1 * (q1_next - q1_val);
    q3_start = lo1 + 2;
  } else {
    miniselect::pdqselect(buf1, buf1 + lo1, buf1 + n);
    q1 = buf1[lo1];
    q3_start = lo1 + 1;
  }

  // Q3: select on buf1[q3_start .. n-1]
  miniselect::pdqselect(buf1 + q3_start, buf1 + lo3, buf1 + n);
  double q3 = robscale::interp_q7(buf1, n, lo3, frac3);

  return (q3 - q1) * IQR_CONSISTENCY;
}

// IQR for pre-sorted data: O(1) direct index reads + Type 7 interpolation.
// OPT-I7: mirrors ensemble.cpp inline IQR block; no copy, no selection.
// Input MUST be sorted ascending. Undefined behaviour otherwise.
inline double iqr_sorted(const double* sorted_x, int n) {
  if (n < 2) return 0.0;
  double h1 = (n - 1.0) * 0.25;
  int lo1 = static_cast<int>(h1);
  double frac1 = h1 - lo1;
  double q1 = sorted_x[lo1];
  if (frac1 > 0.0) q1 += frac1 * (sorted_x[lo1 + 1] - q1);
  double h3 = (n - 1.0) * 0.75;
  int lo3 = static_cast<int>(h3);
  double frac3 = h3 - lo3;
  double q3 = sorted_x[lo3];
  if (frac3 > 0.0 && lo3 + 1 < n) q3 += frac3 * (sorted_x[lo3 + 1] - q3);
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

// Sn: delegates to existing optimized implementation.
// NOTE: guard retained — compute_all_estimators() is called with n-1 in the
// BCa jackknife loop (ensemble.cpp:379); n=2 input reaches sn(x, 1) here.
// sn_sorted() below is safe without a guard (bootstrap always uses full n).
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
// Pre-condition: callers must guarantee n >= 2.
// C_sn_impl_sorted retains its own n < 2 guard as defense-in-depth (returns R_NaReal).
inline double sn_sorted(const double* sorted_x, int n) {
  return robscale::qnsn::C_sn_impl_sorted<double>(sorted_x, static_cast<size_t>(n));
}

// OPT-S7: workspace-accepting overload. workspace must point to a buffer of >= n doubles.
// For n <= sn_stack_threshold, workspace is unused (stack path). For n > threshold,
// workspace replaces the inner_medians heap allocation, eliminating per-bootstrap alloc.
// Pre-condition: n >= 2 (same as sn_sorted above — no guard added here).
inline double sn_sorted(const double* sorted_x, int n, double* workspace) {
  return robscale::qnsn::C_sn_impl_sorted<double>(sorted_x, static_cast<size_t>(n), workspace);
}

inline double qn_sorted(const double* sorted_x, int n) {
  if (n < 2) return 0.0;
  return robscale::qnsn::C_qn_impl_sorted<double>(sorted_x, static_cast<size_t>(n));
}

// OPT-Q6: workspace-accepting overload — eliminates per-bootstrap heap allocation
// inside qn_refinement_kernel. ws == nullptr falls back to heap allocation.
inline double qn_sorted(const double* sorted_x, int n, const robscale::qnsn::QnWorkspace* ws) {
  if (n < 2) return 0.0;
  return robscale::qnsn::C_qn_impl_sorted<double>(sorted_x, static_cast<size_t>(n), ws);
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

  // OPT-E: cache AVX2 flag once; avoids a TLS lookup inside rob_scale_compute.
#if defined(ROBSCALE_HAS_SLEEF) && !defined(ROBSCALE_HAS_ACCELERATE) && \
    defined(ROBSCALE_HAS_AVX2_DISPATCH)
  const bool avx2 = (robscale::qnsn::RuntimeConfig::get().hw.simd_level >=
                     robscale::qnsn::SIMDLevel::AVX2);
#else
  const bool avx2 = false;
#endif
  // Newton-Raphson iteration: buf is reused as scratch (written before read)
  return rob_scale_compute(x, static_cast<size_t>(n), t, s_init, 80, 1.4901161e-8, buf, avx2);
}

// robScale on pre-sorted input (WU-RS9 / OPT-9).
// Requires: sorted_x[0..n-1] sorted ascending. buf must be >= n doubles.
// Uses O(1) median_sorted, O(log n) vshaped_mad, then standard rob_scale_compute.
// Numerically equivalent to rob_scale(sorted_x, buf, n) when input is sorted.
// Called by ensemble_one_replicate where resample is already sorted.
inline double rob_scale_sorted(const double* ROBSCALE_RESTRICT sorted_x,
                                size_t n, double* ROBSCALE_RESTRICT buf) {
  if (n < 4) return 0.0;
  const double t      = robscale::median_sorted(sorted_x, n);
  const double mad_raw = robscale::vshaped_mad(sorted_x, static_cast<int>(n), t, buf);
  const double s_init = robscale::MAD_CONSISTENCY * mad_raw;
  if (s_init <= robscale::IMPLOSION_BOUND) {
    return robscale::adm_core_sorted(sorted_x, static_cast<int>(n), t,
                                     robscale::ADM_CONSISTENCY);
  }
#if defined(ROBSCALE_HAS_SLEEF) && !defined(ROBSCALE_HAS_ACCELERATE) && \
    defined(ROBSCALE_HAS_AVX2_DISPATCH)
  const bool avx2 = (robscale::qnsn::RuntimeConfig::get().hw.simd_level >=
                     robscale::qnsn::SIMDLevel::AVX2);
#else
  const bool avx2 = false;
#endif
  return rob_scale_compute(sorted_x, n, t, s_init, 80, 1.4901161e-8, buf, avx2);
}

}} // namespace robscale::internal

#endif // ROBSCALE_ESTIMATORS_INTERNAL_H
