#ifndef ROBSCALE_ROBUST_CORE_H
#define ROBSCALE_ROBUST_CORE_H

#include <cmath>
#include <cstring>
#include "sort_net.h"

// --- Platform-specific vectorized tanh ---
// On macOS, include only vForce (not full Accelerate) to avoid COMPLEX
// symbol collision between vecLib/vDSP.h and R's Rinternals.h
#if defined(__APPLE__) && defined(__MACH__)
  #include <vecLib/vForce.h>
  #define ROBSCALE_HAS_ACCELERATE 1
#endif

#ifdef _OPENMP
  #pragma omp declare simd notinbranch
#endif
inline double vec_tanh(double x) { return std::tanh(x); }

// Bulk tanh: vectorized on macOS (Accelerate vvtanh), SIMD-hinted elsewhere
inline void bulk_tanh(double* inout, int n) {
#if defined(ROBSCALE_HAS_ACCELERATE)
  vvtanh(inout, inout, &n);
#else
  #ifdef _OPENMP
    #pragma omp simd
  #endif
  for (int i = 0; i < n; ++i) {
    inout[i] = vec_tanh(inout[i]);
  }
#endif
}

// Key constants from Rousseeuw & Verboven (2002)
constexpr double PSI_NORM_CONST       = 0.413241928283814;
constexpr double INV_PSI_NORM_CONST   = 1.0 / 0.413241928283814;
constexpr double RHO_SCALE_CONST      = 0.37394112142347236;
constexpr double INV_RHO_SCALE_CONST  = 1.0 / 0.37394112142347236;
constexpr double MAD_CONSISTENCY      = 1.4826;
constexpr double ADM_CONSISTENCY      = 1.2533141373155001;  // sqrt(pi/2)

// Stack buffer size (512 doubles = 4 KB per segment; 3x arena = 12 KB)
constexpr int STACK_BUF_SIZE = 512;

// Median of a pre-sorted array (caller must ensure n >= 1)
inline double median_sorted(const double* x, int n) {
  if (n & 1) {
    return x[n >> 1];
  } else {
    return (x[(n >> 1) - 1] + x[n >> 1]) * 0.5;
  }
}

// ADM core: constant * mean(|x - center|)
inline double adm_core(const double* x, int n, double center, double constant) {
  double sum = 0.0;
  for (int i = 0; i < n; ++i) {
    sum += std::abs(x[i] - center);
  }
  return constant * sum / n;
}

// --- Floyd-Rivest selection (O(n) average) ---
template <typename Iter>
void floyd_rivest_select(Iter left_in, Iter k, Iter right_in) {
  using T = typename std::iterator_traits<Iter>::value_type;
  size_t n = std::distance(left_in, right_in);
  if (n < 600) {
    std::nth_element(left_in, k, right_in);
    return;
  }

  Iter left = left_in;
  Iter right = right_in - 1;

  while (right > left) {
    if (right - left > 600) {
      size_t nn = static_cast<size_t>(right - left + 1);
      size_t i = static_cast<size_t>(k - left + 1);
      double z = std::log(static_cast<double>(nn));
      double s = 0.5 * std::exp(2.0 * z / 3.0);
      double sd = 0.5 *
          std::sqrt(z * s * (static_cast<double>(nn) - s) /
                    static_cast<double>(nn)) *
          (static_cast<double>(i) - static_cast<double>(nn) / 2.0 >= 0 ? 1.0
                                                                       : -1.0);
      Iter new_left =
          (std::max)(left,
                     k - static_cast<ptrdiff_t>(static_cast<double>(i) * s /
                                                    static_cast<double>(nn) +
                                                sd));
      Iter new_right =
          (std::min)(right,
                     k + static_cast<ptrdiff_t>(static_cast<double>(nn - i) * s /
                                                    static_cast<double>(nn) +
                                                sd));
      floyd_rivest_select(new_left, k, new_right + 1);
    }

    T pivot = *k;
    Iter i = left;
    Iter j = right;
    std::swap(*left, *k);
    if (*right > pivot) std::swap(*left, *right);

    while (i < j) {
      std::swap(*i, *j);
      ++i;
      --j;
      while (*i < pivot) ++i;
      while (*j > pivot) --j;
    }

    if (*left == pivot) {
      std::swap(*left, *j);
    } else {
      ++j;
      std::swap(*j, *right);
    }

    if (j <= k) left = j + 1;
    if (k <= j) right = j - 1;
  }
}

// Median via O(n) selection
inline double median_select(double* buf, int n) {
  if (n <= 8) {
    small_sort(buf, n);
    return median_sorted(buf, n);
  }
  int k = (n - 1) / 2;
  floyd_rivest_select(buf, buf + k, buf + n);
  if (n & 1) {
    return buf[k];
  } else {
    double upper = buf[k + 1];
    for (int i = k + 2; i < n; ++i) {
      if (buf[i] < upper) upper = buf[i];
    }
    return (buf[k] + upper) * 0.5;
  }
}

// MAD via selection: |x_i - med| into dev buffer, then median_select
inline double mad_select(const double* x, int n, double med, double* dev) {
  for (int i = 0; i < n; ++i) {
    dev[i] = std::abs(x[i] - med);
  }
  return MAD_CONSISTENCY * median_select(dev, n);
}

#endif // ROBSCALE_ROBUST_CORE_H
