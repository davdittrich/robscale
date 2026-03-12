#ifndef ROBSCALE_QNSN_KERNELS_H
#define ROBSCALE_QNSN_KERNELS_H

#include <algorithm>
#include <cmath>
#include <type_traits>
#if defined(__x86_64__) || defined(_M_X64) || defined(__i386__) || defined(_M_IX86)
#include <immintrin.h>
#endif

#ifdef __ARM_NEON
#include <arm_neon.h>
#endif

namespace robscale::qnsn {

// Optimized scalar version. 
template <typename T>
void qn_brute_force_scalar(const T * ROBSCALE_RESTRICT sorted_x, size_t n, double * ROBSCALE_RESTRICT diffs) {
  size_t k = 0;
  for (size_t i = 1; i < n; ++i) {
    const double xi = static_cast<double>(sorted_x[i]);
    for (size_t j = 0; j < i; ++j) {
      diffs[k++] = xi - static_cast<double>(sorted_x[j]);
    }
  }
}

#if defined(__AVX2__)
template <typename T>
void qn_brute_force_avx2(const T * ROBSCALE_RESTRICT sorted_x, size_t n, double * ROBSCALE_RESTRICT diffs) {
  size_t k = 0;
  for (size_t i = 1; i < n; ++i) {
    const double xi = static_cast<double>(sorted_x[i]);
    __m256d v_xi = _mm256_set1_pd(xi);
    
    size_t j = 0;
    // Main AVX2 loop: 8 doubles per iteration (2 registers)
    for (; j + 8 <= i; j += 8) {
      __m256d v_xj1 = {}, v_xj2 = {};
      
      if constexpr (std::is_same_v<T, double>) {
        v_xj1 = _mm256_loadu_pd(sorted_x + j);
        v_xj2 = _mm256_loadu_pd(sorted_x + j + 4);
      } else {
        v_xj1 = _mm256_set_pd(static_cast<double>(sorted_x[j+3]), static_cast<double>(sorted_x[j+2]), 
                                      static_cast<double>(sorted_x[j+1]), static_cast<double>(sorted_x[j]));
        v_xj2 = _mm256_set_pd(static_cast<double>(sorted_x[j+7]), static_cast<double>(sorted_x[j+6]), 
                                      static_cast<double>(sorted_x[j+5]), static_cast<double>(sorted_x[j+4]));
      }
      
      _mm256_storeu_pd(diffs + k, _mm256_sub_pd(v_xi, v_xj1));
      _mm256_storeu_pd(diffs + k + 4, _mm256_sub_pd(v_xi, v_xj2));
      k += 8;
    }
    // Remainder
    for (; j < i; ++j) {
      diffs[k++] = xi - static_cast<double>(sorted_x[j]);
    }
  }
}
#endif

#if defined(__ARM_NEON)
template <typename T>
void qn_brute_force_neon(const T * ROBSCALE_RESTRICT sorted_x, size_t n, double * ROBSCALE_RESTRICT diffs) {
  size_t k = 0;
  for (size_t i = 1; i < n; ++i) {
    const double xi = static_cast<double>(sorted_x[i]);
    float64x2_t v_xi = vdupq_n_f64(xi);
    
    size_t j = 0;
    for (; j + 2 <= i; j += 2) {
      float64x2_t v_xj;
      if constexpr (std::is_same_v<T, double>) {
        v_xj = vld1q_f64(sorted_x + j);
      } else {
        double tmp[2] = { static_cast<double>(sorted_x[j]),
                          static_cast<double>(sorted_x[j+1]) };
        v_xj = vld1q_f64(tmp);
      }
      vst1q_f64(diffs + k, vsubq_f64(v_xi, v_xj));
      k += 2;
    }
    for (; j < i; ++j) {
      diffs[k++] = xi - static_cast<double>(sorted_x[j]);
    }
  }
}
#endif

} // namespace robscale::qnsn

#endif // ROBSCALE_QNSN_KERNELS_H
