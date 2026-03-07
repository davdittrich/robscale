#ifndef ROBSCALE_QNSN_KERNELS_H
#define ROBSCALE_QNSN_KERNELS_H

#include <algorithm>
#include <cmath>
#include <immintrin.h>
#include <vector>

#ifdef __ARM_NEON
#include <arm_neon.h>
#endif

namespace robscale::qnsn {

// Optimized scalar version. Modern compilers (GCC 15+) auto-vectorize this better than manual intrinsics
// when provided with proper aliasing and dependence hints.
template <typename T>
void qn_brute_force_scalar(const T * __restrict__ sorted_x, size_t n, double * __restrict__ diffs) {
  size_t k = 0;
  for (size_t i = 1; i < n; ++i) {
    const double xi = static_cast<double>(sorted_x[i]);
    #pragma GCC ivdep
    for (size_t j = 0; j < i; ++j) {
      diffs[k++] = xi - static_cast<double>(sorted_x[j]);
    }
  }
}

} // namespace robscale::qnsn

#endif // ROBSCALE_QNSN_KERNELS_H
