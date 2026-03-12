#ifndef ROBSCALE_QNSN_DISPATCHER_H
#define ROBSCALE_QNSN_DISPATCHER_H

#include "qnsn_kernels.h"
#include "qnsn_runtime_config.h"

namespace robscale::qnsn {

class Dispatcher {
public:
  template <typename T>
  static void qn_brute_force(const T *sorted_x, size_t n, double *diffs, const RuntimeConfig& config) {
    if (n <= 16) {
      qn_brute_force_scalar(sorted_x, n, diffs);
      return;
    }
#if defined(__AVX2__)
    qn_brute_force_avx2(sorted_x, n, diffs);
#elif defined(__ARM_NEON)
    qn_brute_force_neon(sorted_x, n, diffs);
#else
    qn_brute_force_scalar(sorted_x, n, diffs);
#endif
  }
};

} // namespace robscale::qnsn

#endif // ROBSCALE_QNSN_DISPATCHER_H
