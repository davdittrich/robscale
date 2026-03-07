#ifndef ROBSCALE_QNSN_DISPATCHER_H
#define ROBSCALE_QNSN_DISPATCHER_H

#include "qnsn_kernels.h"
#include "qnsn_runtime_config.h"

namespace robscale::qnsn {

class Dispatcher {
public:
  template <typename T>
  static void qn_brute_force(const T *sorted_x, size_t n, double *diffs, const RuntimeConfig& config) {
    if constexpr (std::is_same_v<T, double>) {
      if (config.hw.simd_level == SIMDLevel::AVX2 && n >= 64 && n <= 1024) {
        qn_brute_force_avx2(sorted_x, n, diffs);
        return;
      }
      if (config.hw.simd_level == SIMDLevel::Neon && n >= 64 && n <= 1024) {
        qn_brute_force_neon(sorted_x, n, diffs);
        return;
      }
    }

    qn_brute_force_scalar(sorted_x, n, diffs);
  }
};

} // namespace robscale::qnsn

#endif // ROBSCALE_QNSN_DISPATCHER_H
