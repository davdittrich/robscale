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
#ifdef ROBSCALE_HAS_AVX2_DISPATCH
    // Runtime dispatch: use AVX2 kernel if hardware supports it.
    // config.hw.simd_level is populated once at package load via
    // __builtin_cpu_supports("avx2") — no per-call CPUID overhead.
    if (config.hw.simd_level >= SIMDLevel::AVX2) {
      qn_brute_force_avx2(sorted_x, n, diffs);
      return;
    }
#endif
#if defined(__ARM_NEON)
    qn_brute_force_neon(sorted_x, n, diffs);
#else
    qn_brute_force_scalar(sorted_x, n, diffs);
#endif
  }
};

} // namespace robscale::qnsn

#endif // ROBSCALE_QNSN_DISPATCHER_H
