#ifndef ROBSCALE_QNSN_DISPATCHER_H
#define ROBSCALE_QNSN_DISPATCHER_H

#include "qnsn_kernels.h"
#include "qnsn_runtime_config.h"

namespace robscale::qnsn {

class Dispatcher {
public:
  template <typename T>
  static void qn_brute_force(const T *sorted_x, size_t n, double *diffs, const RuntimeConfig& config) {
    qn_brute_force_scalar(sorted_x, n, diffs);
  }
};

} // namespace robscale::qnsn

#endif // ROBSCALE_QNSN_DISPATCHER_H
