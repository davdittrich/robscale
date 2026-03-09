#ifndef ROBSCALE_QNSN_RUNTIME_CONFIG_H
#define ROBSCALE_QNSN_RUNTIME_CONFIG_H

#include "qnsn_hardware_info.h"
#include "qnsn_thresholds.h"
#include <cmath>
#include <mutex>

namespace robscale::qnsn {

class RuntimeConfig {
public:
  static RuntimeConfig &get() {
    static RuntimeConfig instance;
    return instance;
  }

  // Thresholds
  size_t qn_exact_threshold;
  size_t sn_stack_threshold;
  size_t sn_parallel_threshold;
  size_t qn_parallel_threshold;

  // Hardware info
  HardwareInfo hw;

  RuntimeConfig(const RuntimeConfig &) = delete;
  RuntimeConfig &operator=(const RuntimeConfig &) = delete;

private:
  RuntimeConfig() {
    hw.discover();
    calculate_thresholds();
  }

  void calculate_thresholds() {
    qn_exact_threshold = QN_EXACT_THRESHOLD;
    sn_stack_threshold = SN_STACK_THRESHOLD;

    size_t parallel_thresh = hw.l2_cache_size / sizeof(double);
    if (parallel_thresh < 8192)
      parallel_thresh = 8192;
    if (parallel_thresh > 32768)
      parallel_thresh = 32768;

    sn_parallel_threshold = parallel_thresh;
    qn_parallel_threshold = parallel_thresh;
  }
};

} // namespace robscale::qnsn

#endif // ROBSCALE_QNSN_RUNTIME_CONFIG_H
