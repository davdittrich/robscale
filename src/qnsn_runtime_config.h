#ifndef ROBSCALE_QNSN_RUNTIME_CONFIG_H
#define ROBSCALE_QNSN_RUNTIME_CONFIG_H

#include "qnsn_hardware_info.h"
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
    size_t l2_budget = hw.l2_cache_size / 2;
    qn_exact_threshold = (size_t)std::sqrt(l2_budget / 8);

    qn_exact_threshold = (qn_exact_threshold / 8) * 8;
    if (qn_exact_threshold < 64)
      qn_exact_threshold = 64;
    if (qn_exact_threshold > 128)
      qn_exact_threshold = 128;

    sn_parallel_threshold = (hw.l2_cache_size / sizeof(double));
    if (sn_parallel_threshold < 8192)
      sn_parallel_threshold = 8192;
    if (sn_parallel_threshold > 32768)
      sn_parallel_threshold = 12288;

    qn_parallel_threshold = 8192;
    sn_stack_threshold = 2048;
  }
};

} // namespace robscale::qnsn

#endif // ROBSCALE_QNSN_RUNTIME_CONFIG_H
