#ifndef ROBSCALE_QNSN_RUNTIME_CONFIG_H
#define ROBSCALE_QNSN_RUNTIME_CONFIG_H

#include "robscale_config.h"
#include "qnsn_hardware_info.h"
#include <cmath>
#include <mutex>
#include <algorithm>

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
  size_t sort_tbb_threshold;
  size_t grain_size;

  // Dynamic grain size for very large samples to avoid scheduling overhead
  size_t get_dynamic_grain_size(size_t n) const {
    if (n <= 1000000) return grain_size;
    // Aim for 8 tasks per core for good load balancing
    size_t target_tasks = hw.num_logical_cores * 8;
    if (target_tasks == 0) target_tasks = 1;
    size_t dynamic = n / target_tasks;
    size_t d = (std::max)(grain_size, dynamic);
    // Cap at 32k to avoid scheduling overhead at extreme n
    return (std::min)(d, static_cast<size_t>(32768));
  }

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
    qn_exact_threshold = ROBSCALE_QN_EXACT_THRESHOLD;
    sn_stack_threshold = ROBSCALE_SN_STACK_THRESHOLD;
    sort_tbb_threshold = ROBSCALE_SORT_TBB_THRESHOLD;
    sn_parallel_threshold = ROBSCALE_SN_PARALLEL_THRESHOLD;
    qn_parallel_threshold = ROBSCALE_QN_PARALLEL_THRESHOLD;
    
    // Default grain size
    grain_size = ROBSCALE_TBB_GRAIN_SIZE;
  }
};

} // namespace robscale::qnsn

#endif // ROBSCALE_QNSN_RUNTIME_CONFIG_H
