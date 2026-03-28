#ifndef ROBSCALE_QNSN_RUNTIME_CONFIG_H
#define ROBSCALE_QNSN_RUNTIME_CONFIG_H

#include "robscale_config.h"
#include "qnsn_hardware_info.h"
#include <cmath>
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
  size_t rob_scale_parallel_threshold; // parallel fused kernel for robScale
  size_t sort_boost_threshold;
  size_t sort_tbb_threshold;
  size_t pdq_median_threshold;
  size_t pdq_robscale_threshold;  // adaptive threshold for rob_scale.cpp
  size_t pdq_lowmedian_threshold; // adaptive threshold for Sn lowmedian
  size_t pdq_qn_final_threshold;  // adaptive threshold for Qn final selection
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
    size_t per_core_l2 = hw.l2_per_core;  // Correct on all platforms

    // --- Sort thresholds ---
    sort_boost_threshold = ROBSCALE_SORT_BOOST_THRESHOLD;  // 512, not cache-sensitive

    // sort_tbb: go parallel when data exceeds per-core L2
    // Each element is 8 bytes (double), sorting needs ~2x workspace
    sort_tbb_threshold = (std::max)(size_t(4096), per_core_l2 / (sizeof(double) * 2));

    // --- Parallel thresholds ---
    // Qn: inner loop does binary searches + writes to bounds arrays (3 arrays)
    qn_parallel_threshold = (std::max)(size_t(4096),
        per_core_l2 / (sizeof(double) + 2 * sizeof(int32_t)));

    // Sn: inner loop accesses sorted_x + writes inner_medians (2 arrays)
    sn_parallel_threshold = (std::max)(size_t(4096),
        per_core_l2 / (sizeof(double) * 2));

    // rob_scale_parallel: fused AVX2 NR kernel reads w[] once per iteration.
    // With ~4-5 NR iterations typical, effective cache pressure is ~4-5 × n×8.
    // Divisor=4: threshold where one iteration's data fits in L2 with headroom
    // for prefetch and xp[] residual pressure. Floor at 4096 for TBB amortisation.
    // Note: empirical crossover may differ per machine — this is a first-principles
    // estimate. Benchmark bench/rob_scale_parallel_recalibrate.R to validate.
    rob_scale_parallel_threshold = (std::max)(size_t(4096),
        per_core_l2 / (sizeof(double) * 4));

    // --- Adaptive median-selection thresholds ---
    // Below threshold: FR-based selection wins (lower overhead at medium n).
    // Above threshold: pdqselect wins (better cache locality at large n).
    //
    // Divisor encodes the working-set pressure during selection:
    //   divisor=5: MAD/IQR (2 arrays + constants, empirically tuned)
    //   divisor=2: robScale median/MAD (1–2 warm arrays, lighter pressure)
    //   divisor=2: Sn lowmedian (1 active + 1 residual sorted_x)
    //   divisor=4: Qn final diff-window (1 active + sorted_x/work/bounds warm)
    pdq_median_threshold    = (std::max)(size_t(2048), per_core_l2 / (sizeof(double) * 5));
    pdq_robscale_threshold  = (std::max)(size_t(2048), per_core_l2 / (sizeof(double) * 2));
    pdq_lowmedian_threshold = (std::max)(size_t(2048), per_core_l2 / (sizeof(double) * 2));
    pdq_qn_final_threshold  = (std::max)(size_t(2048), per_core_l2 / (sizeof(double) * 4));

    // --- Grain size ---
    // Each grain block should fit in per-core L2
    // Workers access ~4 arrays of grain_size elements
    grain_size = (std::max)(size_t(512), per_core_l2 / (sizeof(double) * 4));
    // Cap at 8192 to ensure enough tasks for load balancing
    grain_size = (std::min)(grain_size, size_t(8192));

    // --- Fixed thresholds (not cache-sensitive) ---
    qn_exact_threshold = ROBSCALE_QN_EXACT_THRESHOLD;
    sn_stack_threshold = ROBSCALE_SN_STACK_THRESHOLD;
  }
};

} // namespace robscale::qnsn

#endif // ROBSCALE_QNSN_RUNTIME_CONFIG_H
