# Task Plan: Resolve Subtle Performance Regressions

## Goal
Identify, quantify, and resolve all performance regressions where `robscale` v0.2.0 is slower than the v0.1.5 Gold Standard. Every regression must be justified by algorithmic/architectural necessity or resolved to achieve parity/improvement.

## Phase 1: Quantify Regressions (BCa CI Audit) [x]
- [x] Extract all $(n, \text{estimator})$ pairs where `median_speedup < 1.0` vs v0.1.5.
- [x] Determine if the 95% CI overlaps with 1.0. If not, it's a confirmed regression.
- [x] Map these regions to:
    - Tiny ($n < 16$): Sorting network / Constant overhead regime.
    - Small ($16 \le n \le 128$): Exact algorithm crossover regime.
    - Medium ($128 < n \le 51,200$): Serial cache-friendly regime.
    - Large ($n > 51,200$): Parallel overhead regime.

## Phase 2: Root Cause Analysis
### Qn / Sn @ n=100
- **Variable**: `QN_EXACT_THRESHOLD` is 64. 0.1.5 might have been using a different threshold or a faster approximate kernel.
- **Action**: Profile $n=100$ and compare `qn_brute_force_exact` vs the approximate path.

### Qn / Sn @ n=1000
- **Variable**: Proximity to `SN_STACK_THRESHOLD` (2048) and `SORT_BOOST_THRESHOLD` (2048).
- **Action**: Check if Boost Spreadsort overhead is higher than `std::sort` for this specific range.

### Sn @ n > 10^5
- **Variable**: TBB grain size and false sharing in `SnWorker`. 
- **Action**: Profile the parallel worker for cache misses.

## Phase 3: Targeted Remediation [x]
- [x] Align `ROBSCALE_SORT_THRESHOLD` to 512 in `robscale_config.h`.
- [x] Align `ROBSCALE_TBB_GRAIN_SIZE` to 1024/2048 in `sn_estimator.cpp` and `qn_estimator.cpp`.
- [x] Implement dynamic grain size scaling for $n > 10^6$ in `RuntimeConfig`.
- [x] Decouple parallel sorting threshold (6144) from algorithm threshold.
- [x] Streamline `qn_brute_force_exact` by removing redundant checks.

## Phase 4: Final Gold Verification [x]
- [x] Re-run full `targets` pipeline.
- [x] Verify `detailed_gold_figure.png` shows NO regressions (all medians $\ge 1.0$).
