# Project Context (Sn Optimization Execution)

## Tooling
- Language: R package with C++ (Rcpp/RcppParallel)
- Compile gate: `R CMD INSTALL --preclean .`  (must run after each WU before tests)
- Test runner: `Rscript -e "devtools::test()"`  (full suite) or `devtools::test(filter="sn-opt")`
- Coverage: `Rscript -e "covr::package_coverage()"` — 100% threshold (all categories)
- Performance gate: `Rscript bench/sn_gate_check.R` — SOLO mode (C_sn_fast_orig removed after final PR cleanup)

## Key Files
- `src/sn_estimator.cpp` — main C++ implementation
- `src/estimators_internal.h` — internal wrappers used by ensemble
- `src/ensemble.cpp` — ensemble bootstrap (uses sn_sorted at line 145)
- `src/robscale_config.h` — ROBSCALE_RESTRICT, ROBSCALE_NOINLINE macros
- `tests/testthat/test-sn-opt.R` — correctness test suite (127 tests, 3659 total pass)
- `bench/sn_gate_check.R` — performance gate (SOLO mode after C_sn_fast_orig removal)

## Critical Constraints
- ROBSCALE_SN_STACK_THRESHOLD = 2048 (defined in robscale_config.h:34)
- sn is called at ensemble.cpp:145 (BEFORE robScale at line ~162 — work1 is free at line 145)
- SnWorker::results is declared `mutable` — retain this when adding RESTRICT

## Completed Work Units
| WU | Title | Key Changes | Commit |
|----|-------|-------------|--------|
| WU-S0 | Scaffold test + benchmark | test-sn-opt.R (113 tests), bench/sn_gate_check.R | f1b4b85 |
| WU-S1 | NOINLINE C_sn_impl | C_sn_impl_large/medium NOINLINE extracted | ab6f03d |
| WU-S2 | NOINLINE sn_kernel | sn_kernel_large NOINLINE extracted | 720e1f5 |
| WU-S3 | Micro-buffer L1 path | SN_MICRO_SIZE=128 Tier1, SN_MAX_STACK=2048 Tier2 | a617251 |
| WU-S6 | Remove redundant guards | sn_sorted() guard removed; sn() guard retained (jackknife) | 639f83c |
| WU-S4+S5 | RESTRICT + n≤16 fast path | RESTRICT on SnWorker/sn_kernel; small_sort+direct index | e3744ba |
| WU-S7 | Workspace reuse ensemble | sn_kernel/C_sn_impl_sorted workspace overloads; work1 passed | 1c693f8 |
| WU-S8 | Prefetch hints (REVERTED) | I-cache displacement regressed Tier1; reverted | (none) |
| cleanup | Remove C_sn_fast_orig | Diagnostic export removed; RcppExports regenerated | 939e0c6 |

## Established Patterns (from prior OPT sessions)
- H2H gating mandatory for compile-time-only changes (gotcha-bench-001)
- Hand-written reduction loops, never std::algorithm (gotcha-cpp-001)
- Per-estimator NOINLINE thresholds may differ from ROBSCALE_MICRO_BUFFER_SIZE
- 5% gate (ratio ≤ 1.05) is hard maximum — never widen
- Sub-µs measurements (n<200 Sn after WU-S5) are timer-quantization noise: use multi-run median
- Prefetch hints in a multi-tier function can regress OTHER tiers via I-cache displacement

## Post-Optimization Notes
- WU-S8 (prefetch hints) REVERTED: Tier2 prefetch caused I-cache displacement that consistently
  regressed Tier1 (n=17, n=100) by ~10-20%. Hardware prefetcher handles Sn's near-sequential
  access pattern adequately at these sizes.
- WU-S6 asymmetric fix: sn() guard retained for BCa jackknife path (n-1 input reachable when n=2),
  sn_sorted() guard removed (bootstrap always n>=2).
