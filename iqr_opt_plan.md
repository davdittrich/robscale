# IQR Optimization Plan: OPT-I1..I8

## Goal

Implement all identified IQR performance optimizations using TDD discipline.
**Primary:** Speed. **Secondary:** Efficiency.
**Invariant:** No performance regressions from any change in isolation.

## OPT Catalogue

| OPT   | Description                                                       | Priority  | Files            |
|-------|-------------------------------------------------------------------|-----------|------------------|
| I1+I2 | NOINLINE frame split + STACK_SIZE 4096→2048                       | HIGH      | src/iqr.cpp      |
| I4    | `interp_q7` scalar loop → `std::min_element`                      | MED-HIGH  | src/pdq_select.h |
| I5    | `ROBSCALE_RESTRICT` on NOINLINE helper + `interp_q7`              | MEDIUM    | src/iqr.cpp, src/pdq_select.h |
| I6    | n≤16 sort-once-then-index fast path                               | MEDIUM    | src/iqr.cpp, src/estimators_internal.h |
| I3    | Symmetric Q1 selection: O(0.75n)→O(0.25n) scan for frac1>0, n>16 | MED-HIGH  | src/iqr.cpp, src/estimators_internal.h |
| I7    | `iqr_sorted()` variant for pre-sorted data                        | MEDIUM    | src/estimators_internal.h |
| I8    | `ROBSCALE_RESTRICT` on `estimators_internal::iqr()`               | LOW-MED   | src/estimators_internal.h |

**Phase order rationale:** OPT-I6 precedes OPT-I3 deliberately. Once the n≤16 path routes through
`small_sort`, the pdqselect path only handles n>16. For n>16 with frac1>0, the guard
`lo3 > lo1+1` is always satisfied (for n=17: lo3=12 > lo1+1=5, and the gap only widens for larger n), eliminating a boundary edge case from I3.

---

## Benchmark Methodology

- Tool: `bench::mark(min_iterations = N, check = FALSE)` where N=2000 for n≤256, N=500 for n>256; `100L` for ensemble.
- Metric: `bm$median` in nanoseconds.
- Threshold: `1.12` for n≤128 (`.Call`-overhead-dominated); `1.05` for n>128 (signal-dominated). **NEVER CHANGE THESE.**
- Aggregation: **median across seeds per size** (robust to per-seed algorithmic variation at n=129).
- Sizes: `c(16L, 17L, 64L, 100L, 128L, 129L, 1000L, 2049L)` — covers sort/pdqselect boundary (16/17), micro boundary (128), NOINLINE boundary (128/129), large stack path (1000), stack/heap boundary (2049).
- Seeds: `c(42L, 57L, 99L, 123L, 200L, 314L, 628L, 777L, 1024L, 1618L)` — 10 seeds for stable median.
- **Gate mode**: HEAD-TO-HEAD preferred for phases with a diagnostic `*_orig` export.
  At n=16-17 (~3µs), inter-session OS scheduling noise can flip bimodal modes and produce false regressions.
  The gate script (`bench/iqr_gate_check.R`) auto-detects `iqr_impl_orig()` and uses same-session
  back-to-back `bench::mark(orig=..., new=...)` comparison when available.
  Fallback: saved-baseline comparison from `benchmarks/iqr_perf_baseline.rds`.
- All saved-baseline gates compare against `benchmarks/iqr_perf_baseline.rds` (Phase 0, 10 seeds).

---

## Phase 0: Baseline Capture (PREREQUISITE — run exactly once, before any code change)

**0.1** Create `bench/iqr_baseline.R`:
  - Benchmark `iqr_scaled(x)` for all (size × seed) combinations using above methodology.
  - Also benchmark ensemble proxy: `scale_robust(rnorm(10), n_boot=50L, seed=77)` at 100 iterations.
  - Save: `benchmarks/iqr_perf_baseline.rds` with fields `$iqr` (data.frame) and `$ens` (data.frame).

**0.2** Capture bit-identical pin:
  - `set.seed(99); pin$iqr_500 <- iqr_scaled(rnorm(500))`
  - `set.seed(42); pin$iqr_10_c1 <- iqr_scaled(rnorm(10), constant = 1.0)`
  - `set.seed(77); pin$ens_n10 <- scale_robust(rnorm(10), n_boot = 50L)` (ensemble pin, reuses MAD pin value if unchanged)
  - Save: `tests/testthat/fixtures/iqr_baseline_pin.rds`

**0.3** Create `tests/testthat/test-iqr-opt.R` with all gate tests (details below). All pin-dependent
  tests are RED until the fixture is saved; all structural tests are GREEN from the start.

**0.4** Verify: `devtools::test(filter = "iqr-opt")` — 0 failures except pin-dependent tests that
  require the fixture (they will error, not fail, until 0.2 is done).

**0.5** Git commit: `"test+bench: iqr-opt Phase 0 — baseline pin and benchmark"`

---

## Phase 1: OPT-I1+I2 — NOINLINE Frame Split + STACK_SIZE=2048

**Problem:** `src/iqr.cpp` unconditionally allocates `buf_micro[128]` (1 KB) AND `buf_stack[4096]` (32 KB)
in the same stack frame — ~33 KB dead stack on every call, even n=2. Identical to OPT-M1 for MAD
but 2× larger. `STACK_SIZE=4096` is unjustified; IQR needs one copy of data, same as MAD.

**Pattern:** Follows `mad.cpp` exactly: `static ROBSCALE_NOINLINE double iqr_impl_large(...)` +
inline micro path in `iqr_impl`.

**1.1** TDD tests (in `test-iqr-opt.R`):
  - `"NOINLINE boundary: n=128 and n=129 both correct"` — expect correctness at both sides of the
    micro/large split (n=128 ≤ ROBSCALE_MICRO_BUFFER_SIZE; n=129 > it, routes to NOINLINE helper).
    Both must equal `IQR(x, type=7) * K_IQR` within `sqrt(.Machine$double.eps)`.
  - `"Heap boundary: n=2049 correct after STACK_SIZE reduction"` — n=2049 > new STACK_SIZE=2048
    forces heap allocation; must still equal `IQR(x, type=7) * K_IQR`.
  - These tests are GREEN before and after (no behavior change, only stack layout change).

**1.2** Implementation in `src/iqr.cpp`:
  - Extract: `static ROBSCALE_NOINLINE double iqr_impl_large(const double* xp, int n, double constant)`
    with `constexpr int STACK_SIZE = 2048;` (matches MAD convention exactly).
  - Use `static constexpr int IQR_INLINE_LIMIT = 256` (NOT `ROBSCALE_MICRO_BUFFER_SIZE=128`).
    IQR_INLINE_LIMIT=256 pushes the NOINLINE boundary above the n=129 noisy ~6µs region where
    call overhead is disproportionate. ROBSCALE_MICRO_BUFFER_SIZE is shared with gmd/rob_loc/rob_scale.
  - Keep: inline micro path for n≤256 with `double buf_micro[256]` (2KB vs old 33KB dead stack).
  - Entry function `iqr_impl`: dispatch → micro path (n≤256) or `iqr_impl_large()` (n>256).
  - Remove the always-allocated `buf_stack[4096]` from the entry frame.
  - Add temporary diagnostic export `iqr_impl_orig` (original implementation) for head-to-head gate.
    Remove before commit.

**1.3** Performance gate (run `bench/iqr_baseline.R` against saved baseline):
  - For each (size, seed): compute ratio = current_median / baseline_median.
  - Gate: ratio ≤ 1.12 for size≤128; ratio ≤ 1.05 for size>128.
  - Expected direction: improvement (ratio < 1.0) for n≤128; neutral for n>128.
  - If gate fails: `git checkout src/iqr.cpp` and stop.

**1.4** All existing IQR tests GREEN: `devtools::test(filter = "iqr")`.

**1.5** Git commit: `"perf: iqr_scaled() OPT-I1+I2 — NOINLINE frame split, STACK_SIZE 4096→2048"`

---

## Phase 2: OPT-I4 — `interp_q7` Scalar Loop → `std::min_element`

**Problem:** `pdq_select.h:83-84` uses a branch-based scalar `for` loop to find the next order
statistic after a pdqselect. `std::min_element` on contiguous `double*` auto-vectorizes to MINPD
on GCC/Clang — same pattern already used in `robust_core.h:191` (`median_select`). More portable
than `#pragma omp simd reduction(min:)` for R packages (no OpenMP guarantee on CRAN macOS).

**Scope:** `interp_q7` is called ONLY from IQR paths (`src/iqr.cpp` and `src/estimators_internal.h`).
No impact on MAD or other estimators.

**2.1** TDD correctness test (in `test-iqr-opt.R`):
  - `"interp_q7 path: n values with frac>0 match stats::IQR exactly"` — covers n=2,3,4,6,7,8,10,14
    (all values where (n-1)%4 ≠ 0, forcing the interp scan). Tolerance `sqrt(.Machine$double.eps)`.
  - `"interp_q7 frac==0 path: no-scan n values still correct"` — covers n=5, n=9, n=13, n=17
    (all values where (n-1)%4 == 0, frac==0, the scan is bypassed). Must match `IQR(x, type=7) * K_IQR`.
    This guards against the `std::min_element` replacement accidentally executing when frac==0.
  - Test is GREEN before and after (no behavior change, only loop implementation change).

**2.2** Implementation in `src/pdq_select.h::interp_q7`:
  - Replace: `double nv = buf[lo + 1]; for (int i = lo + 2; i < n; ++i) if (buf[i] < nv) nv = buf[i];`
  - With: `double nv = *std::min_element(buf + lo + 1, buf + n);`
  - Add `#include <algorithm>` if not already present in the include chain.

**2.3** Performance gate vs Phase 0 baseline (same thresholds):
  - Expected direction: improvement for n≥64 (where Q1 scan ≥ 48 elements); neutral for n<64.
  - If gate fails: `git checkout src/pdq_select.h` and stop.

**2.4** All existing IQR tests GREEN: `devtools::test(filter = "iqr")`.
  Ensemble correctness: `devtools::test(filter = "ensemble")`.

**2.5** Git commit: `"perf: iqr_scaled() OPT-I4 — interp_q7 scalar loop → std::min_element"`

---

## Phase 3: OPT-I5 — `ROBSCALE_RESTRICT` Annotations

**Problem:** Neither the NOINLINE large helper (created in Phase 1) nor `interp_q7` carry
`ROBSCALE_RESTRICT`. Without RESTRICT, the compiler may conservatively assume pointer aliasing
in the scan loops and memcpy, blocking SIMD code generation. Follows OPT-M3 pattern exactly.

**3.1** TDD correctness test (in `test-iqr-opt.R`):
  - `"RESTRICT path: iqr_scaled results unchanged at n=64, 100, 1000"` — three `expect_identical`
    calls comparing to pin values captured in Phase 0. Tests that RESTRICT doesn't alter output.
  - RED until pin fixture exists; GREEN immediately after pin captured.

**3.2** Implementation:
  - In `src/iqr.cpp`: change signature to
    `static ROBSCALE_NOINLINE double iqr_impl_large(const double* ROBSCALE_RESTRICT xp, ...)`.
  - In `src/pdq_select.h::interp_q7`: change to
    `inline double interp_q7(double* ROBSCALE_RESTRICT buf, int n, int lo, double frac)`.

**3.3** Performance gate vs Phase 0 baseline (same thresholds):
  - Expected direction: marginal (~3-8% for large n); may be below noise floor.
  - If gate fails: `git checkout src/iqr.cpp src/pdq_select.h` and stop.

**3.4** All existing IQR and ensemble tests GREEN:
  - `devtools::test(filter = "iqr")` — correctness of all IQR paths.
  - `devtools::test(filter = "ensemble")` — RESTRICT in `interp_q7` (via `pdq_select.h`) is called
    from `estimators_internal.h::iqr()` which feeds `compute_all_estimators()`; ensemble must pass.

**3.5** Git commit: `"perf: iqr_scaled() OPT-I5 — ROBSCALE_RESTRICT on NOINLINE helper and interp_q7"`

---

## Phase 4: OPT-I6 — Small-n Fast Path (n≤16: sort-once-then-index)

**Problem:** For n≤16, two pdqselect calls plus two `interp_q7` scans have disproportionate
overhead relative to the amount of work. `robscale::small_sort()` (sorting network, branch-free)
is already used in `estimators_internal.h::gmd()` at line 36 and in `ensemble.cpp:87`. After a
single sort, Q1/Q3 are direct index reads — no selection, no min scan.

**Precondition:** `robscale::small_sort` is accessible in `estimators_internal.h` via the
`pdq_select.h` → `sort_net.h` include chain (verified: `estimators_internal.h:36` already calls it).

**4.1** TDD correctness tests (in `test-iqr-opt.R`):
  - `"Sort fast path: n=2..16 match stats::IQR exactly"` — all n in 2..16 must equal
    `IQR(x, type=7) * K_IQR` within `sqrt(.Machine$double.eps)` (Type 7 = direct index after sort).
  - `"Sort/pdqselect boundary: n=16 and n=17 give correct results"` — both must match base R.
    n=16 uses sort path, n=17 uses pdqselect path.
  - Tests are GREEN before and after (existing test-iqr.R already covers many of these sizes; these
    are explicit boundary guards for the new code path).

**4.2** Implementation in `src/iqr.cpp` (inside the n≤ROBSCALE_MICRO_BUFFER_SIZE micro path):
  - After `std::memcpy(buf_micro, ...)`, add:
    ```
    if (n <= 16) {
      robscale::small_sort(buf_micro, n);
      // direct index reads + Type 7 interpolation using buf_micro[lo1], buf_micro[lo1+1],
      // buf_micro[lo3], buf_micro[lo3+1]
      // return (q3 - q1) * constant;
    }
    ```
  - Indices: `lo1 = (int)((n-1)*0.25)`, `frac1 = (n-1)*0.25 - lo1`; analogous for lo3/frac3.
  - Interpolation: `q1 = buf_micro[lo1] + frac1*(buf_micro[lo1+1] - buf_micro[lo1])` if
    `frac1>0 && lo1+1<n`, else `q1 = buf_micro[lo1]`. Same for q3.
  - No `interp_q7` call needed — sorted data allows direct adjacent element read.

**4.3** Apply same fast path to `src/estimators_internal.h::iqr()`:
  - After `std::memcpy(buf1, x, ...)`, add analogous n≤16 branch.
  - `small_sort(buf1, n)` + direct index reads.

**4.4** Performance gate vs Phase 0 baseline for n=16 (threshold 1.12):
  - Also verify no regression at n=17 (first pdqselect path), n=64, n=128, n=1000, n=2049 (first heap).
  - If gate fails at n=16: `git checkout src/iqr.cpp src/estimators_internal.h` and stop.

**4.5** All existing IQR tests GREEN: `devtools::test(filter = "iqr")`.
  Ensemble pin: `devtools::test(filter = "iqr-opt")` — pin test must pass.

**4.6** Git commit: `"perf: iqr_scaled() OPT-I6 — n≤16 sort-once-then-index fast path"`

---

## Phase 5: OPT-I3 — Symmetric Q1 Selection (n>16 pdqselect path)

**Problem:** For n>16 with `frac1>0` (i.e., `(n-1)%4 ≠ 0`, roughly 75% of n values), the current
Q1 path calls `pdqselect(buf, buf+lo1, buf+n)` then `interp_q7` scans `buf[lo1+1..n-1]` (~75% of n
elements) to find the next order statistic. The symmetric approach swaps this: place the (lo1+1)-th
element, then scan only the lower 25% of n with `std::max_element`.

**Post-OPT-I6 simplification:** Since Phase 4 routes n≤16 to the sort path, the pdqselect path
now only handles n>16. For n>16 with frac1>0, the guard `lo3 > lo1+1` is always satisfied
(lo3≥12 ≥ lo1+1+1=6), so no edge-case guard is needed in the symmetric path — just check `frac1>0`.

**Scan comparison for frac1>0 (n>16):**

| Approach   | Q1 select   | Q1 scan           | Q3 select    | Q3 scan           |
|------------|-------------|-------------------|--------------|-------------------|
| Current    | O(n)        | O(0.75n) scalar   | O(0.75n)     | O(0.25n) vector   |
| Symmetric  | O(n)        | O(0.25n) vector   | O(0.75n)     | O(0.25n) vector   |

**5.1** TDD correctness tests (in `test-iqr-opt.R`):
  - `"Symmetric Q1: all four (n-1)%4 residue classes correct"`:
    - Residue 0: n=17 (frac=0, standard path) — must equal IQR(x,type=7)*K_IQR.
    - Residue 1: n=18 (frac=0.25, symmetric) — must equal IQR(x,type=7)*K_IQR.
    - Residue 2: n=19 (frac=0.5, symmetric) — must equal IQR(x,type=7)*K_IQR.
    - Residue 3: n=20 (frac=0.75, symmetric) — must equal IQR(x,type=7)*K_IQR.
    - Also: n=64 (residue 3, frac>0, typical case), n=100 (residue 3), n=1000 (residue 3).
  - `"Symmetric Q1 known result: n=20"`:
    - Fixed vector `x <- c(1,2,3,...,20)` (sorted). IQR type-7 = (15.25-5.75)*K_IQR.
    - `expect_equal(iqr_scaled(x), IQR(x,type=7)*K_IQR, tolerance=sqrt(.Machine$double.eps))`.
  - These tests are GREEN before implementation (correctness-only guards).

**5.2** Implementation in `src/iqr.cpp` (pdqselect path, n>16 after OPT-I6 branch):
  ```
  if (frac1 > 0.0) {
    // Symmetric: place (lo1+1)-th element; max-scan lower partition for Q1
    miniselect::pdqselect(buf, buf + lo1 + 1, buf + n);
    double q1_next = buf[lo1 + 1];
    double q1_val  = *std::max_element(buf, buf + lo1 + 1);  // O(0.25n) vectorized
    q1 = q1_val + frac1 * (q1_next - q1_val);
    // Q3: start from lo1+2 (lo1+1 already placed)
    miniselect::pdqselect(buf + lo1 + 2, buf + lo3, buf + n);
  } else {
    miniselect::pdqselect(buf, buf + lo1, buf + n);
    q1 = buf[lo1];
    miniselect::pdqselect(buf + lo1 + 1, buf + lo3, buf + n);
  }
  // Q3 interpolation (both paths):
  q3 = buf[lo3];
  if (frac3 > 0.0 && lo3 + 1 < n)
    q3 += frac3 * (*std::min_element(buf + lo3 + 1, buf + n) - q3);
  return (q3 - q1) * constant;
  ```
  Note: when frac1==0 (residue 0), frac3==0 too — no Q3 scan needed.

**5.3** Apply same symmetric logic to `src/estimators_internal.h::iqr()` for n>16:
  - After the n≤16 sort branch from Phase 4, the remaining pdqselect path uses the same symmetric logic.

**5.4** Performance gate vs Phase 0 baseline:
  - Primary measurement: n=64, n=128, n=1000 (all residue 3, frac>0, symmetric path active).
  - Expected: improvement (ratio < 1.0) at n≥64.
  - Also verify n=16, n=17 (sort path + standard pdqselect — should be unchanged).
  - If gate fails: `git checkout src/iqr.cpp src/estimators_internal.h` and stop.

**5.5** All existing IQR tests GREEN: `devtools::test(filter = "iqr")`.
  Pin test: `devtools::test(filter = "iqr-opt")`.

**5.6** Git commit: `"perf: iqr_scaled() OPT-I3 — symmetric Q1 selection, O(0.75n)→O(0.25n) scan"`

---

## Phase 6: OPT-I7 — `iqr_sorted()` Variant

**Problem:** `estimators_internal.h` exposes `sn_sorted()` and `qn_sorted()` but no `iqr_sorted()`.
The ensemble already implements O(1) IQR directly inline (`ensemble.cpp:123-140`); this phase
creates a named function mirroring that block for API completeness and future callers.

**6.1** TDD correctness test (RED until implementation):
  - `"iqr_sorted exists and matches iqr_scaled on sorted input"`:
    ```R
    # Verify via diagnostic export (temporary iqr_sorted_impl, removed in Phase 8 or at end of phase)
    set.seed(42); x <- sort(rnorm(100))
    expect_equal(iqr_sorted_impl(x), iqr_scaled(x), tolerance = sqrt(.Machine$double.eps))
    ```
  - Also verify on n=2,3,4,5,16,17,100,1000 sorted inputs.

**6.2** Diagnostic export (temporary — same pattern as `mad_from_data_bench` in MAD opt plan):
  - Add `// [[Rcpp::export]]` to a thin wrapper `iqr_sorted_impl` in `src/iqr.cpp` that calls
    `robscale::internal::iqr_sorted(x.begin(), n)`. Remove in Phase 8.

**6.3** Implementation of `iqr_sorted()` in `src/estimators_internal.h`:
  - New function mirroring `ensemble.cpp:123-140`:
    ```cpp
    inline double iqr_sorted(const double* sorted_x, int n) {
      if (n < 2) return 0.0;
      double h1 = (n - 1.0) * 0.25; int lo1 = (int)h1; double frac1 = h1 - lo1;
      double q1 = sorted_x[lo1];
      if (frac1 > 0.0 && lo1 + 1 < n) q1 += frac1 * (sorted_x[lo1+1] - q1);
      double h3 = (n - 1.0) * 0.75; int lo3 = (int)h3; double frac3 = h3 - lo3;
      double q3 = sorted_x[lo3];
      if (frac3 > 0.0 && lo3 + 1 < n) q3 += frac3 * (sorted_x[lo3+1] - q3);
      return (q3 - q1) * IQR_CONSISTENCY;
    }
    ```

**6.4** No performance regression gate needed (new function; no change to existing paths).
  Verify: all existing IQR and ensemble tests GREEN:
  - `devtools::test(filter = "iqr")` — no existing IQR path changed.
  - `devtools::test(filter = "ensemble")` — no ensemble path changed.
  - If any test fails: `git checkout src/iqr.cpp src/estimators_internal.h` and diagnose before retrying.

**6.5** Remove diagnostic `iqr_sorted_impl` export from `src/iqr.cpp` before commit.

**6.6** Git commit: `"feat: iqr_sorted() — O(1) IQR variant for pre-sorted data"`

---

## Phase 7: OPT-I8 — RESTRICT on `estimators_internal::iqr()`

**Problem:** `estimators_internal.h::iqr(const double* x, double* buf1, double* buf2, int n)` lacks
`ROBSCALE_RESTRICT` on its pointer parameters. In the `compute_all_estimators()` call, `x` and
`buf1` are always distinct allocations, but the compiler cannot prove this without RESTRICT.

**7.1** TDD test (in `test-iqr-opt.R`):
  - `"ensemble iqr RESTRICT: pin unchanged"` — `expect_identical` on `scale_robust(rnorm(10), n_boot=50)`
    against pin$ens_n10. Verifies RESTRICT doesn't alter output.

**7.2** Implementation in `src/estimators_internal.h::iqr()`:
  - Change signature:
    `inline double iqr(const double* ROBSCALE_RESTRICT x, double* ROBSCALE_RESTRICT buf1, double* buf2, int n)`

**7.3** Performance gate (ensemble proxy):
  - `scale_robust(rnorm(10), n_boot=50L)` — compare median to Phase 0 ensemble baseline.
  - Threshold: 1.12 (small n ensemble is noisy).
  - If gate fails: `git checkout src/estimators_internal.h` and stop.

**7.4** All tests GREEN: `devtools::test(filter = "iqr")` + `devtools::test(filter = "ensemble")`.

**7.5** Git commit: `"perf: iqr_scaled() OPT-I8 — ROBSCALE_RESTRICT on estimators_internal::iqr()"`

---

## Phase 8: Final Regression Check + NEWS Update

**8.1** Run full correctness suite: `devtools::test()` — all tests GREEN, 0 failures.

**8.2** Pin verification: `devtools::test(filter = "iqr-opt")` — all 7 test blocks GREEN.

**8.3** Full performance sweep vs Phase 0 baseline across all sizes/seeds — document speedup table.

**8.4** Update `NEWS.md`: add entry for `iqr_scaled()` OPT-I1..I8 with expected speedup ranges.

**8.5** Git commit: `"docs: IQR OPT-I1..I8 complete — NEWS.md update"`

---

## `tests/testthat/test-iqr-opt.R` Structure

Tests listed in order (all must be GREEN after Phase 8):

```
Test I.1:  Pin unchanged after OPT-I1..I8                        [RED until Phase 0.2]
Test I.2:  NOINLINE boundary: n=128 and n=129 correct             [GREEN throughout]
Test I.3:  Heap boundary: n=2049 correct after STACK_SIZE=2048    [GREEN throughout]
Test I.4a: interp_q7 frac>0 path: n with (n-1)%4≠0 correct       [GREEN throughout]
Test I.4b: interp_q7 frac==0 path: n=5,9,13,17 no-scan correct   [GREEN throughout]
Test I.5:  RESTRICT path: iqr_scaled results identical at pins    [RED until Phase 0.2]
Test I.6:  Sort fast path: n=2..16 match stats::IQR exactly       [GREEN throughout]
Test I.7:  Sort/pdqselect boundary: n=16 and n=17 correct         [GREEN throughout]
Test I.8:  Symmetric Q1: all four (n-1)%4 residue classes         [GREEN throughout]
Test I.9:  Symmetric Q1 known result: n=20                        [GREEN throughout]
Test I.10: iqr_sorted exists and is correct (via diagnostic export)[RED until Phase 6.2]
Test I.11: Ensemble pin unchanged after all IQR changes            [RED until Phase 0.2]
```

---

## Rollback Protocol

If any performance gate fails:
1. `git checkout <affected files>` — revert to pre-phase state.
2. Document the failure in `progress.md` with actual ratio and threshold.
3. Investigate root cause before attempting the phase again.
4. Do NOT proceed to the next phase with a failing gate.

---

## Files Created/Modified Summary

| Action | File |
|--------|------|
| CREATE | `bench/iqr_baseline.R` |
| CREATE | `tests/testthat/test-iqr-opt.R` |
| CREATE | `tests/testthat/fixtures/iqr_baseline_pin.rds` (generated) |
| CREATE | `benchmarks/iqr_perf_baseline.rds` (generated) |
| MODIFY | `src/iqr.cpp` (Phases 1, 3, 4, 5, 6) |
| MODIFY | `src/pdq_select.h` (Phases 2, 3) |
| MODIFY | `src/estimators_internal.h` (Phases 4, 5, 6, 7) |
| MODIFY | `NEWS.md` (Phase 8) |
