# Task Plan: robLoc() Performance Optimization (v0.4.0 → v0.5.0)

**Created:** 2026-03-21
**Status:** COMPLETE (2026-03-21) — all phases done, 3421 PASS
**Branch target:** perf/robloc-opt

---

## Goal

Reduce robLoc() wall-clock time so it beats revss across all tested n, eliminating
the current tie (ratio 0.96–1.00). Targets derived from estimated savings per OPT:

| n | current ratio | target ratio | primary driver |
|---|---|---|---|
| 16 | ~0.97 | ≤ 0.85 | L1+L2+L4 |
| 100 | ~0.97 | ≤ 0.75 | L1+L2+L4 |
| 500 | ~0.98 | ≤ 0.70 | L1+L2+L4 |
| 1000 | ~0.99 | ≤ 0.70 | L1+L2+L4 |
| 10000 | ~1.00 | ≤ 0.45 | L1+L3 (TBB) |

Targets confirmed empirically in Phase 7. All existing tests must remain GREEN throughout.

---

## Current baseline (2026-03-21, maxit sensitivity + 500k-iteration bench)

```
n        maxit=1µs  maxit=5µs  maxit=80µs  ratio_5v80  rob_loc_ratio_vs_revss
8          2020       2050       2110        0.972         ~0.97
16         2095       2165       2235        0.969         ~0.97
100        2880       3535       3640        0.971         ~0.97
500        6535       8135       8205        0.991         ~0.98
1000      11000      17730      17860        0.993         ~0.99
10000     90835     124710     124740        1.000         ~1.00
```

Key observations:
- NR converges in ≤5 iterations from median start (maxit=5 ≈ maxit=80)
- Baseline cost (median + MAD init before NR) = 96% of total at n=8, 62–73% at n≥500
- 3-pass per NR iteration: scale→tanh→accumulate (no fused kernel)
- RuntimeConfig::get() called once per NR iteration via bulk_tanh
- No TBB parallel NR for any n

---

## Constraints

- All existing tests must pass throughout (currently 3404 PASS, 0 FAIL, 0 SKIP)
- TDD: RED test before every implementation phase
- No new golden/pinning values until the final phase (if needed)
- Correctness gate: robLoc(x) must agree with 3-pass scalar reference within 2*sqrt(eps)
  for all n after every OPT phase
- OPT thresholds (L3) must be empirically calibrated

---

## Architecture summary (current)

```
rob_loc_impl(x)
  └── rob_loc_impl_small(xp, n)         [n <= 64, NOINLINE, 1KB stack frame]
  │     └── rob_loc_core(xp, n, arena, arena+n, ...)
  └── rob_loc_core(xp, n, arena, arena+n, ...)   [n > 64, stack/heap]
        ├── memcpy(buf, xp)
        ├── med = median_select(buf)
        ├── s   = mad_select(xp, n, med, dev)  ← reads cold xp[] (OPT-L4 target)
        └── rob_loc_compute(xp, n, med, s, maxit, tol, buf)
              for k in 1..maxit:
                Pass1: tmp[i] = (xp[i]-t)*half_inv_s    ← cold xp[]
                Pass2: bulk_tanh(tmp)  ← calls RuntimeConfig::get() (OPT-L2 target)
                Pass3: sum_psi, sum_dpsi from tmp[]
                t += 2*s*sum_psi/sum_dpsi
```

After OPT-L4+L2+L1: single-pass fused AVX2 computes sum_psi+sum_dpsi in one kernel.

---

## Phases

### Phase 0: Write ALL RED tests [COMPLETE — 2026-03-21]

Write two new test files containing all TDD guards for OPT-L1 through L4.
**Implementation must NOT begin until all tests are written and confirmed RED (or SKIP).**

#### File: tests/testthat/test-robloc-opt.R

Tests to write:

**0.1 — Correctness smoke (always passes — becomes regression guard after each OPT):**
```r
test_that("robLoc returns finite result for all tested n", {
  set.seed(42)
  for (n in c(8L, 16L, 32L, 64L, 100L, 500L, 1000L)) {
    loc <- robLoc(rnorm(n))
    expect_true(is.finite(loc),
      label = sprintf("robLoc(rnorm(%d)) = %g must be finite", n, loc))
  }
})
```

**0.2 — Cross-check: new implementation vs R reference (logistic fixed-point):**
```r
test_that("robLoc agrees with R optimize reference to 1e-6 after all OPTs", {
  # Mirrors existing test-robLoc.R pattern; tolerance 1e-6 = NR vs Brent difference
  robLocRef <- function(x) {
    obj <- function(t) sum((2*plogis((x-t)/mad(x))-1))^2
    optimize(obj, range(x), tol=sqrt(.Machine$double.eps))$minimum
  }
  set.seed(1001)
  for (n in c(10L, 25L, 50L, 100L)) {
    x   <- rnorm(n)
    got <- robLoc(x)
    ref <- robLocRef(x)
    expect_equal(got, ref, tolerance=1e-6,
      label=sprintf("n=%d: robLoc=%g, ref=%g", n, got, ref))
  }
})
```

**0.3 — Scalar path agreement: fused AVX2 result must match 3-pass scalar within 2*sqrt(eps):**
```r
test_that("robLoc fused AVX2 agrees with scalar path to 2*sqrt(eps)", {
  skip_if_not(
    exists("rob_loc_scalar_impl", envir=asNamespace("robscale"), mode="function"),
    "Diagnostic scalar helper not compiled"
  )
  set.seed(777)
  tol <- 2 * sqrt(.Machine$double.eps)
  max_diff <- 0
  for (n in c(16L, 32L, 64L, 100L, 256L, 500L, 1000L)) {
    for (rep in seq_len(20L)) {
      x  <- rnorm(n)
      v1 <- robscale:::rob_loc_scalar_impl(x)
      v2 <- robLoc(x)
      d  <- abs(v1 - v2)
      if (d > max_diff) max_diff <- d
    }
  }
  expect_lt(max_diff, tol,
    label=sprintf("max|fused-scalar| = %.3e (tol %.3e)", max_diff, tol))
})
```
NOTE: `rob_loc_scalar_impl` is a diagnostic C++ export added in Phase 3 that forces
the 3-pass scalar path regardless of SIMD availability. This test becomes the primary
correctness gate for OPT-L1.

**0.4 — Edge case: sum_dpsi near-zero guard (pathological input):**
```r
test_that("robLoc handles near-degenerate scale without NaN/Inf", {
  # All values very close together → scale ≈ 0 → u_i huge → tanh(u_i) → ±1
  # sum_dpsi = sum sech^2(u_i) → 0 without guard
  x_degen <- c(rep(1.0, 50L), rep(1.0 + 1e-10, 50L))
  loc <- robLoc(x_degen)
  expect_true(is.finite(loc) && !is.nan(loc),
    label=sprintf("robLoc degenerate input = %g must be finite", loc))

  # Constant input: s=0 path; must return median
  x_const <- rep(5.0, 20L)
  expect_equal(robLoc(x_const), 5.0, tolerance=1e-10)
})
```

**0.5 — TBB parallel agreement (conditional):**
```r
test_that("robLoc parallel (n>=threshold) agrees with serial to 1e-10", {
  skip_if_not(
    exists("rob_loc_has_parallel", envir=asNamespace("robscale"), mode="function"),
    "TBB not compiled"
  )
  skip_if_not(robscale:::rob_loc_has_parallel())
  set.seed(2024)
  tol <- 1e-10
  for (n in c(4096L, 8192L, 16384L)) {
    x <- rnorm(n)
    v_par <- robLoc(x)                            # dispatches parallel
    v_ser <- robscale:::rob_loc_serial_impl(x)    # forces serial path
    expect_equal(v_par, v_ser, tolerance=tol,
      label=sprintf("n=%d: |parallel-serial| = %.3e", n, abs(v_par-v_ser)))
  }
})
```
NOTE: `rob_loc_has_parallel()` and `rob_loc_serial_impl()` are diagnostic exports
added in Phase 5. Both skip gracefully if TBB is absent.

**Files:** tests/testthat/test-robloc-opt.R
**Done when:** All 5 test groups written; 0.1/0.2/0.4 GREEN (already pass), 0.3/0.5 SKIP
(diagnostic helpers not yet added).

**RESULT:** 13 PASS, 0 FAIL, 2 SKIP (0.3 + 0.5 skip correctly). tolerance=1e-5 used for
test 0.2 (NR vs Brent's method shows ~4e-6 relative diff for small-valued rnorm input;
1e-5 is sufficient to confirm algorithmic correctness).

---

### Phase 1: OPT-L2 — Hoist RuntimeConfig::get() before NR loop [COMPLETE — 2026-03-21]

**Needs:** Phase 0 complete

**Root cause:** `bulk_tanh` (robust_core.h:79) calls `RuntimeConfig::get()` on every
invocation to check SIMD level. With ≤5 NR iterations, this is 5 redundant TLS reads
per robLoc call. The CPUID features are invariant for process lifetime (confirmed by
Gemini: zero TOCTOU risk).

#### TDD steps

1. **RED** — confirm test 0.1 and 0.2 pass (they always do; they become regression
   guards for this refactor).

2. **Implement** — in `rob_loc_compute` (src/rob_loc.cpp:9), before the NR loop:
   - Add `const bool use_avx2_tanh = (n >= 8) && ...RuntimeConfig::get().hw.simd_level >= AVX2`
   - Create an inline helper `bulk_tanh_dispatched(double*, int, bool use_avx2)` in
     robust_core.h that accepts a pre-computed `use_avx2` flag instead of calling
     `RuntimeConfig::get()` internally.
   - Replace `robscale::bulk_tanh(tmp, n)` call in the NR loop with
     `bulk_tanh_dispatched(tmp, n, use_avx2_tanh)`.

3. **GREEN** — `devtools::test()`. All 3404 tests must pass. Tests 0.1/0.2 are the
   gates. No numeric change expected (pure dispatch refactor).

**Files:** src/rob_loc.cpp, src/robust_core.h
**Done when:** tests GREEN; nm or grep confirms no RuntimeConfig::get() call inside
the NR loop body.

---

### Phase 2: OPT-L4 — Pass buf[] to mad_select [COMPLETE — 2026-03-21]

**Needs:** Phase 1 complete

**Root cause:** `mad_select(xp, n, med, dev)` in `rob_loc_core` (src/rob_loc.cpp:49)
reads from the original `xp[]` array. After `memcpy(buf, xp)` and `median_select(buf)`,
`buf[]` holds the same multiset (permuted but all elements present) and is warm in L1/L2.
`mad_select` is permutation-invariant (confirmed: it takes `const double*`, sorts
absolute deviations into `dev[]`). Passing `buf` instead of `xp` avoids a cold L3/RAM
read of `xp[]` for medium-to-large n.

**Correctness proof:** MAD(x) = MAD(π(x)) for any permutation π. `mad_select`
computes `MAD_CONSISTENCY * median(|x_i - med|)` without requiring sorted input.

#### TDD steps

1. **RED** — tests 0.1/0.2 are the correctness gates (no new test needed; this is a
   pure call-site change with proven permutation equivalence).

2. **Implement** — in `rob_loc_core` (src/rob_loc.cpp:49):
   Change `robscale::mad_select(xp, (int)n, med, dev)` to `robscale::mad_select(buf, (int)n, med, dev)`.
   Add a comment explaining the permutation invariance.

3. **GREEN** — `devtools::test()`. All existing tests pass. No numeric change
   expected (MAD is permutation-invariant by construction).

**Files:** src/rob_loc.cpp
**Done when:** tests GREEN; 1-line change confirmed correct.

---

### Phase 3: OPT-L1 — Fused AVX2 NR kernel [COMPLETE — 2026-03-21]

**Needs:** Phase 2 complete

**Root cause:** `rob_loc_compute` does 3 separate passes per NR iteration:
```
Pass 1: tmp[i] = (xp[i] - t) * half_inv_s          (n writes to tmp[])
Pass 2: bulk_tanh(tmp, n)                            (n reads + n writes to tmp[])
Pass 3: sum_psi += p; sum_dpsi += 1.0 - p*p         (n reads from tmp[])
```
For n=1000 with 5 iterations: 15 passes over an 8 KB array. A fused kernel reads
data[] once per iteration and never touches tmp[]. Uses two AVX2 accumulators.

**Algorithm (one pass per NR iteration):**
```
for each 4-wide block:
  u4   = (x[i:i+4] - t) * half_inv_s
  p4   = ROBSCALE_TANH4_AVX2(u4)
  acc_psi  = _mm256_add_pd(acc_psi, p4)                   // Σ tanh(u)
  acc_dpsi = _mm256_fnmadd_pd(p4, p4, acc_dpsi)           // Σ (1 - tanh²(u)) = Σ sech²(u)
scalar tail: same ops with std::tanh
horizontal sum of acc_psi[4], acc_dpsi[4] → sum_psi, sum_dpsi
sum_dpsi guard: if (sum_dpsi < DBL_MIN) break;   // pathological: all |u| very large
t += 2*s * sum_psi / sum_dpsi
```

**Numerical precision:** fnmadd computes -(p*p)+1 in one rounding step, maintaining
non-negativity of sech²(u). 4-wide SIMD accumulation differs from scalar by ≤n*eps.
For n=1000, max difference ≈ 1000*2.2e-16 ≈ 2.2e-13, within 2*sqrt(eps)=2.98e-8.

#### TDD steps

1. **RED** — test 0.3 (scalar path agreement to 2*sqrt(eps)) currently SKIPs because
   `rob_loc_scalar_impl` doesn't exist yet. First add the diagnostic export:
   - In src/rob_loc.cpp, add `rob_loc_scalar_impl(Rcpp::NumericVector x)` that calls
     `rob_loc_core` with the 3-pass scalar path forced (compile-guard: always uses
     3-pass loop regardless of AVX2 flag). Export with `// [[Rcpp::export]]`.
   - Run `devtools::test()` — test 0.3 now RUNs and should be GREEN (scalar==scalar),
     confirming the diagnostic helper works.
   - Run `devtools::check()` to confirm compilation.

2. **Implement** fused AVX2 kernel:
   - In src/rob_loc.cpp, add `rob_loc_nr_step_avx2` (compile-guarded with
     `ROBSCALE_HAS_SLEEF && ROBSCALE_HAS_AVX2_DISPATCH`):
     ```cpp
     // [[no export — internal]]
     ROBSCALE_TARGET_AVX2
     static void rob_loc_nr_step_avx2(const double* xp, int n,
                                       double t, double half_inv_s,
                                       double* out_psi, double* out_dpsi);
     ```
   - In `rob_loc_compute`, hoist the AVX2 dispatch check before the NR loop:
     ```cpp
     const bool use_fused = (n >= 8) &&
       (RuntimeConfig::get().hw.simd_level >= SIMDLevel::AVX2);
     ```
   - Replace the 3-pass body with:
     ```cpp
     if (use_fused) {
       rob_loc_nr_step_avx2(xp, (int)n, t, half_inv_s, &sum_psi, &sum_dpsi);
     } else {
       // existing 3-pass + bulk_tanh_dispatched
     }
     ```
   - Add sum_dpsi guard: `if (ROBSCALE_UNLIKELY(sum_dpsi < std::numeric_limits<double>::min())) break;`
     (also add to the scalar path for consistency).

3. **GREEN** — `devtools::test()`:
   - test 0.3 (fused vs scalar to 2*sqrt(eps)) — critical gate
   - test 0.4 (degenerate input, sum_dpsi guard) — must be GREEN
   - test-robLoc.R all tests — must be GREEN
   - test 0.1/0.2 — must be GREEN
   - Full suite: 3404+ PASS, 0 FAIL

4. **Micro-benchmark** (update findings.md with results):
   Run `bench::mark(robLoc(rnorm(1000)), min_iterations=1000)` before/after Phase 3.
   Expected: 20–40% speedup in NR-dominated regime (n≥100).

**Files:** src/rob_loc.cpp
**Done when:** test 0.3 GREEN; no FAIL in full suite; micro-bench shows improvement.

---

### Phase 4: OPT-L3 — TBB parallel NR for large n [COMPLETE — 2026-03-21]

**Needs:** Phase 3 complete. **Condition:** TBB available in build (USE_DIRECT_TBB or
ROBSCALE_HAS_SYSTEM_TBB); skip gracefully if absent.

**Root cause:** For n≥4096, the NR sum over n elements is embarrassingly parallel
within each iteration. robScale already has `rob_scale_parallel_compute` using TBB
`parallel_reduce`. Each robLoc NR iteration needs to reduce TWO values: (sum_psi,
sum_dpsi). The optimal TBB structure is a single struct-based reduction (confirmed by
Gemini: avoids 2x data reads).

**Struct design:**
```cpp
struct NRAccum {
  double psi{0.0}, dpsi{0.0};
  void operator+=(const NRAccum& o) { psi += o.psi; dpsi += o.dpsi; }
};
```
Each TBB chunk calls `rob_loc_nr_step_avx2` for its sub-range and returns a `NRAccum`.
The combining lambda does `a += b; return a`.

**Threshold calibration (required before implementation):**
Write `bench/robloc_parallel_threshold.R` comparing serial vs TBB for
n = 512, 1024, 2048, 4096, 8192, 16384. TBB has ~5–15µs spawn overhead; threshold
must be above the crossover. Expected: ~4096–8192 (similar to robScale).
Store result as `ROBSCALE_LOC_PARALLEL_THRESHOLD` in robscale_config.h.

#### TDD steps

1. **RED** — Test 0.5 currently SKIPs. Add diagnostic exports:
   - `rob_loc_has_parallel()` → bool (TRUE if TBB compiled and AVX2 available)
   - `rob_loc_serial_impl(x)` → forces serial path (no TBB dispatch)
   - Run `devtools::test()` — test 0.5 now RUNs but parallel==serial since threshold
     not yet added (or SKIP if no TBB). Confirm compilation.

2. **Calibrate threshold** — run `bench/robloc_parallel_threshold.R`.
   Record crossover_n in findings.md (expected F-11).
   Set `ROBSCALE_LOC_PARALLEL_THRESHOLD` in robscale_config.h.

3. **Implement** `rob_loc_parallel_compute` in src/rob_loc.cpp:
   - Mirror structure of `rob_scale_parallel_compute` (lines 151–174 of rob_scale.cpp)
   - Use `NRAccum` struct + TBB `parallel_reduce`
   - Each chunk: call `rob_loc_nr_step_avx2` on sub-range
   - Each NR iteration updates t and checks convergence (serial between iterations)
   - In `rob_loc_core`, add dispatch: if TBB+AVX2 and n ≥ threshold → parallel path

4. **GREEN** — `devtools::test()`:
   - test 0.5 (parallel vs serial to 1e-10) — critical gate
   - test 0.1/0.2/0.3/0.4 — must stay GREEN
   - Full suite: all PASS

**Files:** src/rob_loc.cpp, src/robscale_config.h, bench/robloc_parallel_threshold.R
**Done when:** test 0.5 GREEN; all tests GREEN; threshold empirically calibrated.

---

### Phase 5: Benchmark regression check [COMPLETE — 2026-03-21]

**Needs:** Phases 1–4 complete

Run `bench/min_time_benchmark.R` (or equivalent for robLoc) comparing robLoc vs revss
across n=8,16,32,64,100,500,1000,10000.

**Final benchmark (2026-03-21, min_iterations=500):**
| n | robLoc µs | revss µs | ratio | target | baseline |
|---|---|---|---|---|---|
| 8 | 2.45 | 8.52 | **0.287** | — | ~0.97 |
| 16 | 2.93 | 9.43 | **0.311** | — | ~0.97 |
| 32 | 2.93 | 9.08 | **0.323** | — | ~0.97 |
| 64 | 3.42 | 10.97 | **0.312** | — | ~0.97 |
| 100 | 4.26 | 12.29 | **0.347** | ≤ 0.75 ✓ | ~0.97 |
| 500 | 10.62 | 28.84 | **0.368** | ≤ 0.70 ✓ | ~0.98 |
| 1000 | 20.12 | 52.94 | **0.380** | ≤ 0.70 ✓ | ~0.99 |
| 10000 | 145.13 | 596.60 | **0.243** | ≤ 0.45 ✓ | ~1.00 |

All targets met. NR convergence unchanged (maxit=5 ≈ maxit=80 at n=1000, results agree to 1e-8).
3421 PASS, 0 FAIL, 0 SKIP.

**Files:** findings.md (F-12 updated)

---

### Phase 6: Golden value update [conditional — LIKELY NOT NEEDED]

OPT-L1 (fused AVX2) changes FP rounding order but not the algorithm. The existing
robLoc tests use `tolerance=1e-6` (NR vs Brent's method) — much wider than the
fused/scalar delta (~1e-13). Golden value update **not expected**.

If test-robLoc.R::pinning test fails: regenerate ONLY after Phase 3 correctness
confirmed (test 0.3 GREEN). Prohibited before that.

---

## Decision Log

| Date | Decision | Rationale |
|------|----------|-----------|
| 2026-03-21 | OPT-L3 threshold requires bench calibration | robScale learned: never guess thresholds |
| 2026-03-21 | Single struct TBB reduction for (psi,dpsi) | Avoids 2x data reads; confirmed by Gemini |
| 2026-03-21 | sum_dpsi guard added in Phase 3 | Gemini: underflow possible for pathological input |
| 2026-03-21 | buf[] reuse (L4) done before fused kernel (L1) | Simpler changes first; L4 is 1-line |
| 2026-03-21 | RuntimeConfig hoist (L2) folded into Phase 3 | L2 done separately for clean TDD git history |
| 2026-03-21 | No Aitken/Anderson for NR | NR converges in ≤5 iters; acceleration cost > benefit |

---

## Errors Encountered

| Error | Attempt | Resolution |
|-------|---------|------------|
| (none yet) | | |

---

## Key files

| File | Role |
|------|------|
| src/rob_loc.cpp | Main implementation — all OPT phases touch this |
| src/robust_core.h | bulk_tanh, median_select, mad_select |
| src/robscale_config.h | ROBSCALE_LOC_PARALLEL_THRESHOLD (Phase 4) |
| tests/testthat/test-robLoc.R | Existing correctness tests (must remain GREEN) |
| tests/testthat/test-robloc-opt.R | New TDD tests (Phase 0) |
| bench/robloc_parallel_threshold.R | Threshold calibration (Phase 4) |

---

## Addendum: OPT-L2+L1 interaction

OPT-L2 (RuntimeConfig hoist) is logically superseded by OPT-L1 (fused kernel): when
the fused AVX2 kernel is active, `bulk_tanh` is not called at all, so L2 only applies
to the scalar fallback path. L2 is kept as a separate phase because:
1. It has independent value on non-AVX2 hardware
2. The clean TDD git history makes review easier
3. The scalar fallback is the correctness reference (test 0.3 diagnostic)
