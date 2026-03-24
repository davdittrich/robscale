# Findings: robScale() Performance Analysis (2026-03-21)

## Previous findings preserved in git history. This section covers the 2026-03-21
## performance analysis for the n<=64 erratic scaling and n=500-1000 regression.

---

## F-1: Benchmark Evidence

### Min-time benchmark (bench/min_time_benchmark.R, 2026-03-21)
revss = Fortran package, robscale = this package (C++ with AVX2/SLEEF/TBB)

| n | robscale µs | revss µs | ratio (rob/rev) |
|---|---|---|---|
| 8 | 5.135 | 12.210 | **0.421** (2.4× faster) |
| 16 | 6.585 | 12.990 | **0.507** (2.0× faster) |
| 32 | 11.710 | 16.735 | **0.700** (1.4× faster) |
| 64 | 16.540 | 20.340 | **0.813** (1.2× faster) |
| 100 | 20.400 | 23.360 | **0.873** |
| 500 | 90.640 | 80.520 | **1.126** ← REGRESSION |
| 1000 | 157.280 | 140.220 | **1.122** ← REGRESSION |
| 10000 | 1233.600 | 1199.800 | **1.028** (tied) |
| 100000 | 2362.000 | 14556.000 | **0.162** (TBB kicks in, 6×) |
| 1000000 | 17480.000 | 145480.000 | **0.120** (TBB, 8×) |

**Two distinct regimes:** n≤64 (speedup degrades monotonically) and n=500-1000 (regression).

---

## F-2: Root Cause — n≤64 Erratic Scaling

### Mechanism
`rob_scale_core` dispatches to `median_net(x, n)` for n ≤ `ROBSCALE_SORT_MEDIAN_THRESHOLD = 64`.
`median_net` is called **twice per call** (once for median of x, once for median of |x-t| for MAD).

The partial-sort network comparator counts (from src/sort_net.h, measured 2026-03-21):

| n | median_net swaps | ×2 (MAD total) | swaps/element |
|---|---|---|---|
| 8 | 16 | 32 | 2.0 |
| 16 | 46 | 92 | 2.9 |
| 32 | 128 | 256 | 4.0 |
| 64 | 337 | 674 | 5.3 |

Growth is **O(n^1.5)** — superlinear. Floyd-Rivest comparisons scale O(n):
- FR(n=8): ~12 comparisons ×2 = 24 total (vs 32 for net) — FR already faster
- FR(n=32): ~48 comparisons ×2 = 96 total (vs 256 for net) — 2.7× advantage
- FR(n=64): ~96 comparisons ×2 = 192 total (vs 674 for net) — 3.5× advantage

**Conclusion:** `ROBSCALE_SORT_MEDIAN_THRESHOLD = 64` is far too high.
The optimal crossover is estimated at n=16 (FR overhead breaks even with network at ~46 swaps).
This must be confirmed empirically by Phase 2 benchmark.

### Secondary: PLT overhead for median_net<double>
`extern template double median_net<double>` in sort_net.h causes PLT indirect call on Linux.
The comment in sort_net.h:12522 explicitly acknowledges this.
At n=8 where median_net_8 is only 16 swaps, a PLT overhead of ~5-10 cycles is ~2-3% of the
function body. Affects 2 calls per rob_scale_core. Fix: visibility("hidden") in sort_net_inst.cpp.

---

## F-3: Root Cause — n=500–1000 Regression

### What's known
- At n=500, rob_scale_core uses adaptive_robscale_median_select (Floyd-Rivest, O(n))
- M-scale iteration: 30-50 steps documented in aitken_findings.md
- revss is Fortran (confirmed via nm: robScale_f_, __sorting_MOD_oqs)
- Aitken Δ² acceleration is implemented in aitken_iterate; fires when sequence is monotone+contracting

### Hypotheses (ranked by likelihood)
1. **High iteration count at medium n**: if robscale does 40 iterations and revss does 25
   (different convergence criterion or initial estimate), the 12% gap is explained.
2. **Aitken guard too strict**: `d1*d2 > 0` requires monotone sequence. At n=500-1000,
   early iterations may oscillate slightly (non-monotone) → Aitken doesn't fire →
   no iteration reduction. Revss may implicitly benefit from a smoother iteration path.
3. **revss uses a better initial scale estimate**: revss v3.0.0 changed to bias-corrected
   "AA" constants (madn) as starting value. This may start closer to the fixed point.

### What needs measurement (Phase 0+1)
- Actual mean iteration count at n=100, 200, 500, 1000
- Aitken fire rate (fraction of outer iterations where Aitken extrapolation was accepted)
- Total rho_sum evaluations vs expected

**Status:** not yet measured. Phase 0 must add diagnostic instrumentation first.

---

## F-4: Iteration Count Context (from aitken_findings.md)
- Mean 32–54 iterations for M-scale, frequent maxit=80 cap hits at small n
- Aitken reduces iteration count by 30–50% for small n, ~20% for large n
- At medium n (500-1000), the benefit is expected ~20%, but may be less if guard rejects often
- The aitken_findings.md documents maximum observed deviation: ~2×sqrt(eps) (within tolerance)

---

## F-5: revss Architecture
revss v3.1.0 is **Fortran-based** (from nm -D revss.so):
- `robScale_f_` (Fortran entry point)
- `robScale_c` (C wrapper)
- `__sorting_MOD_oqs` (Fortran order-statistics module)
- `__roblocscale_MOD_*` (Fortran module internals)

revss v3.0.0 changed algorithm defaults (confirmed in test-cross-check.R comment):
> revss v3.0.0 changed robScale() and robLoc() defaults:
>   - robScale() uses bias-corrected "AA" constants (madn) as starting value
This means revss's initial scale estimate changed in v3.0.0. This may be why revss
converges faster at n=500-1000 — it starts closer to the fixed point.

---

## F-6: Test Landscape (from inspection of tests/testthat/)

### Tests that are algorithm-sensitive (must stay green throughout)
- `test-robScale.R::pinning test` — golden values with tolerance `2*sqrt(eps)`, Aitken-aware
- `test-robScale.R::robScale converges to true fixed point` — correctness certificate
- `test-cross-check.R::robLoc golden values` — tolerance `1e-12` (very strict for robLoc only)

### Tests that are algorithm-independent (structural)
- `test-robscale-pdqselect.R` — non-NA positive result, determinism, boundary tests
- `test-build-config.R` — config thresholds, build metadata

### Key tolerance information
- Pinning test: `2*sqrt(eps)` ≈ 3e-8 — can accommodate Aitken trajectory changes
- Cross-check: `1e-12` for robLoc, but revss comparison is SKIPPED for revss >= 3.0.0
- Fixed-point test: `2*sqrt(eps)` — algorithm-independent correctness certificate

### IMPORTANT: What threshold changes do NOT affect
Changing `ROBSCALE_SORT_MEDIAN_THRESHOLD` from 64 to ~16 changes WHICH selection
algorithm is used for the MAD initialization (n=17..64: FR vs net). Both algorithms
find the exact same median value (within floating point). The M-scale fixed-point
iteration itself is unchanged. Therefore:
- The fixed-point test will still pass (same algorithm, same convergence)
- The pinning test should still pass (same trajectory, different initialization path
  that produces the same s_init with negligible floating-point difference)
- **Risk**: floating-point differences in FR vs median_net for the same input could
  cause median to differ by ±1 ULP (last bit). This would change s_init by ~1e-16,
  which after 30-50 iterations propagates to < 1e-14 in final s. Within 2*sqrt(eps).

---

## F-7: Comparator count methodology
Counts measured by parsing SWAP_IF_GREATER macro calls in sort_net.h using Python:
```python
# See task_plan.md Phase 2 for the full script
# n=8: 16, n=16: 46, n=32: 128, n=64: 337 (partial sort networks, not full sort)
```
These are PARTIAL sort networks (optimized for median extraction only), not full
sort_net_N networks. They are already significantly pruned from the full sort.

---

## F-8: Gemini Collaboration (2026-03-21)
**SESSION_ID:** 3ce30e4d-cb08-4696-b94d-30e7ab3cc0a0 (pending full response)
**Key questions asked:**
1. Is OPT-A (lower threshold ~16) correctly identified as primary fix?
2. Additional opportunities missed?
3. Other causes for n=500-1000 regression?
4. TDD strategy for threshold changes in numerical packages?

**Status:** Awaiting full Gemini response. Will update this section when complete.

---

## F-9: Iteration Diagnostics (Phase 0+1 — MEASURED 2026-03-21)
Run: bench/iteration_diagnostics.R (n_reps=200, set.seed(42))

| n | mean_oi | sd_oi | max_oi | aitken_rate | rho_evals | converge% |
|---|---|---|---|---|---|---|
| 10 | 11.96 | 11.63 | 40 | 47.5% | 23.1 | 97.0% |
| 20 | 7.42 | 7.68 | 35 | 57.6% | 14.0 | 100% |
| 32 | 6.46 | 6.67 | 30 | 61.5% | 12.0 | 100% |
| 64 | 4.25 | 2.88 | 22 | 68.6% | 7.6 | 100% |
| 100 | 3.97 | 1.97 | 24 | 70.4% | 7.0 | 100% |
| 200 | 3.50 | 0.53 | 5 | 70.6% | 6.0 | 100% |
| 500 | 3.23 | 0.46 | 4 | 68.4% | 5.5 | 100% |
| 1000 | 3.08 | 0.38 | 4 | 67.0% | 5.2 | 100% |
| 2000 | 3.02 | 0.21 | 4 | 66.7% | 5.0 | 100% |
| 5000 | 2.95 | 0.22 | 3 | 65.8% | 4.9 | 100% |

**Key finding:** Only ~5 rho_sum calls at n=500-1000. Aitken fires 67-68%.
**The iteration count is NOT the bottleneck for the n=500-1000 regression.**
Previously hypothesised 30-50 iterations (from aitken_findings.md) is WRONG for current code.

---

## F-10: Threshold Crossover (Phase 2 — MEASURED 2026-03-21)
Run: bench/threshold_crossover_bench.R (n_iters=500000, n_runs=3, min-of-runs)

| n | net_ns | fr_ns | ratio | winner |
|---|---|---|---|---|
| 8 | 910.0 | 912.0 | 0.998 | tie |
| 10 | 928.0 | 932.0 | 0.996 | tie |
| 12 | 934.0 | 946.0 | 0.987 | tie |
| 14 | 920.0 | 930.0 | 0.989 | tie |
| 16 | 946.0 | 964.0 | 0.981 | tie |
| 18 | 960.0 | 950.0 | 1.011 | tie |
| 20 | 946.0 | 946.0 | 1.000 | tie |
| 24 | 968.0 | 952.0 | 1.017 | tie |
| 28 | 1004.0 | 984.0 | 1.020 | tie |
| 32 | 1026.0 | 1034.0 | 0.992 | tie |
| 40 | 1042.0 | 1014.0 | 1.028 | tie |
| 48 | 1066.0 | 1030.0 | 1.035 | tie |
| 56 | 1100.0 | 1104.0 | 0.996 | tie |
| 64 | 1160.0 | 1062.0 | 1.092 | **FR** |

**Caveat:** R→Rcpp call overhead (~900 ns baseline) dominates algorithm cost.
The "tie" results at n=8..56 are inflated by this overhead, not meaningful.
True C++ crossover (no R overhead) is estimated at n~16-22 from the marginal
cost analysis: net increments 250 ns from n=8→64, FR increments 150 ns — diverging
from n~24-32. Theoretical analysis (F-2 comparator counts) puts crossover at n~16.

**Decision:** Set `ROBSCALE_SORT_MEDIAN_THRESHOLD = 16` (theoretical optimum,
Gemini-confirmed, conservative re benchmark noise). Numerical agreement confirmed:
max absolute difference across all n=8..64: 0.000e+00 (bit-exact agreement).

---
| 16 | TBD | TBD | TBD | TBD |
| 20 | TBD | TBD | TBD | TBD |
| 24 | TBD | TBD | TBD | TBD |
| 32 | TBD | TBD | TBD | TBD |
| 48 | TBD | TBD | TBD | TBD |
| 64 | TBD | TBD | TBD | TBD |

**crossover_n = TBD** (to be set after Phase 2 benchmark)

---

## F-8 (Updated): Gemini Analysis (SESSION_ID: 3ce30e4d-cb08-4696-b94d-30e7ab3cc0a0)

### Confirmations
- **OPT-A threshold=16 confirmed** as the primary fix. Gemini: "A threshold of 16 is ideal.
  It aligns with the `pdq_median_select` logic in `src/pdq_select.h` and is the point where
  the branch-predictability of sorting networks typically stops outweighing the O(n^1.5)
  comparator growth."
- **PLT fix (OPT-B) confirmed** as a worthwhile secondary optimization.

### New opportunity: OPT-F — Pass dev[] to rob_scale_compute (eliminate per-iter subtraction)
`rob_scale_core` already computes `dev[i] = |x[i] - t|` (absolute deviations from median).
`rob_scale_compute` then re-computes `(data[i] - data_offset)` per element per iteration.
**Key insight:** since rho_sum = sum(tanh(...)^2) and tanh is odd, we have:
  tanh((x[i] - t) / (2Cs))^2 = tanh(|x[i] - t| / (2Cs))^2 = tanh(dev[i] / (2Cs))^2
So passing `dev[]` as the data buffer with `data_offset = 0.0` gives IDENTICAL rho_sum.
Saves 1 subtraction per element per iteration → for n=1000, 40 iters: 40,000 subtractions.
Estimated gain: ~5–15% at n=500–1000 (subtraction is cheap but compounds over iterations).
**Risk:** This changes `rob_scale_compute`'s API slightly (caller must pre-compute dev[]).
Must verify that the AVX2 fused kernel handles data_offset=0 correctly (it does — it's
just `_mm256_sub_pd(d, off4)` with off4 = 0).

### New opportunity: OPT-H — Fuse median+MAD for small n
For n ≤ threshold (post-OPT-A: n ≤ 16), one full sort_net pass gives sorted x.
From sorted x: median is x[n/2] in O(1). MAD = median(|x[i] - median|) can be
computed via a merge-style O(n) algorithm on the sorted array (deviations form a
"V-shape" around the median position, merging from both ends gives sorted deviations).
Cost: sort_net_16 (60 comparators) + O(n) MAD = 60 + 16 = ~76 operations total
vs current: median_net_16 (46) + 16 deviations + median_net_16 (46) = 108 operations
Saves ~30% for n=16. At n ≤ 8 the gain is smaller. **Modest benefit, medium effort.**

### Gemini's n=500–1000 regression hypotheses
1. **tanh accuracy**: SLEEF provides 1.0 ULP accuracy (computationally expensive).
   revss uses Fortran's DTANH which may use a faster approximation.
   → Verify: time tanh alone vs glibc/sleef at n=500.
2. **Function call overhead**: the lambda `rho_sum` in rob_scale_compute may not inline
   perfectly in all compilers. Check assembly to confirm inlining.
3. Note: Gemini initially suggested "parallel overhead" but this is wrong — the parallel
   threshold is max(4096, L2/32) >> 1000, so TBB is NOT triggered at n=500-1000.

### Gemini's TDD recommendations
1. **Differential verification**: run with both thresholds (16 and 64), assert outputs
   agree to within 1e-15 for n ∈ [17, 64]. This is STRICTER than 2*sqrt(eps) and confirms
   FR gives the same median as median_net (same floating-point order statistic).
2. **Property-based invariance**: verify robScale(x+k) = robScale(x) and
   robScale(a*x) = |a|*robScale(x) across threshold boundary n=15,16,17.
3. **Iteration count tracking**: add debug mode to aitken_iterate to return iteration count.

### Additional OPT ranking (updated with Gemini input)
| Priority | OPT | Expected gain | Notes |
|---|---|---|---|
| 1 | OPT-A: lower threshold 64→16 | 2–3× MAD init at n=32..64 | Phase 3 |
| 2 | OPT-F: pass dev[] to rho_sum | 5–15% at n=500–1000 | Phase 5 addition |
| 3 | OPT-C: Aitken guard or init | 5–20% at n=500–1000 | Phase 5 |
| 4 | OPT-B: visibility("hidden") | ~10 ns per call | Phase 4 |
| 5 | OPT-H: fuse median+MAD (n≤16) | ~30% for n=16 | Phase 6 (conditional) |
| 6 | OPT-E: cache RuntimeConfig | ~5 ns per call | Phase 7 |

---

# Findings: robLoc() Performance Analysis (2026-03-21)

---

## F-11: robLoc() Baseline Measurement

### maxit sensitivity (500k iterations each, 2026-03-21)

Measuring: overhead of initialization (median + MAD) vs NR iteration cost.

| n | maxit=1 µs | maxit=5 µs | maxit=80 µs | ratio_5v80 |
|---|---|---|---|---|
| 8 | 2020 | 2050 | 2110 | 0.972 |
| 16 | 2095 | 2165 | 2235 | 0.969 |
| 32 | 2290 | 2470 | 2560 | 0.965 |
| 100 | 2880 | 3535 | 3640 | 0.971 |
| 500 | 6535 | 8135 | 8205 | 0.991 |
| 1000 | 11000 | 17730 | 17860 | 0.993 |
| 10000 | 90835 | 124710 | 124740 | 1.000 |

**Key conclusions:**
- maxit=5 ≈ maxit=80 at all n → NR converges in ≤5 iterations from median start
- Baseline (maxit=1) = 96% of total at n=8, decreasing to 62–73% at n≥500
- NR iterations are NOT the bottleneck at small n; they dominate at large n
- Per-iteration cost: n=8→7.5ns, n=1000→1683ns, n=10000→8469ns
- robLoc ratio vs revss: 0.96–1.00 (tied) at all tested n

### robLoc vs revss comparison

| n | robLoc µs | revss µs | ratio |
|---|---|---|---|
| 8 | ~2.1 | ~2.2 | ~0.97 |
| 100 | ~3.6 | ~3.7 | ~0.97 |
| 500 | ~8.2 | ~8.4 | ~0.98 |
| 1000 | ~17.9 | ~18.1 | ~0.99 |
| 10000 | ~124.7 | ~124.5 | ~1.00 |

---

## F-12: robLoc() Root Cause Analysis

### Current 3-pass NR loop (src/rob_loc.cpp:16-30)

```
for k in 1..maxit:
  Pass 1: tmp[i] = (xp[i] - t) * half_inv_s         [n writes to tmp[]]
  Pass 2: bulk_tanh(tmp, n)                           [n reads + n writes to tmp[]; calls RuntimeConfig::get()]
  Pass 3: sum_psi += p; sum_dpsi += 1-p*p             [n reads from tmp[]]
  t += 2*s * sum_psi / sum_dpsi
```

For n=1000, 5 iterations: 15 passes over an 8KB array (fits L1, but is still 3x more bandwidth than needed).

### Identified opportunities (ranked by expected impact)

| OPT | Description | Expected gain | Effort |
|-----|-------------|--------------|--------|
| L1 | Fused AVX2 NR kernel (3 passes → 1) | 30–50% of NR part | Medium |
| L2 | Hoist RuntimeConfig::get() before NR loop | ~5–15 ns total/call | Trivial |
| L3 | TBB parallel_reduce for NR sums at n≥threshold | 4–8× at n=10000 | Medium |
| L4 | Pass buf[] (warm) to mad_select instead of xp[] | 5–15% at n≥100 | Trivial |

### OPT-L1 fused kernel math

Two accumulators in one pass:
```
u_i = (x_i - t) * half_inv_s
p_i = tanh(u_i)
sum_psi  += p_i                  [= Σ tanh(u)]
sum_dpsi += 1 - p_i²             [= Σ sech²(u) via fnmadd]
```
NR update: t += 2*s * sum_psi / sum_dpsi

fnmadd precision: -(p*p)+1 in one FMA rounding step — strictly better than
sub(one, mul(p,p)) which has two roundings and can underflow to exactly 0.

### sum_dpsi near-zero risk (Gemini F-12a)

sum_dpsi = Σ sech²(u_i) → 0 when ALL |u_i| >> 1 (data spread ≫ scale estimate).
This can happen when:
- Scale is severely underestimated at initialization
- Near-degenerate data (many tied values)

Current scalar code has NO guard. Need: `if (sum_dpsi < DBL_MIN) break;` in both
the fused kernel and the 3-pass fallback path.

---

## F-13: Gemini Technical Review (SESSION_ID: 85c6b412-8a72-4050-ba9c-6d377ccde621)

Confirmed 2026-03-21 on questions about robLoc fused AVX2 kernel design.

**Q1: Two-accumulator AVX2 (acc_psi += p, acc_dpsi += fnmadd(p,p,one)) numerically stable?**
→ YES. p = tanh(u) ∈ (-1,1), so 1-p² ∈ (0,1]. fnmadd preferred over sub(one,mul(p,p))
  as it performs the operation with a single rounding step, reducing drift.

**Q2: TBB structure for (psi, dpsi) reduction?**
→ Single `struct { double psi; double dpsi; }` with operator+= in one parallel_reduce.
  Avoids 2x data reads vs separate reductions. Preserves cache locality for xp[].

**Q3: mad_select(buf,n,med,dev) correctness?**
→ CORRECT and beneficial. MAD is permutation-invariant. buf[] is warm in L1/L2 after
  memcpy+median_select. Passing buf instead of xp is both correct and cache-friendly.

**Q4: sum_dpsi guard?**
→ NEEDED. sech²(u) can underflow to 0 for extreme u_i values. Guard:
  `if (sum_dpsi < std::numeric_limits<double>::min()) break;`

**Q5: RuntimeConfig hoist TOCTOU risk?**
→ ZERO RISK. CPUID features are invariant for process lifetime. Hoisting to static
  const bool is the standard pattern.

---

## F-14: robScale() Profile — n=500–1000 (WU-RS7, 2026-03-24)

### Diagnostic data (rob_scale_diag_impl, seed=42)

| n     | outer_iters | aitken_fires | rho_evals | converged |
|-------|-------------|--------------|-----------|-----------|
| 100   | 4           | 3            | 7         | 1         |
| 200   | 3           | 2            | 5         | 1         |
| 500   | 4           | 3            | 7         | 1         |
| 1000  | 3           | 2            | 5         | 1         |
| 5000  | 3           | 2            | 5         | 1         |
| 10000 | 3           | 2            | 5         | 1         |

**Finding**: rho_eval count (5–7) is constant w.r.t. n. Aitken Δ² acceleration is
working. Iteration count is NOT the bottleneck for large n.

### Timing (C_rob_scale_fast, bench::mark, seed=42)

| n    | median (µs) |
|------|-------------|
| 100  | ~3.8        |
| 200  | ~5.9        |
| 500  | ~13.8       |
| 1000 | ~20.5       |

Scaling is roughly O(n): 1000/500 ratio ≈ 1.5×, consistent with sort + linear
rho_sum dominating.

### Sort fraction

Compared C_rob_scale_fast on random vs pre-sorted input:
- Differences < 1 µs at all n (within timer noise ~500 ns).
- Conclusion: sort cost is **not separable** from rho_sum cost via this technique.
  The C++ introsort (for random data) and linear-scan variant (for sorted) are
  both fast enough that timer jitter swamps the difference at the measured n range.

### R-level overhead

- `robScale(x)` vs `C_rob_scale_fast(x)`: ~1–2 µs difference (R dispatch, `anyNA`,
  `length` checks). Verified via Rprof: 96% of time in "robscale::robScale" (opaque C
  block at 10 ms sampling interval — all C++ time is below profiler resolution).

### Implications for future work

1. **rho_sum dominates** the non-sort portion at large n. With 5–7 evals × n elements
   using the 8-wide AVX2 kernel, further gains require either (a) reducing rho_evals
   (Aitken already at minimum), or (b) improving kernel throughput (e.g., 16-wide
   AVX-512, out of scope for CRAN target).
2. **Sort elimination (OPT-9)** is the next measurable opportunity: for the ensemble
   path, resamples are already sorted before entering the robScale path. A
   `rob_scale_sorted` entry point that skips the sort step could save 20–40% at
   n=500–1000 within the ensemble (WU-RS9).
3. **No bias-correction opportunity** identified from this profile. OPT-8 is a
   correctness/accuracy change, not a speed change.
