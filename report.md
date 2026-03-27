# Adversarial Documentation Audit Report

Ground truth: C++ source code (executable lines only, not comments).
Reviewed: README.qmd, vignettes/robscale-intro.Rmd, R/*.R roxygen, NEWS.md,
DESCRIPTION, cran-comments.md.

---

## README.qmd

### WRONG-R1. robScale API section describes the wrong iteration algorithm

**Location**: README.qmd lines 423–428
**Claim**: "M-estimator for scale using multiplicative iteration with the rho
function" with formula $S^{(k+1)} = S^{(k)} \sqrt{2n^{-1}\sum\rho}$.
**Reality**: `src/rob_scale.cpp:202–243` (`nr_scale_compute`) implements
Newton–Raphson: `delta_s = s * numer / denom` where
`numer = sum_tanh2/n − 0.5`, `denom = 2/n * sum_u_tanh_sech2`, then
`s += delta_s`. The multiplicative formula appears ONLY as a degenerate-
denominator fallback (line 228–229). The performance section at line 627
correctly says "Newton–Raphson." The API section contradicts it.

### WRONG-R2. Micro-buffer is 128 doubles, not 64

**Location**: README.qmd line 620
**Claim**: "A 64-double micro-buffer covers the small-sample regime"
**Reality**: `robscale_config.h:80`: `#define ROBSCALE_MICRO_BUFFER_SIZE 128`.
Confirmed by `rob_scale.cpp:357`: `double arena[ROBSCALE_MICRO_BUFFER_SIZE]`
with comment "128 doubles = 1 KB". Off by 2×.

### WRONG-R3. Cross-check tolerance and scope are stale

**Location**: README.qmd line 831
**Claim**: "100% pass rate at tolerance √ε_mach ≈ 1.49×10⁻⁸" across "5,400
inputs (n=3..20, 100 reps, all three functions adm/robLoc/robScale)."
**Reality**: `test-cross-check.R:3` uses `tol <- 1e-4` (6,700× more lenient).
Lines 8–17 show robLoc/robScale comparisons are conditional on
`packageVersion("revss") < "3.0.0"`. For current revss (≥3.0.0), only adm
is compared. Effective scope: ~1,800 inputs (adm only), tolerance 1e-4.

### MISLEADING-R4. Denominator formula applies only to robLoc, not both

**Location**: README.qmd line 630
**Claim**: "the denominator Σ(1 − ψ_i²) requires only squaring and
subtraction — no additional transcendental calls"
**Reality**: True for robLoc (`rob_loc.cpp:186`: `sum_dpsi += 1.0 - p*p`).
False for robScale, whose denominator is `2/n * Σ(u_i·tanh(u_i)·sech²(u_i))`
(`rob_scale.cpp:49–51`), requiring multiplication by u_i and tanh(u_i). The
sentence says "Both M-estimators" use this pattern; only robLoc does.

### IMPRECISE-R5. "halving memory reads" understates the improvement

**Location**: README.qmd lines 610–611
**Claim**: "halving memory reads relative to the standard three-pass approach"
**Reality**: The old approach (`rob_loc_scalar_impl:316–321`) does 3 passes.
The fused kernel does 1 pass. 3→1 is a 3× reduction, not "halving" (2→1).

### IMPRECISE-R6. ARE_ADM, ARE_SN, ARE_QN code constants round down

**Location**: README.qmd line 295 (API table)
**Claim**: ADM 88.3%, Sn 58.2%, Qn 82.3%
**Reality**: `robust_core.h` uses `ARE_ADM=0.88`, `ARE_SN=0.58`, `ARE_QN=0.82`.
These rounded constants affect analytical CIs. The README displays the standard
literature values; the code rounds slightly.

### IMPRECISE-R7. XorShift32 seed is `r + 12345`, not `r`

**Location**: README.qmd line 577
**Claim**: "XorShift32 PRNG seeded by replicate index"
**Reality**: `ensemble.cpp:68`: `XorShift32 rng(static_cast<uint32_t>(r + 12345))`.
Seed is r + 12345, not r. Determinism is preserved; the description simplifies.

---

## vignettes/robscale-intro.Rmd

### WRONG-V1. robScale flowchart describes old multiplicative iteration

**Location**: robscale-intro.Rmd lines 255–281 (Mermaid flowchart)
**Claim**: Flowchart shows `v = sqrt(2 * sum_rho / n)`, `s *= v`, convergence
test `|v − 1| ≤ tol`.
**Reality**: `rob_scale.cpp:202–243` implements NR: `delta_s = s * numer/denom`,
`s += delta_s`, convergence test `|delta_s|/s ≤ tol`. The multiplicative step
is only the degenerate fallback (line 229). The prose at lines 406–414 correctly
describes NR, but the flowchart contradicts it.

### WRONG-V2. Section header says "Multiplicative iteration"

**Location**: robscale-intro.Rmd lines 246–247, 255
**Claim**: "M-estimator for scale via multiplicative iteration"; flowchart title
"Multiplicative iteration for `robScale()`"
**Reality**: Code uses NR. The prose in the same section (lines 406–414) correctly
says NR. Internal self-contradiction.

### WRONG-V3. Ensemble weighting formula uses σ̂_j² instead of bootstrap variance

**Location**: robscale-intro.Rmd lines 126, 503
**Claim**: $\hat\sigma = \frac{\sum \hat\sigma_j / \hat\sigma_j^2}{\sum 1/\hat\sigma_j^2}$
— denominator is the square of each estimator's point value.
**Reality**: `ensemble.cpp:246–264` computes `weights[j] = (1/vars[j]) / weight_sum`
where `vars[j]` is the **bootstrap variance** (line 241:
`vars[j] = sq_sum / (count - 1.0)`), not the square of the point estimate.
The formula should use Var_boot(σ̂_j), not σ̂_j².

### WRONG-V4. Micro-buffer is 128 doubles, not 64

**Location**: robscale-intro.Rmd line 363
**Claim**: "A 64-double micro-buffer covers n ≤ 64"
**Reality**: `robscale_config.h:80`: `ROBSCALE_MICRO_BUFFER_SIZE 128`. The n≤64
dispatch threshold is a separate concept (it controls which entry point is
called, not the buffer size). The buffer is 128 doubles = 1 KB.

### WRONG-V5. ADM breakdown point is 1/n, not 25%

**Location**: robscale-intro.Rmd line 94
**Claim**: `adm(x)` breakdown = 25%
**Reality**: The ADM's explosion breakdown point is approximately 1/n — a single
extreme outlier inflates the mean deviation. The README.qmd API table (line 296)
correctly says `$1/n$`. The vignette's 25% is inconsistent with both the README
and standard results.

### MISLEADING-V6. SIMD hierarchy presented as 6 runtime levels

**Location**: robscale-intro.Rmd lines 309–320
**Claim**: 6-level hierarchy: Accelerate > libmvec-8 > libmvec-4 > SLEEF >
OpenMP > scalar, presented as a runtime cascade.
**Reality**: `robust_core.h:41–49` shows libmvec-4 and SLEEF are **mutually
exclusive at compile time** — they define the same macro `ROBSCALE_TANH4_AVX2`
via `#elif`. A binary contains at most one of them. At runtime, `bulk_tanh()`
dispatches at most 4 levels deep. The vignette should say "5 compile-time
backends, up to 4 runtime dispatch levels."

### MISLEADING-V7. Ensemble parallel claim unreachable with defaults

**Location**: robscale-intro.Rmd lines 419–420
**Claim**: "`scale_robust()` parallelizes the n_boot bootstrap iterations
across cores."
**Reality**: `ensemble.cpp:198` gates parallel dispatch on
`int64_t(n) * nboot >= ENSEMBLE_PARALLEL_THRESHOLD` where the threshold is
10,000 (line 31). Default ensemble mode: n < 20, n_boot = 200, max product
= 19 × 200 = 3,800 < 10,000. The parallel path is **unreachable** with
default settings. It fires only for large n or large n_boot.

### MISLEADING-V8. "2,048-double buffer" conflates threshold with buffer size

**Location**: robscale-intro.Rmd lines 363–364
**Claim**: "a 2,048-double buffer handles n ≤ 2,048"
**Reality**: `rob_scale.cpp:376–377`: `SCALE_STACK_SIZE = 2048` and
`double buf_stack[SCALE_STACK_SIZE * 2]` = **4,096 doubles** (32 KB).
The 2,048 is the n threshold, not the buffer size.

### IMPRECISE-V9. robLoc iteration count table claims exactly 3 for all n

**Location**: robscale-intro.Rmd lines 384–391
**Claim**: NR column shows 3 iterations for n = 4, 8, 20, 100.
**Reality**: `rob_loc.cpp:143` says "typically 2–4"; line 425 says
"≤ 4 evaluations on typical data." The table's blanket "3" is an idealization.

### IMPRECISE-V10. "follows original R&V 2002 definitions" overstates

**Location**: robscale-intro.Rmd lines 540–542
**Claim**: "robscale follows the original @rousseeuw2002 definitions"
**Reality**: The constants (c, consistency factors) match R&V 2002. The
iteration algorithm does not — R&V 2002 Eq. 27 is the multiplicative
fixed-point; robscale uses NR. Should say "follows R&V 2002 for the
estimating equations and constants; the solver uses Newton–Raphson."

### IMPRECISE-V11. robScale ARE "~50–70%" vs code constant 0.55

**Location**: robscale-intro.Rmd line 99
**Claim**: "~50–70%¹" with footnote "ARE varies with sample size"
**Reality**: `robust_core.h:62`: `ARE_ROBSCALE = 0.55`. The code uses a fixed
55%; the "50–70%" range is unsourced.

### STALE-V12. Alpha constant 0.413241928283814 is unused

**Location**: robscale-intro.Rmd line 510
**Claim**: α = 0.413241928283814, described as "scoring normalization constant"
**Reality**: Grep for `0.413241928283814` returns zero matches in src/. The
constant is a relic of the scoring iteration that was replaced by NR.

### STALE-V13. "1/α" constexpr claim — no such constant exists

**Location**: robscale-intro.Rmd line 435
**Claim**: "Reciprocal constants (1/α, 1/c) are constexpr"
**Reality**: `robscale_config.h:123` has `INV_RHO_SCALE_CONST = 1/c`. There is
no 1/α anywhere in the source.

---

## R/robScale.R (→ man/robScale.Rd)

### WRONG-RS1. Iteration formula shows Eq. 27 (multiplicative), code does NR

**Location**: R/robScale.R lines 51–54
**Claim**: Shows $S^{(k+1)} = S^{(k)}\sqrt{\frac{2}{n}\sum\psi^2}$
**Reality**: `rob_scale.cpp:224–238`: `delta_s = s * numer / denom; s += delta_s`.
The multiplicative formula is only the degenerate fallback (line 229).

### MISLEADING-RS2. "linearly converging" then "Quadratic convergence"

**Location**: R/robScale.R lines 63–70
**Claim**: Line 63 says "linearly converging fixed-point"; line 69 says
"Quadratic convergence yields 3–4 iterations."
**Reality**: The code implements NR (quadratic). The "linearly converging"
sentence describes the OLD algorithm. Both cannot be true simultaneously.

---

## R/robLoc.R (→ man/robLoc.Rd)

No errors found. SIMD dispatch, NR iteration, fallback behavior, and the
fused AVX2 kernel description all verified against `rob_loc.cpp`.

---

## R/qn.R (→ man/qn.Rd)

### WRONG-QN1. "Johnson-Mizoguchi selection" does not exist in codebase

**Location**: R/qn.R line 41
**Claim**: "A specialized Johnson-Mizoguchi selection algorithm for
medium-sized datasets."
**Reality**: Grep for "Johnson" and "Mizoguchi" returns zero matches in src/.
The actual algorithm is the Croux–Rousseeuw weighted-median refinement
(`qn_refinement_kernel` in `qn_estimator.cpp`).

---

## R/adm.R (→ man/adm.Rd)

### IMPRECISE-ADM1. "C++17 introselect" — actual algorithms differ

**Location**: R/adm.R line 51
**Claim**: "a C++17 introselect algorithm for larger datasets"
**Reality**: The code uses Floyd–Rivest selection (which wraps
`std::nth_element` for n < 600) and pdqselect above the L2-derived
threshold. Neither is "introselect" specifically, though `std::nth_element`
may be implemented as introselect by the standard library.

---

## R/gmd.R (→ man/gmd.Rd)

### IMPRECISE-GMD1. "29.3% breakdown point" uncited

**Location**: R/gmd.R line 33
**Claim**: "breakdown point is approximately 29.3%"
**Reality**: Standard result for GMD is 1 − 1/√2 ≈ 29.3%. Correct value but
no reference cited for it, and it's not derivable from code.

---

## NEWS.md

### No factual errors found.

Performance entries (AVX-512 dispatch, GMD FMA kernel, estimator-major layout)
match the current source. Historical entries (v0.4.0 Aitken guard, etc.) are
correct descriptions of what happened at that version.

---

## DESCRIPTION

### No factual errors found.

The Description field's claims ("SIMD vectorization", "TBB parallelism",
"fused single-buffer algorithms", "variance-weighted ensemble") all verified.

---

## cran-comments.md

### No factual errors found after recent edits.

DOI note updated to reference only Rousseeuw & Croux 1993. Method references
list Rousseeuw & Croux 1993 and Rousseeuw & Verboven 2002.

---

## Summary

| Severity | Count | IDs |
|----------|------:|-----|
| **WRONG** | 9 | R1, R2, R3, V1, V2, V3, V4, V5, RS1 |
| **MISLEADING** | 6 | R4, V6, V7, V8, RS2, QN1 |
| **IMPRECISE** | 8 | R5, R6, R7, V9, V10, V11, ADM1, GMD1 |
| **STALE** | 2 | V12, V13 |
| **OK** | 40+ | All other verified claims |

### Critical pattern

The dominant failure mode is the **multiplicative-vs-NR split**: the
algorithm was changed from multiplicative fixed-point to Newton–Raphson,
but the documentation was only partially updated. The flowchart,
formulas, section headers, and some prose still describe the old algorithm
while the convergence paragraphs describe the new one. This creates
internal self-contradictions within the same document.
