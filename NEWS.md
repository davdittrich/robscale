# robscale 0.5.2

## Bug fixes

* `scale_robust(x, method = "qn")` with `n >= threshold` now correctly returns
  the Qn estimator instead of silently returning `gmd(x)`. The `auto_switch`
  mechanism now only applies to `method = "ensemble"`; named methods are always
  dispatched as requested.

## New features

* `scale_robust(..., ci = TRUE)` now supports bootstrap confidence intervals for
  individual named methods via `boot_method = "bca"`, `"percentile"`, or
  `"parametric"`. The new C++ function `cpp_single_estimator_ci_bounds` runs
  the bootstrap on the single requested estimator.

* `boot_method = "analytical"` is now accepted explicitly by `scale_robust()`.
  It selects the chi-squared interval for `sd` and the ARE-based normal
  approximation for all other named methods.

* `print.robscale_ci` now shows the CI method in the output header
  (e.g., `95% CI (analytical): [0.6123, 1.1456]`).

## Portability

* **AVX2 detection guarded for non-GCC/Clang compilers** (WU-PORT-1):
  `__builtin_cpu_supports("avx2")` is now wrapped in
  `#if defined(__GNUC__) || defined(__clang__)`. MSVC and ICC builds
  previously hit an undefined-identifier error here; on unsupported
  compilers, SIMD capability defaults conservatively to `None`.

* **`hardware_concurrency()` zero return clamped** (WU-PORT-1):
  `std::thread::hardware_concurrency()` returns 0 on platforms that
  cannot determine the thread count. The result is clamped to
  `std::max(..., 1u)`, preventing a zero-divisor in the TBB
  grain-size calculation on such platforms.

* **sysfs L2-cache parse hardened** (WU-PORT-1): The Linux sysfs
  integer read now wraps `std::stoul` in `try { } catch (...)`.
  Malformed or out-of-range values fall back to
  `L2_FALLBACK_KB * 1024` and `num_logical_cores = 1` instead of
  propagating an uncaught exception to the R session.

* **`mad.cpp` micro-buffer tied to `ROBSCALE_MICRO_BUFFER_SIZE`**
  (WU-PORT-2): Both `mad_impl_auto` and `mad_impl_center` previously
  hard-coded `buf_micro[64]` and `n <= 64`. Both now reference
  `buf_micro[ROBSCALE_MICRO_BUFFER_SIZE]` and
  `n <= static_cast<int>(ROBSCALE_MICRO_BUFFER_SIZE)`. Changing the
  package-wide constant in `robscale_config.h` now propagates to both
  MAD paths without manual adjustment.

## Performance

* **Ensemble bootstrap singleton overhead eliminated** (WU-PERF-1):
  `ensemble_one_replicate` called
  `robscale::qnsn::RuntimeConfig::get().sort_boost_threshold` on
  every bootstrap replicate. The call is replaced with the
  compile-time constant `ROBSCALE_SORT_BOOST_THRESHOLD`, removing the
  static-local singleton access from the hot path. Gate results:
  ~59% speedup at n = 1,000--5,000 (stable-zone ratio ≈ 0.41).

* **AVX2 diff-fill enabled in TBB Qn refinement path** (WU-PERF-2):
  The `qn_fill_diffs_avx2` kernel previously ran only in the
  brute-force path. It now dispatches inside the TBB `parallel_for`
  lambda when `use_avx2 && len >= 4`; the `use_avx2` flag is hoisted
  before the parallel region to avoid repeated CPUID lookups inside
  the lambda. Gate results: ~11% speedup at n = 1,000--10,000
  (stable-zone ratio ≈ 0.89).

* **AVX-512 tanh dispatch**: On glibc ≥ 2.35 systems where libmvec provides
  `_ZGVeN8v_tanh`, the bulk tanh kernel processes 8 doubles per iteration
  instead of 4. CPUID dispatch ensures the 8-wide path runs only on AVX-512F
  hardware; other systems fall through to AVX2 with zero overhead.

* **GMD AVX2 FMA kernel**: The Gini mean difference weighted sum uses a shared
  `gmd_weighted_sum` kernel with `_mm256_fmadd_pd` for n ≥ 8 on AVX2
  hardware, replacing three duplicate scalar loops.

* **Estimator-major bootstrap layout**: `boot_results` is stored as
  7 × n_boot (estimator-major) so the mean/variance reduction pass reads each
  estimator's replicates contiguously (stride-1 instead of stride-7).

## Internal

* **Dead `tmp` parameter removed from `rob_scale_compute`**
  (WU-CLEAN-1): The function accepted an unused
  `double* ROBSCALE_RESTRICT tmp` pointer left over from a planned
  optimization that was never implemented. Removed from the
  definition, its forward declaration, and all three call sites.

* **`STACK_SIZE` unified to `ROBSCALE_SN_STACK_THRESHOLD`**
  (WU-CLEAN-1): Four independent `constexpr int STACK_SIZE = 2048`
  definitions in `gmd.cpp`, `iqr.cpp`, `mad.cpp`, and
  `sn_estimator.cpp` are replaced with the package-wide constant in
  `robscale_config.h`.

* **`iqr_select_and_interp` helper extracted; dead `buf2` removed**
  (WU-CLEAN-2): The pdqselect + max-scan + Type-7 interpolation block
  was duplicated inline across the Q1 and Q3 selection paths in
  `iqr()`. Extracting it as `static inline iqr_select_and_interp` in
  `estimators_internal.h` eliminated the duplication and enabled
  compiler CSE across both call sites. Side effect: ~75% speedup at
  n = 1,000 (ratio ≈ 0.25). The now-unused `buf2` parameter is
  removed from `iqr()` and its call site in `ensemble.cpp`.

* **`sn_inner_serial` template extracted; Sn inner-loop triplicated no
  more** (WU-CLEAN-3): The Sn inner loop appeared identically in
  `SnWorker::operator()`, the Tier-1 serial fallback, and the Tier-2
  large-$n$ branch. All three now delegate to
  `sn_inner_serial<T>`, tagged `__attribute__((always_inline))`. The
  attribute is load-bearing: without it, GCC 15 declines to inline
  across all three call sites at `-O2`, producing a ~5× regression at
  n = 1,000.

* **`ROBSCALE_SORT_MEDIAN_THRESHOLD` renamed
  `ROBSCALE_SORT_NETWORK_THRESHOLD`** (WU-CLEAN-4): The macro guards
  the sort-network dispatch path, not median calculation. The old name
  was misleading. Renamed across all six files; no value change.

* **Documentation additions** (WU-CLEAN-4): `ensemble_one_replicate`
  receives a full doc comment covering its parameters, the XorShift32
  PRNG bit-shifts, and the sorted-once design pattern.
  `discover_linux`, `discover_macos`, and `discover_windows` in
  `qnsn_hardware_info.h` receive brief function-level comments. The
  `diag.cpp` section labels are replaced with semantic names
  ("M-scale iteration diagnostics", "Median crossover benchmark
  helpers").

# robscale 0.5.1

## Performance

* **`iqr_scaled()` dead stack eliminated** (OPT-I1+I2): The entry frame previously
  allocated both a 32 KB `buf_stack[4096]` and a 2 KB `buf_micro[256]` unconditionally,
  even for n=2. `buf_stack` is now declared only inside a `ROBSCALE_NOINLINE` large-n
  helper (`n > 256`); `buf_stack` size reduced 4096 → 2048 (IQR needs one copy).
  Stack at small n reduced from ~34 KB to ~2 KB.

* **`ROBSCALE_RESTRICT` on `iqr_impl_large` and `interp_q7`** (OPT-I5): Pointer
  parameters carry `ROBSCALE_RESTRICT`, allowing the compiler to generate stronger
  SIMD code in the selection and scan loops.

* **`iqr_scaled()` n≤16 sort fast path** (OPT-I6): For n≤16, two pdqselect calls
  plus two min-scan passes are replaced by a single `small_sort` (branch-free sorting
  network) followed by direct index reads. Q1/Q3 are O(1) after the sort.

* **Symmetric Q1 selection for n>16 with frac1>0** (OPT-I3): When the Q1 quantile
  requires interpolation (`(n-1) % 4 ≠ 0`, ~75% of n values), the lower-partition
  scan is reduced from O(0.75n) to O(0.25n). Instead of finding the minimum of
  `buf[lo1+1..n-1]`, pdqselect places `lo1+1` exactly and a max-scan over
  `buf[0..lo1]` recovers Q1. Applies to both the micro path (17≤n≤256) and the
  large-n NOINLINE helper (n>257). Gate results: n=100 −4.3%, n=1000 −3.3%.

* **`ROBSCALE_RESTRICT` on `estimators_internal::iqr()`** (OPT-I8): Pointer
  parameters of the inline `iqr()` function (used in the ensemble compute path) carry
  `ROBSCALE_RESTRICT`, allowing alias-free code generation in the pdqselect path.

* **`iqr_sorted()` added** (OPT-I7): New `robscale::internal::iqr_sorted(sorted_x, n)`
  for pre-sorted data — O(1) Type 7 quantile reads, no copy, no selection. Completes
  the sorted-variant family alongside `sn_sorted()` and `qn_sorted()`.

* **`mad_scaled()` small-n stack frame reduced 16 KB → 512 B** (OPT-M1): `mad_impl_auto`
  and `mad_impl_center` now route `n > 64` through `ROBSCALE_NOINLINE` helpers.
  The 16 KB `buf_stack[2048]` is declared only in those helpers, never in the entry frame.

* **`#pragma omp simd` on all deviation loops** (OPT-M2): All absolute-deviation loops
  in `mad_impl_auto`, `mad_impl_center`, and their large-n helpers are SIMD-annotated
  via the established `ROBSCALE_HAS_OMP_SIMD` guard pattern.

* **`ROBSCALE_RESTRICT` on NOINLINE helper parameters** (OPT-M3): `mad_impl_auto_large`
  and `mad_impl_center_large` carry `ROBSCALE_RESTRICT` on their input pointers,
  allowing the compiler to generate stronger SIMD code.

* **`bulk_abs_diff` / `bulk_abs_diff_inplace` shared SIMD kernels** (OPT-M4): Two
  inline SIMD kernels added to `src/robust_core.h`. All six deviation loops in
  `mad.cpp` are replaced with kernel calls. n=1000 shows ~2–3.5% improvement.

* **`mad_from_data` uses `bulk_abs_diff_inplace`** (OPT-M5): The deviation loop in
  `estimators_internal.h::mad_from_data` (used by the ensemble compute path) is
  replaced with the shared kernel.

* **Ensemble MAD block uses `bulk_abs_diff`** (OPT-M6): The per-replicate deviation
  loop in `ensemble_one_replicate` is replaced with `bulk_abs_diff`.

* **V-shaped deviation median in ensemble MAD block** (OPT-M7): Since `resample` is
  sorted, the deviation array is V-shaped. `vshaped_mad` finds the median of the two
  sorted halves via O(log n) binary search (`kth_of_two_sorted`), replacing the O(n)
  `median_select` call and eliminating the `work2` write-read cycle.

# robscale 0.5.0

## Performance

* **`gmd()` small-n stack frame reduced 17 KB → 1 KB** (OPT-G1): `gmd_impl` now
  routes `n ≤ 128` through an `ROBSCALE_NOINLINE` helper that allocates only a
  128-double arena. The large `buf_stack[2048]` is declared only in the `n > 128`
  path, eliminating 16 KB of unnecessary stack allocation on every small-sample
  call.

* **Ensemble bootstrap sort upgraded to `boost::float_sort`** (OPT-G4): The
  per-replicate sort in `ensemble_one_replicate` now dispatches to
  `boost::sort::spreadsort::float_sort` (serial radix sort) for `n` above the
  runtime-configured `sort_boost_threshold`. Expected speedup ~1.5–2× for
  bootstrap samples of size ≥ 512. `tbb::parallel_sort` is not used here to
  avoid nested-TBB oversubscription inside the outer `tbb::parallel_for`.

* **`gmd()` internal sort path uses `optimized_sort`** (OPT-G3): `internal::gmd()`
  in `estimators_internal.h` now uses `robscale::qnsn::optimized_sort` (tiered
  `std::sort` → `boost::float_sort` → `tbb::parallel_sort`) instead of plain
  `std::sort`. Call sites (`compute_all_estimators` at lines 216 and 330 of
  `ensemble.cpp`) are serial, so the optional TBB dispatch is safe.

* **GMD kernel scale factor precomputed** (OPT-G5): `constant * 2 / (n*(n-1))`
  is now computed once before the accumulation loop across all three call sites
  (`gmd_sorted`, `internal::gmd()`, ensemble inline block), reducing per-element
  work to a single multiply.

* **`ROBSCALE_RESTRICT` annotations on `gmd_sorted` and `gmd_core`** (OPT-G6):
  No-aliasing hints allow the compiler to generate tighter SIMD code for the
  accumulation loop.

---

# robscale 0.4.0

## Bug fixes

* **Aitken Δ² step-reduction guard removed** (`robScale()`): The previous
  guard rejected a valid extrapolation whenever the Aitken jump exceeded the
  two-step pair displacement (`|cand - s2| < |s2 - s0|`). For convergence
  rates $r \gtrsim 0.73$ this condition is never satisfied, so the
  acceleration never fired and the M-scale iteration required up to 80
  `rho_sum` evaluations per call at small $n$. Removing the guard — keeping
  only the minimal safety `candidate > 0` — cuts typical `rho_evals` from
  50–80 down to 7–12 for the previously-slow cases (n=5 seed=47/247,
  n=7 seed=49/149). The case n=7 seed=49, which previously hit `maxit` without
  converging, now converges in fewer than 40 evaluations.

## Performance

* **Fused AVX2 kernel threshold lowered from n≥8 to n≥4** (`robScale()`):
  The `use_fused` flag in `rob_scale_compute` now activates at $n \geq 4$
  instead of $n \geq 8$. Inputs of length 4, 5, 6, and 7 now use the
  single-pass fused AVX2 kernel rather than the three-pass scalar fallback,
  eliminating the remaining per-call overhead for very small samples.

* **Fused single-pass NR kernel for `robLoc()`** (OPT-L1): `rob_loc_nr_step_avx2`
  collapses the three-pass Newton--Raphson loop (scale residuals, `bulk_tanh`,
  accumulate) into a single AVX2 pass. Two accumulators advance in lockstep:
  `acc_psi` via `addpd` and `acc_p2` via `fmadd` ($p_i^2 = \tanh^2(u_i)$).
  The derivative sum follows from $\sum \text{sech}^2(u_i) = n - \sum \tanh^2(u_i)$,
  avoiding a second transcendental evaluation. A degenerate-scale guard
  (`sum_dpsi < DBL_MIN`) prevents NaN on near-constant inputs.

* **RuntimeConfig hoist for `robLoc()`** (OPT-L2): The SIMD dispatch flag
  was re-read from thread-local storage on every Newton--Raphson iteration.
  CPUID features are invariant over process lifetime; the flag is now hoisted
  once before the loop. Added `bulk_tanh_dispatched()` to `robust_core.h` to
  accept the pre-hoisted flag.

* **Warm-cache MAD for `robLoc()`** (OPT-L4): `rob_loc_core` copies input
  into `buf[]` and selects the median in place. MAD is permutation-invariant,
  so `mad_select` now reads from `buf[]`---warm in L1/L2 cache---rather than
  the cold input array.

* **TBB parallel NR for `robLoc()`** (OPT-L3): For $n \geq \max(4096, L_2/32)$,
  `rob_loc_parallel_compute` partitions the array across TBB threads via
  `parallel_reduce`. Each chunk calls `rob_loc_nr_step_avx2`; an `NRAccum`
  struct accumulates partial sums with `operator+=`.

* **Combined speedup:** These four optimizations yield **2.6--4.1x** speedup
  over `revss` for `robLoc()` at $n = 100$ to $10{,}000$ on x86\_64 with AVX2
  (benchmark ratios: 0.347 at $n = 100$, 0.243 at $n = 10{,}000$).

* **Sorting-network threshold corrected for `robScale()`** (OPT-A): Lowered
  `ROBSCALE_SORT_MEDIAN_THRESHOLD` from 64 to 16. The median-net comparator
  count grows as $O(n^{1.5})$ ($n = 64$: 337 swaps; $n = 16$: 46 swaps);
  Floyd--Rivest is $O(n)$, crossing over at $n \approx 16$. Using sorting
  networks up to $n = 64$ incurred ~3.5× unnecessary overhead per median call.
  This change eliminates the $n = 500$--$1{,}000$ regression in `robScale()`:
  benchmarks show ratios 0.201 ($n = 500$) and 0.174 ($n = 1{,}000$), versus
  1.126 and 1.122 before the fix (a 5--6× reversal).

* **PLT elimination for `median_net<T>`** (OPT-B): Added
  `__attribute__((visibility("hidden")))` to the `median_net<T>` template
  definition in `sort_net.h`. On Linux `-fPIC`, this eliminates PLT
  indirection for all intra-DSO calls, replacing the W (weak) symbol with a
  direct $t$ (local) binding.

* **`adm()` hot-buffer pass** (OPT-A): `adm_impl_auto()` now passes the
  already-copied working buffer to `adm_core()` rather than the cold SEXP
  pointer. The $n$ doubles remain in L1/L2 cache from the initial copy
  through the deviation sum, avoiding an additional cold read.

* **`adm()` stack-frame split** (OPT-B): Extracted `adm_large_n()` as a
  `ROBSCALE_NOINLINE` helper. This reduces the `adm_impl_auto()` stack frame
  from 8,256 to 624 bytes (objdump verified), preventing the 32 KB large-$n$
  buffer from penalizing small calls.

## Internal

* **Benchmark seed pool widened to 7 seeds** (`benchmarks/run_benchmarks.R`):
  `BENCH_SEEDS` increased from 3 to 7 seeds (spaced 50 apart: 42, 92, 142,
  192, 242, 292, 342). At small $n$ ($n \leq 16$) the previous 3-seed pool
  had high variance because 2 of 3 seeds could land on the previously-slow
  Aitken path, biasing the pooled-median timing toward the slow cluster.

* Added `rob_loc_scalar_impl`, `rob_loc_has_parallel`, and `rob_loc_serial_impl`
  diagnostic exports for TDD correctness gates.
* Added compile-time diagnostic helpers (`src/diag.cpp`) for sorting-network
  threshold validation; not user-facing.

# robscale 0.3.0

## Performance

* **Fused median-then-MAD**: `mad_scaled()` now computes absolute deviations
  in-place on the median selection buffer, eliminating the second scratch array.
  Memory drops from $2n$ to $n$ doubles. The same fusion applies to `robScale()`
  and the ensemble's internal MAD and M-scale paths, reducing total ensemble
  workspace from $3n$ to $2n$ doubles. Benchmarks show 7--15% speedup at medium
  $n$ on x86\_64, with no regressions on ARM64 or at small $n$.

* **Noinline small-$n$ dispatch for `robScale()` and `robLoc()`**: extracted
  shared core logic (`rob_scale_core`, `rob_loc_core`) and separated the
  small-$n$ path ($n \leq 64$) into `ROBSCALE_NOINLINE` helpers with a 1 KB
  stack frame. The large-$n$ path retains its 32 KB buffer independently.
  Benchmarks showed no measurable timing improvement (the R-to-C++ `.Call()`
  boundary dominates at small $n$), but the refactoring clarifies the code
  structure and prevents the compiler from penalizing small calls with the
  large-$n$ stack reservation.

## Internal

* Added `ROBSCALE_NOINLINE` portability macro to `robscale_config.h`
  (`__attribute__((noinline))` on GCC/Clang, `__declspec(noinline)` on MSVC).
* `estimators_internal.h`: `mad_from_data()` and `rob_scale()` signatures
  reduced from three to two buffer arguments.
* `ensemble.cpp`: workspace allocation reduced from $3n$ to $2n$ doubles.

# robscale 0.2.2

* Tests: Adapted cross-validation tests for compatibility with `revss` v3.0.0,
  which changed `robScale()` and `robLoc()` defaults to use bias-corrected "AA"
  constants. The `adm()` cross-check remains active for all `revss` versions.
* Tests: Added frozen golden reference tests for `robScale()` and `robLoc()`,
  decoupled from upstream `revss`. Catches regressions in our own algorithm
  independent of `revss` version changes.

# robscale 0.2.1

* Performance: Refined runtime SIMD dispatch for `qn()` and `sn()` kernels,
  achieving optimal vectorization across x86_64 (AVX2/FMA) and ARM64 (NEON).
* New feature: `scale_robust()` provides a unified dispatcher that automatically
  selects between a variance-weighted ensemble of 7 estimators (for small
  samples) and the Gini Mean Difference (GMD) for larger samples (n >= 20).
* New estimators: Added `gmd()`, `iqr_scaled()`, `mad_scaled()`, and `sd_c4()`
  to the public API.
* Confidence Intervals: Added `ci = TRUE` to all scale estimators, providing
  analytical intervals (GMD, MAD, Sn, Qn, bias-corrected SD) or bootstrap-based
  intervals (ensemble).
* Performance: Replaced `stats::mad()` and `stats::IQR()` with optimized C++
  implementations (`mad_scaled()`, `iqr_scaled()`) using O(n) selection via the
  pdqselect algorithm.
* Documentation: Expanded README with detailed benchmarking vs. `robustbase`,
  `Hmisc`, `GiniDistance`, and `collapse`.
* Citations: Consolidated and updated all package citations in `inst/CITATION`
  and documentation.

# robscale 0.1.6

* Performance: Extended optimal sorting networks from n <= 8 to n <= 16 using
  Dobbelaere's verified optimal networks. Cross-platform benchmarking confirmed
  2-4x speedups over `std::sort` for n = 9-16 on both ARM64 (Apple Silicon)
  and x86_64 (AMD Zen 3).
* Bug fix: Corrected the n = 7 sorting network comparator sequence (was
  producing incorrect sort order for certain inputs).
* Bug fix: Restored ADM consistency constant to sqrt(pi/2) = 1.2533 (was
  incorrectly changed to 1.3926).
* Bug fix: Restored small-sample fallback logic in `robScale()` that was
  inadvertently removed, causing incorrect results when MAD collapses for
  n < 4.

# robscale 0.1.5

* CRAN: Removed `#pragma GCC diagnostic ignored "-Wdeprecated-volatile"` from
  `src/qn_estimator.cpp`, `src/sn_estimator.cpp`, and `src/qnsn_sort_utils.h`
  to resolve the "pragmas suppressing diagnostics" NOTE requested by CRAN
  maintainers.
* Windows: Added `configure.win` and `cleanup.win` to suppress spurious
  win-builder warnings about missing Windows configuration.
* CI: Added GitHub Actions R-CMD-check workflow with multi-platform matrix
  (Windows, macOS, Ubuntu) as independent evidence of Windows compilation.
* Portability: Fixed uninitialized NEON register in `qnsn_kernels.h` that
  Clang 17 on macOS promoted to a compilation WARNING.
* Build: Switched `RcppParallel::LdFlags()` to the canonical
  `RcppParallel::RcppParallelLibs()` in `Makevars.win` and `Makevars.in`.

# robscale 0.1.4

* Portability: Fixed macOS (Apple Silicon) compilation error by switching to
  the more portable `<Accelerate/Accelerate.h>` include and resolving a
  symbol conflict with R's `COMPLEX` type.
* Portability: Restored explicit `CXX_STD = CXX17` in Makevars (required by
  CRAN policy for packages using C++17 features).
* Portability: Removed non-portable GCC pragmas to silence warnings on Clang
  (win-builder Debian).
* CRAN: Retained `GNU make` in `SystemRequirements` (`$(shell)` is a GNU
  extension used for RcppParallel linking).
* CRAN: Added `inst/WORDLIST` to whitelist technical terms from `aspell`.
* CRAN: Quoted 'Qn' and 'Sn' in `DESCRIPTION` to satisfy metadata checks.
* Build: Resolved duplicate `-ltbb` flag from `Makevars` which caused warnings
  on some macOS configurations.

# robscale 0.1.3

* Added high-efficiency scale estimators: `qn()` and `sn()` (Rousseeuw &
  Croux, 1993) with specialized sorting network kernels and cache-aware
  parallelization.
* Added `robustbase` to `Suggests` for reference and benchmarking.
* Added documentation examples for all five exported functions.
* Performance optimizations: `configure` script auto-detects SIMD
  capabilities (AVX2/FMA on x86_64, NEON on ARM64) without requiring
  any environment variables.

# robscale 0.1.1

* Consistent NA handling: `adm()` now raises an error when `na.rm = FALSE` and
  NAs are present, matching `robLoc()` and `robScale()` behavior.
* Replaced manual `new[]`/`delete[]` with `std::unique_ptr` for RAII memory
  safety in all C++ estimator functions.
* Removed unused scoring-iteration constants (`PSI_NORM_CONST`,
  `INV_PSI_NORM_CONST`) from `robust_core.h`.
* Hoisted loop-invariant `half_inv_s` above the Newton--Raphson loop in
  `rob_loc.cpp`.
* Added Mermaid algorithm flowcharts for `robLoc` and `robScale` to README.
* Expanded test suite: explicit-center tests for `adm()`, `maxit` edge cases
  for `robLoc()` and `robScale()`, sorting network verification for n = 2--8.
* Added `cph` role to `Authors@R` in DESCRIPTION.

# robscale 0.1.0

Initial release.

* Three functions: `adm()`, `robLoc()`, and `robScale()` implementing
  the robust location and scale M-estimators of Rousseeuw & Verboven (2002)
  for very small samples.
* API-compatible drop-in replacement for the `revss` package.
* C++17 implementation via Rcpp with:
  * Newton--Raphson iteration for location (quadratic convergence).
  * Algebraic `tanh(x/2)` identity for the logistic psi function with
    platform-vectorized bulk evaluation (Apple Accelerate on macOS,
    OpenMP SIMD hints on Linux).
  * Floyd--Rivest O(n) selection for median computation.
  * Optimal sorting networks for n <= 8.
  * Stack-allocated arena buffers (zero heap allocation for n <= 512).
* Numerical equivalence with `revss` verified across 5,400 systematic
  comparisons (tolerance: sqrt(.Machine$double.eps)).
