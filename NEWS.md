# robscale 0.2.0

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
