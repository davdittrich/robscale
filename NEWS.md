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
* Performance optimizations: Fixed `configure` script to ensure that
  hardware-specific SIMD flags (e.g., `-mavx2`) are strictly opt-in
  via `ROBSCALE_FAST=1`.

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
  - Newton--Raphson iteration for location (quadratic convergence).
  - Algebraic `tanh(x/2)` identity for the logistic psi function with
    platform-vectorized bulk evaluation (Apple Accelerate on macOS,
    OpenMP SIMD hints on Linux).
  - Floyd--Rivest O(n) selection for median computation.
  - Optimal sorting networks for n <= 8.
  - Stack-allocated arena buffers (zero heap allocation for n <= 512).
* Numerical equivalence with `revss` verified across 5,400 systematic
  comparisons (tolerance: sqrt(.Machine$double.eps)).
