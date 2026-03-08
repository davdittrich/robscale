# robscale 0.1.2

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
* Added high-efficiency scale estimators: `qn()` and `sn()` (Rousseeuw &
  Croux, 1993) with specialized sorting network kernels and cache-aware
  parallelization.
* Added documentation examples for all five exported functions.
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
