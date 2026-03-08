# CRAN submission comments for robscale 0.1.2

This is a new submission for version 0.1.2. The previous version (0.1.1) is
currently on CRAN. 

Version 0.1.2 introduces two major new features and several refinements:
- **New estimators**: Added `qn()` and `sn()` (Rousseeuw & Croux, 1993) with RcppParallel and cache-aware kernels.
- **Portability**: Fixed the `configure` script to ensure that hardware-specific flags (e.g., `-mavx2`, `-mfma`) are strictly **opt-in** via `ROBSCALE_FAST=1`, satisfying CRAN's portability requirements.
- **Documentation**: Added missing `@examples` for all exported functions.

## Test environments

* macOS Tahoe 26.3 (aarch64), R 4.5.2, Apple clang 17.0.0

## R CMD check results

0 errors | 0 warnings | 0 notes

(Note: System-specific notes for missing `aspell` or `tidy` may appear in some environments but do not reflect package defects.)

## Notes

* The package reimplements and extends the robust estimators from the CRAN package `revss` (listed in Suggests) in C++17 via Rcpp.
* Numerical equivalence with `revss` is verified by 5,400 cross-comparisons in the test suite.
