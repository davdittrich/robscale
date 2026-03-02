## Resubmission

This is a resubmission (v0.1.1). In the previous submission (v0.1.0), CRAN
requested that 'Rcpp' be single-quoted in the DESCRIPTION. This has been
addressed. Additional changes in v0.1.1: consistent NA error handling across
all exported functions, RAII memory management in C++ code, expanded test
suite, and added `cph` role to Authors@R.

## Test environments

* macOS Tahoe 26.3 (aarch64), R 4.5.2, Apple clang 17.0.0

## R CMD check results

0 errors | 0 warnings | 1 note

* NOTE: New submission

## Notes

* This is the first release of the `robscale` package.
* The package reimplements the robust estimators from the CRAN package `revss`
  (listed in Suggests) in C++17 via Rcpp. The API is intentionally identical to
  `revss` to serve as a drop-in replacement.
* Numerical equivalence with `revss` is verified by 5,400 cross-comparisons in
  the test suite (`inst/tinytest/test_cross_check.R`).
* The package uses a `configure` script to detect `-fopenmp-simd` compiler
  support and the macOS Accelerate framework. A `cleanup` script removes the
  generated `src/Makevars`. On Windows, `src/Makevars.win` provides a static
  fallback (C++17 only, no platform-specific flags).
