# CRAN submission comments for robscale 0.1.4

This submission (0.1.4) addresses all points raised in pre-test logs for version 
0.1.2/0.1.3 across Windows, Linux (Debian/Clang), and macOS platforms.

## Changes in 0.1.4
- Fixed macOS ARM64 compilation failure by switching to `<Accelerate/Accelerate.h>` and resolving a `COMPLEX` symbol conflict with R's internal types.
- Restored explicit `CXX_STD = CXX17` in Makevars (required by CRAN policy).
- Removed non-portable GCC-specific pragmas that caused warnings on Clang.
- Resolved duplicate `-ltbb` linking flag warning on macOS.
- Retained `GNU make` in `SystemRequirements` (required for `$(shell)` in Makevars).
- Added `inst/WORDLIST` to address technical term spelling notes.
- Quoted 'Qn' and 'Sn' in metadata.

Version 0.1.2 introduces two major new features and several refinements:

- **New estimators**: Added `qn()` and `sn()` (Rousseeuw & Croux, 1993) with RcppParallel and cache-aware kernels.
- **Portability**: Fixed the `configure` script to ensure that hardware-specific flags (e.g., `-mavx2`, `-mfma`) are
  strictly **opt-in** via `ROBSCALE_FAST=1`, satisfying CRAN's portability requirements.
- **Documentation**: Added missing `@examples` for all exported functions.

## Test environments

- Arch Linux (x86_64), R 4.5.2, GCC 15.2.1

## R CMD check results

0 errors | 0 warnings | 0 notes

(Note: System-specific notes for missing `aspell` or `tidy` may appear in some environments but do not reflect package defects.)

## URL Checks

`urlchecker` flags three DOI URLs as `403 Forbidden`. These have been manually
verified to work correctly in a browser. The access denials are likely due to
bot-protection on the part of the publishers (Taylor & Francis, ACM, and JSTOR).

Affected DOI links:
- [doi:10.1080/01621459.1993.10476408](https://doi.org/10.1080/01621459.1993.10476408)
- [doi:10.1145/360680.360691](https://doi.org/10.1145/360680.360691)
- [doi:10.2307/2332448](https://doi.org/10.2307/2332448)

## Notes

- The package reimplements and extends the robust estimators from the CRAN
  packages `revss` and `robustbase` (both listed in Suggests) in C++17 via Rcpp.
- Numerical equivalence with `revss` is verified by 5,400 cross-comparisons in
  the test suite.
