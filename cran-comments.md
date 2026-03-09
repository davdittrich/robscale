# CRAN submission comments for robscale 0.1.5

## Changes in 0.1.5

- Added `configure.win` and `cleanup.win` to suppress win-builder's
  "this package has a configure script / It probably needs manual configuration"
  warning.
- Added GitHub Actions CI workflow with a multi-platform matrix including
  Windows (R-release and R-devel), macOS, and Ubuntu.
- Switched `RcppParallel::LdFlags()` to the canonical
  `RcppParallel::RcppParallelLibs()` in `Makevars.win` and `Makevars.in`.

## Win-builder incoming_pretest note

The **Windows** incoming_pretest for 0.1.4 failed with an installation error.
The root cause is a toolchain mismatch on the pretest server: R was compiled with
MinGW g++ (GCC 14.3.0), but the build invokes Cygwin's g++ (GCC 15.2.0). When
Cygwin's g++ resolves standard headers via MinGW's `_mingw_stdarg.h`, that
header raises `#error Only Win32 target is supported!`. The error occurs at the
very first `#include <Rcpp.h>` before any package-specific code is reached — no
source-level change can fix this.

The **Debian/Clang** incoming_pretest passed cleanly (0 errors, 0 warnings,
2 NOTEs).

As independent evidence that robscale compiles and passes `R CMD check` on
Windows, we now include a GitHub Actions CI workflow that tests on
`windows-latest` with both R-release and R-devel (proper Rtools45/MinGW
toolchain). Results are visible at:
https://github.com/davdittrich/robscale/actions

## Test environments

- Arch Linux (x86_64), R 4.5.2, GCC 15.2.1
- GitHub Actions: Windows Server 2022 (x86_64), R-release + R-devel
- GitHub Actions: macOS (ARM64), R-release
- GitHub Actions: Ubuntu 22.04 (x86_64), R-release + R-devel

## R CMD check results

0 errors | 0 warnings | 0 notes

(System-specific notes for missing `aspell` or `tidy` may appear in some
environments but do not reflect package defects.)

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
