# CRAN submission comments for robscale 0.2.0

## Changes in 0.2.0

- New feature: `scale_robust()` unified dispatcher with variance-weighted
  ensemble and automatic switching to GMD for n >= 20.
- New estimators: `gmd()`, `iqr_scaled()`, `mad_scaled()`, and `sd_c4()` added
  to exported functions.
- Statistical inference: Added analytical and bootstrap confidence intervals
  to all scale estimators via a new `ci = TRUE` argument.
- Optimization: `mad_scaled()` and `iqr_scaled()` now leverage the `pdqselect`
  algorithm for O(n) selection, outperforming base R's O(n log n) sorting.
- Citations: Synchronized `inst/CITATION` and documentation with the new
  metadata and foundational references (Gini, 1912).

## Test environments

- Arch Linux (x86_64), R 4.5.3, GCC 15.2.1
- GitHub Actions: Windows Server 2022 (x86_64), R-release + R-devel
- GitHub Actions: macOS (ARM64), R-release
- GitHub Actions: Ubuntu 24.04 (x86_64), R-release + R-devel

## R CMD check results

0 errors | 0 warnings | 3 notes

- **CRAN incoming feasibility**: New submission; maintainer address confirmed.
- **Compilation flags used**: NOTE regarding non-portable flags (`-march=x86-64`,
  etc.) provided by the R core build environment for performance and security.
- **HTML math rendering**: NOTE regarding V8 unavailability on the test system;
  docs render correctly where V8 is present.

## URL Checks

`urlchecker` flags five DOI URLs as `403 Forbidden`. These have been manually
verified to work correctly in a browser. The access denials are likely due to
bot-protection on the part of the publishers (Taylor & Francis, ACM, and JSTOR).

Affected DOI links:

- [doi:10.1080/00401706.1962.10490022](https://doi.org/10.1080/00401706.1962.10490022)
- [doi:10.1080/01621459.1993.10476408](https://doi.org/10.1080/01621459.1993.10476408)
- [doi:10.1145/360680.360691](https://doi.org/10.1145/360680.360691)
- [doi:10.2307/2332448](https://doi.org/10.2307/2332448)
- [doi:10.2307/2333958](https://doi.org/10.2307/2333958)

## Notes

- The package reimplements and extends the robust estimators from the CRAN
  packages `revss` and `robustbase` (both listed in Suggests) in C++17 via Rcpp.
- Numerical equivalence with `revss` is verified by 5,400 cross-comparisons in
  the test suite.
