# CRAN submission comments for robscale 0.5.2

## Changes in 0.5.2

### Bug fix

`scale_robust(x, method = "qn")` with `n >= threshold` previously returned
`gmd(x)` silently. The `auto_switch` guard now only fires for
`method = "ensemble"`; named methods are always dispatched as requested.

### New features

- `scale_robust(..., ci = TRUE)` supports bootstrap confidence intervals for
  individual named methods via `boot_method = "bca"`, `"percentile"`, or
  `"parametric"`. A new C++ bootstrap kernel (`cpp_single_estimator_ci_bounds`)
  runs the resampling on the single requested estimator.
- `boot_method = "analytical"` is accepted explicitly by `scale_robust()`.
- `print.robscale_ci` now shows the CI method in the output header.

## Test environments

- Arch Linux (x86_64), R 4.5.x, GCC 15.1.0
- macOS (ARM64), R 4.5.x, Apple Clang (rhub)
- Fedora Linux (x86_64), R-devel, GCC (rhub)
- Windows (x86_64), R-devel (rhub)

## R CMD check results

0 errors | 0 warnings | 1 note

- **Non-portable compilation flag** (`-march=native`): Used only when the
  building user's compiler supports it, detected at configure time. CRAN
  binary builds do not use this flag; it is present only in local
  `src/Makevars`. The configure script writes a portable fallback for
  environments where `-march=native` is not supported.

- **GNU make**: `SystemRequirements` field lists GNU make (needed for
  `$(shell)` in Makevars). Already declared in DESCRIPTION.

## DOI URLs

Several reference DOIs (Rousseeuw & Croux 1993, Aitken 1926, Steffensen 1933)
return HTTP 403 during automated URL checks because the publishers
(Taylor & Francis, JSTOR) block crawler requests. The URLs resolve correctly
in a browser. They are correct and permanent DOIs.

## Method references

Key references for the methods implemented:

- Rousseeuw, P.J. & Croux, C. (1993). Alternatives to the Median Absolute
  Deviation. *Journal of the American Statistical Association*, 88, 1273-1283.
  doi:10.1080/01621459.1993.10476408
- Aitken, A.C. (1926). On Bernoulli's numerical solution of algebraic equations.
  *Proceedings of the Royal Society of Edinburgh*, 46, 289-305.
  doi:10.2307/2333958
- Steffensen, J.F. (1933). Remarks on iteration. *Skandinavisk Aktuarietidskrift*,
  16, 64-72. doi:10.1080/03461238.1933.10419209
