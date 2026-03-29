# CRAN submission comments for robscale 0.5.2

## Changes in 0.5.2

### Bug fix

`scale_robust(x, method = "qn")` with `n >= threshold` previously returned
`gmd(x)` silently. The `auto_switch` guard now only fires for
`method = "ensemble"`; named methods are always dispatched as requested.

### New features

- `scale_robust(..., ci = TRUE)` supports bootstrap confidence intervals for
  individual named methods via `boot_method = "bca"`, `"percentile"`, or
  `"parametric"`.
- `boot_method = "analytical"` is accepted explicitly by `scale_robust()`.
- `print.robscale_ci` now shows the CI method in the output header.

### Performance

- Input validation moved from R to C++ (`validate_finite`): eliminates the
  per-call `any(!is.finite(x))` heap allocation that consumed 20–65% of
  call time for simple estimators.
- `Rcpp::RNGScope` removed from all exports (`rng = false`): no estimator
  function uses R's RNG.
- Newton–Raphson convergence criterion for `robLoc` scaled by location
  magnitude, avoiding 80 wasted iterations on large-valued data.
- BCa jackknife estimator-index mapping corrected for `gmd` and `sd_c4`.

## Test environments

- Arch Linux (x86_64), R 4.5.0, GCC 15.1.0 (local)
- Ubuntu 24.04 (x86_64), R-release (GitHub Actions): 0 errors, 0 warnings
- macOS latest (ARM64), R-release (GitHub Actions): 0 errors, 0 warnings
- Windows (x86_64), R-release (GitHub Actions): 0 errors, 0 warnings

## R CMD check results

0 errors | 0 warnings | 1 note (local only)

The single NOTE (`-march=native`) appears only in local builds where the
user's `~/.R/Makevars` sets this flag. CRAN binary builds use their own
compiler flags and will not see this NOTE. The package's `configure` script
writes portable defaults.

## DOI URLs

Three DOIs (Rousseeuw & Croux 1993, Shamos 1976, David & Nagaraja 2003)
return HTTP 403 during automated URL checks because the publishers block
crawler User-Agents. All resolve correctly in a browser. They are correct
and permanent DOIs.

## Method references

- Rousseeuw, P.J. & Croux, C. (1993). Alternatives to the Median Absolute
  Deviation. *Journal of the American Statistical Association*, 88, 1273–1283.
  doi:10.1080/01621459.1993.10476408
- Rousseeuw, P.J. & Verboven, S. (2002). Robust estimation in very small
  samples. *Computational Statistics & Data Analysis*, 40(4), 741–758.
  doi:10.1016/S0167-9473(02)00078-6
