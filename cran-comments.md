# CRAN submission comments for robscale 0.5.3

## Changes in 0.5.3

### Build fixes (resubmission from 0.5.2)

- Added `target("avx2,fma")` / `target("avx512f")` attributes to
  `extern "C"` declarations of libmvec/SLEEF vectorized tanh functions.
  Clang 21 (CRAN Debian) requires target attributes on `__m256d`/`__m512d`
  function signatures even in declarations, not just definitions.
- Unified TBB preprocessor guards across all source files: `qn_estimator.cpp`,
  `sn_estimator.cpp`, `ensemble.cpp`, and `worker_compat.h` now use the
  combined `(ROBSCALE_HAS_SYSTEM_TBB || USE_DIRECT_TBB)` pattern, matching
  `rob_scale.cpp` and `rob_loc.cpp`. This restores direct TBB parallelism
  for Qn, Sn, and ensemble bootstrap on system-oneTBB builds (CRAN Linux).
- Eliminated unused-variable warning (`use_avx2` in `qn_estimator.cpp`).

## Test environments

- Arch Linux (x86_64), R 4.5.0, GCC 15.1.0 (local)
- Ubuntu 24.04 (x86_64), R-release (GitHub Actions): 0 errors, 0 warnings
- macOS latest (ARM64), R-release (GitHub Actions): 0 errors, 0 warnings
- Windows (x86_64), R-release (GitHub Actions): 0 errors, 0 warnings

## R CMD check results

0 errors | 0 warnings | 1 note (local only)

The single NOTE (`-march=native`) appears only in local builds where the
user's `~/.R/Makevars` sets this flag. CRAN binary builds use their own
compiler flags and will not see this NOTE.

## DOI URLs

Three DOIs (Rousseeuw & Croux 1993, Shamos 1976, David & Nagaraja 2003)
return HTTP 403 during automated URL checks because the publishers block
crawler User-Agents. All resolve correctly in a browser.

## Method references

- Rousseeuw, P.J. & Croux, C. (1993). Alternatives to the Median Absolute
  Deviation. *Journal of the American Statistical Association*, 88, 1273–1283.
  doi:10.1080/01621459.1993.10476408
- Rousseeuw, P.J. & Verboven, S. (2002). Robust estimation in very small
  samples. *Computational Statistics & Data Analysis*, 40(4), 741–758.
  doi:10.1016/S0167-9473(02)00078-6
