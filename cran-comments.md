# CRAN submission comments for robscale 0.5.4

## Changes in 0.5.4

### Performance

- SIMD median selection networks (AVX2): 23–58% faster at n = 8–32 for
  `robScale`, `robLoc`, `MAD`, and `ADM`.
- Raised sort network threshold (16 → 56) and median network threshold
  (16 → 36), benefiting `Qn`, `Sn`, `GMD`, and `IQR`.
- Earlier ensemble parallelism: TBB fires for all bootstrap calls with
  ≥ 2 cores and ≥ 16 replicates. Ensemble is 4–8× faster for n = 3–49.
- Restored AVX2 auto-vectorization in `adm_core` (6–9% regression fix).

### Architecture

- Dropped AVX-512 dispatch (8-wide tanh). The tanh backend hierarchy is
  now: Apple Accelerate → glibc libmvec AVX2 → SLEEF AVX2 → scalar.
- Runtime AVX2 dispatch for ADM via explicit `use_avx2` flag.
- L2 cache plausibility guard (< 64 KB → fallback to 256 KB).

### Build fixes

- Fixed sed delimiter collision in `configure` (| → !).
- Added `RcppParallel.h` before system TBB headers for correct include
  guard ordering.
- Added `libtbb.dylib` check for macOS RcppParallel TBB detection.
- Removed `TBB` from `SystemRequirements` (provided by RcppParallel).
- Explicit `<cmath>` include for portability.

## Test environments

- Arch Linux (x86_64), R 4.5.3, GCC 15.2.1 (local)
- Ubuntu 24.04 (x86_64), R-release, GCC (GitHub Actions): PASS
- Ubuntu 24.04 (x86_64), R-devel, GCC (GitHub Actions): PASS
- macOS latest (ARM64), R-release (GitHub Actions): PASS
- Windows (x86_64), R-release (GitHub Actions): PASS
- Windows (x86_64), R-devel (GitHub Actions): FAIL — renv/untar
  infrastructure bug in dependency install step; package never compiled

## R CMD check results

0 errors | 0 warnings | 0 notes

The `-march=native` NOTE seen in local builds is from the user's
`~/.R/Makevars`, not the package. CRAN binary builds use their own
compiler flags and will not see this NOTE.

## DOI URLs

Three DOIs (Rousseeuw & Croux 1993, Shamos 1976, David & Nagaraja 2003)
return HTTP 403 during automated URL checks because the publishers block
crawler User-Agents. All resolve correctly in a browser.

## Downstream dependencies

None (new package).

## Method references

- Rousseeuw, P.J. & Croux, C. (1993). Alternatives to the Median Absolute
  Deviation. *Journal of the American Statistical Association*, 88, 1273–1283.
  doi:10.1080/01621459.1993.10476408
- Rousseeuw, P.J. & Verboven, S. (2002). Robust estimation in very small
  samples. *Computational Statistics & Data Analysis*, 40(4), 741–758.
  doi:10.1016/S0167-9473(02)00078-6
