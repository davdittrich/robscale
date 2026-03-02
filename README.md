# robscale

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18828607.svg)](https://doi.org/10.5281/zenodo.18828607)

Fast robust estimation of location and scale in very small samples.

## Overview

Classical estimators of location (mean) and scale (standard deviation) can be
destroyed by a single outlier when the sample size is small. In experimental
sciences---analytical chemistry, materials testing, clinical trials---measurements
are often repeated only 3 to 8 times per condition, making this vulnerability a
practical concern rather than a theoretical one.

Rousseeuw & Verboven (2002) developed M-estimators using a logistic psi function
that remain well-behaved at these sample sizes. Their estimators decouple location
and scale estimation to avoid the positive-feedback instability of Huber's
Proposal 2, use the MAD as a fixed auxiliary scale, and fall back gracefully when
the sample is too small even for iteration. The
[revss](https://CRAN.R-project.org/package=revss) package by Avraham Adler
provides a faithful pure-R implementation of these estimators.

`robscale` is a C++17/Rcpp reimplementation of the same three estimators---`adm`,
`robLoc`, and `robScale`---with an API-compatible interface. It produces
numerically identical results to `revss` (verified across 5,400 systematic
comparisons) while running 6--33x faster, with the largest gains in the target
regime of $n \leq 20$.

Beyond compiled execution, `robscale` departs from `revss` in three algorithmic
choices: it replaces the scoring fixed-point iteration for location with true
Newton--Raphson, it exploits the algebraic identity $\psi_{\log}(x) = \tanh(x/2)$
to enable platform-vectorized transcendental evaluation, and it uses the
Floyd--Rivest selection algorithm for $O(n)$ median computation instead of
$O(n \log n)$ sorting.

## Installation

```r
# Requires a C++17 compiler (GCC >= 7, Clang >= 5, MSVC >= 19.14)
remotes::install_github("davdittrich/robscale")
```

## Motivating example

```r
library(robscale)

x <- c(2.0, 3.1, 2.7, 2.9, 3.3)   # clean measurements

mean(x)                              # 2.8000
sd(x)                                # 0.5000
robLoc(x)                            # 2.8471
robScale(x)                          # 0.3837

x[5] <- 100                         # recording error

mean(x)                              # 22.1400  -- destroyed
sd(x)                                # 43.5270  -- destroyed
robLoc(x)                            # 2.9184   -- stable
robScale(x)                          # 0.4729   -- stable
```

The classical estimators are pulled toward the outlier. The robust estimators
barely move.

## API reference

| Function | Purpose | Key arguments |
|:---|:---|:---|
| `adm(x)` | Average distance to the median | `center`, `constant` |
| `robLoc(x)` | Robust M-estimate of location | `scale` |
| `robScale(x)` | Robust M-estimate of scale | `loc`, `implbound` |

All three functions accept `na.rm` (default `FALSE`).

### `adm(x, center, constant = 1.2533141373155001, na.rm = FALSE)`

Computes the mean absolute deviation from the median, scaled by a consistency
constant for asymptotic normality under the Gaussian model:

$$\text{ADM}(x) = C \cdot \frac{1}{n}\sum_{i=1}^{n} |x_i - \text{med}(x)|$$

where $C = \sqrt{\pi/2} \approx 1.2533$ (Nair, 1947). When `center` is
supplied, it replaces the median.

The ADM is not itself robust against outliers (explosion breakdown $1/n$), but
it is highly resistant to implosion (implosion breakdown $(n{-}1)/n$). It serves
as the fallback scale estimator when the MAD collapses to zero.

```r
adm(c(1, 2, 3, 5, 7, 8))
adm(c(1, 2, 3, 5, 7, 8), constant = 1)   # without consistency correction
```

### `robLoc(x, scale = NULL, na.rm = FALSE, maxit = 80L, tol = sqrt(.Machine$double.eps))`

M-estimator for location using Newton--Raphson iteration with the logistic psi
function (Rousseeuw & Verboven 2002, Eq. 21). Starting value:
$T^{(0)} = \text{median}(x)$. Auxiliary scale: $S = \text{MAD}(x)$ (or the
user-supplied `scale`). See [Algorithmic innovations](#algorithmic-innovations)
for the iteration formula.

**Fallback logic:** When `scale` is unknown and $n < 4$, or when `scale` is
known and $n < 3$, the function returns `median(x)` without iteration.
Providing a known `scale` lowers the minimum sample size from 4 to 3 because the
MAD (which is unreliable at $n = 3$) is no longer needed.

```r
robLoc(c(1, 2, 3, 5, 7, 8))
robLoc(c(1, 2, 3), scale = 1.5)   # known scale enables n = 3
```

### `robScale(x, loc = NULL, implbound = 1e-4, na.rm = FALSE, maxit = 80L, tol = sqrt(.Machine$double.eps))`

M-estimator for scale using multiplicative iteration with the rho function---the
square of the logistic psi (Rousseeuw & Verboven 2002, Eq. 27):

$$S^{(k)} = S^{(k-1)} \cdot \sqrt{2 \cdot \frac{1}{n}\sum \psi_{\log}^2\!\left(\frac{x_i - T}{c \cdot S^{(k-1)}}\right)}$$

where $c = 0.37394112142347236$ and $T = \text{median}(x)$ is held fixed.
Starting value: $S^{(0)} = \text{MAD}(x)$.

**Fallback logic:** When $n$ is below the minimum for iteration (4 for unknown
location, 3 for known), the function checks for MAD implosion:

- If $\text{MAD}(x) \leq$ `implbound`: return `adm(x)` (implosion fallback)
- Otherwise: return `MAD(x)`

Providing a known `loc` centers the data at that value and uses the
median-distance-to-zero ($1.4826 \cdot \text{median}(|x_i - \mu|)$) as the
initial scale, lowering the minimum sample size from 4 to 3.

```r
robScale(c(1, 2, 3, 5, 7, 8))
robScale(c(1, 2, 3, 5, 7, 8), loc = 5)   # known location
```

## Algorithmic innovations

Both `robscale` and `revss` implement the estimators of Rousseeuw & Verboven
(2002). They solve the same estimating equations and produce the same numerical
results. The two packages differ in *how* the computation is carried out. This
section describes the algorithmic choices in `robscale` that diverge from the
`revss` implementation, and their performance consequences.

### 1. The tanh identity for the logistic psi function

The logistic psi function central to both estimators is:

$$\psi_{\log}(x) = \frac{e^x - 1}{e^x + 1}$$

`revss` evaluates this as `2 * plogis(u) - 1`, which calls R's `plogis`
(the logistic CDF, $1/(1 + e^{-x})$). The computation requires one call to
`exp()` followed by two arithmetic operations, plus the overhead of R's
vectorised dispatch, intermediate vector allocation, and garbage-collection
pressure.

The algebraic identity

$$\psi_{\log}(x) = \tanh(x/2)$$

is immediate from the definition of the hyperbolic tangent. `robscale` exploits
this identity to reduce $\psi_{\log}$ to a single `tanh` call. This is not
merely a cosmetic rewrite:

- **Branch elimination.** A direct implementation of $(e^x - 1)/(e^x + 1)$
  overflows for large $|x|$, requiring a sign-based branch to keep intermediate
  values bounded. The `tanh` function handles this internally with a single code
  path.

- **Platform vectorization.** Standard library `tanh` implementations are
  available in vectorized form on major platforms. On macOS, `robscale` calls
  Apple Accelerate's `vvtanh`, which processes the entire argument array in a
  single SIMD pass. On x86_64 Linux with GCC or Clang, `#pragma omp simd` hints
  allow the compiler to generate SIMD-vectorized `tanh` without requiring
  `-ffast-math` or the OpenMP runtime library. On Windows, the scalar
  `std::tanh` fallback still benefits from branch elimination.

### 2. Newton--Raphson iteration for location

`revss` iterates the location estimator using the scoring fixed-point iteration
(Rousseeuw & Verboven 2002, Eq. 21):

$$T^{(k+1)} = T^{(k)} + S \cdot \frac{\frac{1}{n}\sum \psi_{\log}\!\left(\frac{x_i - T^{(k)}}{S}\right)}{\alpha}$$

where $\alpha = \int \psi_{\log}'(u)\,d\Phi(u) \approx 0.4132$ is the
normalization constant. This is a fixed-point iteration with *linear*
convergence: each step reduces the error by a constant factor.

`robscale` instead applies Newton--Raphson to the estimating equation
$\sum \psi_{\log}((x_i - T)/S) = 0$, yielding:

$$T^{(k+1)} = T^{(k)} + \frac{2\,S\sum \psi\!\left(\frac{x_i - T^{(k)}}{2S}\right)}{\sum \left[1 - \psi^2\!\left(\frac{x_i - T^{(k)}}{2S}\right)\right]}$$

where $\psi(\cdot) = \tanh(\cdot)$. This follows from observing that the
derivative of the logistic psi satisfies $\psi_{\log}'(x) = 1 - \psi_{\log}^2(x)
= 1 - \tanh^2(x/2)$. Since $\tanh$ values have already been computed for the
numerator, the denominator requires only squaring and subtraction---no additional
transcendental function calls.

Newton--Raphson achieves *quadratic* convergence near the solution: the number
of correct digits approximately doubles per iteration. The practical effect is
a reduction from 4--8 iterations (scoring) to 3 iterations (Newton--Raphson)
for reaching the same tolerance of $\sqrt{\epsilon_{\text{mach}}} \approx 1.49
\times 10^{-8}$:

| $n$ | Scoring iterations | Newton--Raphson iterations |
|---:|---:|---:|
| 4 | 7 | 3 |
| 5 | 8 | 3 |
| 8 | 7 | 3 |
| 20 | 6 | 3 |
| 100 | 5 | 3 |

At small $n$, the scoring iteration count is higher because the starting value
(the median) can be far from the M-estimate in units of the auxiliary scale.
Newton--Raphson absorbs this gap in fewer steps.

### 3. Floyd--Rivest selection for median computation

Both the median and the MAD require computing a quantile---the median of the
data, and the median of the absolute deviations. `revss` uses R's `median()`
and `mad()`, which call `sort()` internally: an $O(n \log n)$ operation.

`robscale` uses the Floyd--Rivest selection algorithm (Floyd & Rivest, 1975), an
average-case $O(n)$ algorithm that partitions the array around the $k$th smallest
element without fully sorting it. For the even-$n$ case, a single linear scan
over the upper partition locates the $(k{+}1)$th element needed for averaging.

For $n \leq 8$, the selection step is replaced by a full sort using optimal
sorting networks (Knuth, TAOCP Vol. 3, Sec. 5.3.4). These are branchless
compare-and-swap sequences with the minimum number of comparisons for each $n$:

| $n$ | Comparators |
|---:|---:|
| 3 | 3 |
| 4 | 5 |
| 5 | 9 |
| 6 | 12 |
| 7 | 16 |
| 8 | 19 |

For $n > 8$, the Floyd--Rivest algorithm dispatches to `std::nth_element` for
partitions below 600 elements.

### 4. Arena allocation on the stack

Each estimator requires working arrays: a copy of the input (for destructive
selection), absolute deviations (for MAD), and a temporary buffer (for bulk
`tanh` arguments). `revss` allocates these as R vectors, incurring R's SEXPREC
header overhead and adding garbage-collection pressure.

`robscale` uses a stack-allocated arena of 512 doubles per segment. For
$n \leq 512$---which covers the target regime and far beyond---there is zero
heap allocation. For $n > 512$, the code falls back to `new[]`/`delete[]`.

### 5. Compile-time reciprocal constants

The constants $1/\alpha = 1/0.413241928283814$ and
$1/c = 1/0.37394112142347236$ are declared `constexpr`, allowing the compiler to
replace divisions in the iteration loop with multiplications. On ARM64, this
avoids the ~10-cycle `fdiv` instruction in favour of a ~3-cycle `fmul`.

### 6. Loop-invariant hoisting

Values that are constant across iterations---`inv_s = 1.0 / s`,
`half_inv_s = 0.5 / s`, `inv_n = 1.0 / n`---are computed once before the loop.
The `revss` implementation recomputes `(x - t) / s` as a fresh R vector each
iteration, traversing the interpreter for every vectorised operation.

## Benchmarks

Platform: R 4.5.2, Apple clang 17.0.0, macOS (Darwin 25.3.0), Apple Silicon
(aarch64). Median of 10,000 `microbenchmark` iterations per cell ($n \leq 100$);
2,000 iterations for $n > 100$.

| $n$ | Function | `revss` ($\mu$ s) | `robscale` ($\mu$ s) | Speedup |
|---:|:---|---:|---:|---:|
| 3 | `adm` | 8.86 | 0.86 | 10x |
| 3 | `robLoc` | 21.77 | 1.56 | 14x |
| 3 | `robScale` | 48.50 | 1.68 | 29x |
| 4 | `adm` | 10.37 | 0.86 | 12x |
| 4 | `robLoc` | 56.05 | 1.72 | 33x |
| 4 | `robScale` | 93.28 | 3.03 | 31x |
| 5 | `adm` | 8.53 | 0.82 | 10x |
| 5 | `robLoc` | 43.95 | 1.64 | 27x |
| 5 | `robScale` | 95.94 | 3.20 | 30x |
| 8 | `adm` | 10.33 | 0.86 | 12x |
| 8 | `robLoc` | 36.24 | 1.68 | 22x |
| 8 | `robScale` | 79.62 | 2.87 | 28x |
| 20 | `adm` | 10.58 | 0.90 | 12x |
| 20 | `robLoc` | 37.56 | 1.80 | 21x |
| 20 | `robScale` | 90.53 | 3.94 | 23x |
| 100 | `adm` | 11.07 | 1.03 | 11x |
| 100 | `robLoc` | 53.18 | 2.58 | 21x |
| 100 | `robScale` | 113.41 | 8.61 | 13x |
| 500 | `adm` | 14.74 | 1.76 | 8x |
| 500 | `robLoc` | 103.14 | 6.19 | 17x |
| 500 | `robScale` | 288.64 | 35.10 | 8x |
| 1000 | `adm` | 18.16 | 3.16 | 6x |
| 1000 | `robLoc` | 163.14 | 10.82 | 15x |
| 1000 | `robScale` | 484.78 | 65.56 | 7x |

**Interpretation.** Speedups are largest in the target regime ($n \leq 20$):
10--33x for individual functions. The `robLoc` speedup exceeds what compiled
execution alone would predict because the Newton--Raphson iteration converges in
3 steps where the scoring iteration requires 4--8 (see
[Newton--Raphson iteration](#2-newtonraphson-iteration-for-location)).

At larger $n$, the advantage narrows. The `adm` function has no iteration loop,
so its speedup (6--12x) reflects only the compiled loop and $O(n)$ median
selection. The `robScale` narrowing at large $n$ (7--8x) is expected: the
iteration count is the same in both implementations (multiplicative iteration is
retained), so the per-iteration cost of vectorized `tanh` vs interpreted
`plogis()` becomes the dominant factor, and both scale as $O(n)$.

The `revss` timings include R function-call dispatch, vectorised `plogis()`
allocation, and garbage-collection pressure. `robscale` eliminates all three.

## Numerical equivalence

The test suite (`inst/tinytest/test_cross_check.R`) compares `robscale` and
`revss` across 5,400 randomly generated inputs:

- Sample sizes $n = 3, 4, \ldots, 20$
- 100 replicates per $n$ (uniform on $[-100, 100]$, seed 42)
- All three functions: `adm`, `robLoc`, `robScale`
- Tolerance: $\sqrt{\epsilon_{\text{mach}}} \approx 1.49 \times 10^{-8}$

All 5,400 comparisons pass. The Newton--Raphson iteration converges to the same
fixed point as the scoring iteration---it solves the same estimating equation---so the results
differ only by rounding at the level of the convergence tolerance.

## Mathematical background

The estimators follow Rousseeuw & Verboven (2002). A brief summary of the key
definitions:

**Logistic psi function** (Eq. 23):

$$\psi_{\log}(x) = \frac{e^x - 1}{e^x + 1} = \tanh(x/2)$$

Bounded in $(-1, 1)$, smooth ($C^\infty$), strictly monotone. Boundedness
provides robustness; smoothness avoids the corner artifacts of Huber's psi at
small $n$.

**Decoupled estimation.** Location and scale are estimated separately with a
fixed auxiliary estimate, breaking the positive-feedback loop of Huber's
Proposal 2. `robLoc` fixes scale at $\text{MAD}(x)$; `robScale` fixes
location at $\text{median}(x)$.

**Rho function for scale** (Eq. 26):

$$\rho_{\log}(x) = \psi_{\log}^2(x / c)$$

where $c = 0.37394112142347236$ is the constant that yields 50% breakdown point.

**Key constants** (full double precision):

| Symbol | Value | Definition |
|:---|:---|:---|
| $\alpha$ | `0.413241928283814` | $\int \psi_{\log}'(u)\,d\Phi(u)$; scoring normalization constant |
| $c$ | `0.37394112142347236` | Solution to $\int \rho_{\log}(u)\,d\Phi(u) = 0.5$; scale rho constant |
| $C_{\text{ADM}}$ | `1.2533141373155001` | $\sqrt{\pi/2}$; ADM consistency constant |
| $C_{\text{MAD}}$ | `1.4826` | $1/\Phi^{-1}(3/4)$; MAD consistency constant |

## Relation to revss

This package reimplements the algorithms from the
[revss](https://CRAN.R-project.org/package=revss) package by Avraham Adler.
The API is intentionally identical: `adm()`, `robLoc()`, and `robScale()` accept
the same arguments and return the same values. Code that uses `revss` can switch
to `robscale` by changing only the `library()` call.

Users who do not need compiled performance---or who prefer a dependency-free
pure-R package---should use `revss` directly. It is mature, well-tested, and
available on CRAN.

## References

Rousseeuw, P.J. and Verboven, S. (2002). Robust estimation in very small
samples. *Computational Statistics & Data Analysis*, **40**(4), 741--758.
[doi:10.1016/S0167-9473(02)00078-6](https://doi.org/10.1016/S0167-9473(02)00078-6)

Floyd, R.W. and Rivest, R.L. (1975). Expected time bounds for selection.
*Communications of the ACM*, **18**(3), 165--172.
[doi:10.1145/360680.360691](https://doi.org/10.1145/360680.360691)

Nair, K.R. (1947). A Note on the Mean Deviation from the Median. *Biometrika*,
**34**(3/4), 360--362.
[doi:10.2307/2332448](https://doi.org/10.2307/2332448)

## Author

Dennis Alexis Valin Dittrich
([ORCID](https://orcid.org/0000-0002-4438-8276))

## License

MIT. Copyright 2026 Dennis Alexis Valin Dittrich.
