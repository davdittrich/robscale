# robscale: A Comprehensive Toolkit for Robust Scale Estimation


[![](https://zenodo.org/badge/DOI/10.5281/zenodo.18828607.svg)](https://doi.org/10.5281/zenodo.18828607)

## Overview

Outliers compromise the reliability of classical estimators even in
moderate samples. A single recording error can destroy the standard
deviation, yet the most widely used robust implementations in R—`revss`
for small-sample M-estimation and `robustbase` for the $Q_n$ and $S_n$
scale estimators—carry significant computational overhead that blocks
real-time genomics pipelines, high-frequency finance, and other
time-critical workflows.

`robscale` v0.2.0 provides a comprehensive suite of 11 exported
functions spanning the full robustness–efficiency spectrum: from the
non-robust but maximally efficient bias-corrected standard deviation
(`sd_c4`, 100% ARE) through the near-optimal Gini mean difference
(`gmd`, 98% ARE, 29.3% breakdown), to maximum-breakdown estimators
(`qn`, `sn`, `mad_scaled`, all 50% breakdown). The unified
`scale_robust()` dispatcher automatically selects the optimal strategy—a
variance-weighted bootstrap ensemble of all 7 estimators for small
samples, with an automatic switch to the GMD for large samples.

All estimators are implemented as C++17 kernels with SIMD-accelerated
transcendental functions, Newton–Raphson iteration, and parallelized
$O(n \log n)$ algorithms. Against `revss`, the package achieves
**10.6–31.0x** speedups for the small-sample M-estimators. Against
`robustbase`, it achieves **1.8–6.9x** for $S_n$ and **1.7–6.8x** for
$Q_n$—with gains peaking near **6.9x** at $n = 10^7$ as TBB parallelism
reduces the computational bottleneck for massive datasets. The new
estimators (`gmd`, `iqr_scaled`, `mad_scaled`) outperform their base R
and CRAN counterparts by **2.4–17.5x** (GMD vs `Hmisc`), **1.0–30.2x**
(IQR vs `stats::IQR`), and **2.0–21.4x** (MAD vs `stats::mad`).

## Installation

``` r
install.packages("robscale")

# Development version:
# install.packages("remotes")
# remotes::install_github("davdittrich/robscale")
```

## Motivating example

``` r
library(robscale)

x <- c(2.0, 3.1, 2.7, 2.9, 3.3)   # clean measurements

# Recommended entry point: scale_robust()
scale_robust(x)                      # ensemble (small n < 20)
scale_robust(rnorm(50))              # auto-switches to gmd (n >= 20)

# Individual estimators span the robustness-efficiency frontier
gmd(x)                               # 98% ARE, 29.3% breakdown
qn(x)                                # 82.3% ARE, 50% breakdown
mad_scaled(x)                        # 36.8% ARE, 50% breakdown

# Outlier resistance
x[5] <- 100                          # recording error
sd(x)                                # destroyed
gmd(x)                               # stable (29.3% breakdown)
qn(x)                                # stable (50% breakdown)
mad_scaled(x)                        # stable (50% breakdown)
```

For small samples ($n < 20$), `scale_robust()` combines all 7 estimators
via a variance-weighted bootstrap ensemble, giving each estimator
influence proportional to its precision. For larger samples, it
automatically switches to the GMD, which achieves 98% ARE at negligible
computational cost.

## API reference

<div id="tbl-api-scale">

Table 1

**Table 1: Scale estimators** (sorted by decreasing ARE)

| Function | Purpose | ARE | Breakdown | Reference |
|:---|:---|:---|:---|:---|
| `sd_c4(x)` | Bias-corrected standard deviation | **100%** | 0% | Welford (1962) |
| `gmd(x)` | Gini mean difference | **98%** | 29.3% | Gini (1912); Nair (1936) |
| `adm(x)` | Average deviation from median | **88.3%** | $1/n$ | Nair (1947) |
| `qn(x)` | $Q_n$ scale estimator | **82.3%** | 50% | Rousseeuw & Croux (1993) |
| `sn(x)` | $S_n$ scale estimator | **58.2%** | 50% | Rousseeuw & Croux (1993) |
| `robScale(x)` | M-estimate of scale | **55.0%** | 50% | Rousseeuw & Verboven (2002) |
| `iqr_scaled(x)` | Scaled interquartile range | **37%** | 25% | Bickel & Lehmann (1976) |
| `mad_scaled(x)` | Scaled median absolute deviation | **36.8%** | 50% | Rousseeuw & Croux (1993) |

</div>

<div id="tbl-api-dispatch">

Table 2

**Table 2: Dispatcher and utilities**

| Function | Purpose |
|:---|:---|
| `scale_robust(x)` | Unified dispatcher: ensemble for small $n$, auto-switches to GMD for large $n$ |
| `get_consistency_constant(method, n)` | Returns the consistency constant or finite-sample correction for a given estimator |

</div>

Additionally, `robLoc(x)` provides an M-estimate of location (98.4% ARE,
50% breakdown; Rousseeuw & Verboven 2002). All functions accept `na.rm`
(default `FALSE`).

### `sd_c4(x, na.rm = FALSE)`

Computes the sample standard deviation corrected by $c_4(n)$ to remove
the small-sample bias of the square-root estimator:

$$\hat\sigma = \frac{s}{c_4(n)} = \frac{s}{\sqrt{2/(n{-}1)} \cdot \Gamma(n/2) / \Gamma((n{-}1)/2)}$$

where $s$ is the usual sample standard deviation. Uses Welford’s online
algorithm for numerically stable variance computation, avoiding
catastrophic cancellation. This is a non-robust estimator (0% breakdown)
with 100% ARE by construction—it serves as the efficiency anchor in the
ensemble.

``` r
sd_c4(c(1, 2, 3, 5, 7, 8))
```

### `gmd(x, constant = 0.886226925452758, na.rm = FALSE)`

Computes the Gini mean difference (Gini, 1912), scaled by a consistency
constant for asymptotic normality under the Gaussian model:

$$\text{GMD}(x) = C \cdot \frac{2}{n(n{-}1)}\sum_{i=1}^{n} (2i - n - 1)\, x_{(i)}$$

where $x_{(1)} \le \ldots \le x_{(n)}$ are the order statistics and
$C = \sqrt{\pi}/2 \approx 0.8862$. The computation requires a full sort
($O(n \log n)$), with sorting networks applied for $n \le 16$.

The GMD achieves **98% ARE** (Nair, 1936) with a **29.3% breakdown
point**, making it the most statistically efficient robust alternative
in this package. It is the estimator `scale_robust()` auto-switches to
for large samples.

``` r
gmd(c(1, 2, 3, 5, 7, 8))
gmd(c(1, 2, 3, 5, 7, 8), constant = 1)   # raw (unscaled)
```

### `adm(x, center, constant = 1.2533141373155001, na.rm = FALSE)`

Computes the mean absolute deviation from the median, scaled by a
consistency constant for asymptotic normality under the Gaussian model:

$$\text{ADM}(x) = C \cdot \frac{1}{n}\sum_{i=1}^{n} |x_i - \text{med}(x)|$$

where $C = \sqrt{\pi/2} \approx 1.2533$ (Nair, 1947). When `center` is
supplied, it replaces the median. The ADM achieves **88.3% ARE** but is
not itself robust against outliers (explosion breakdown $1/n$). It
serves as the fallback scale estimator when the MAD collapses to zero.

``` r
adm(c(1, 2, 3, 5, 7, 8))
adm(c(1, 2, 3, 5, 7, 8), constant = 1)   # without consistency correction
```

### `robLoc(x, scale = NULL, na.rm = FALSE, maxit = 80L, tol = sqrt(.Machine$double.eps))`

M-estimator for location using Newton–Raphson iteration with the
logistic psi function (Rousseeuw & Verboven 2002, Eq. 21). Starting
value: $T^{(0)} = \text{median}(x)$. Auxiliary scale:
$S = \text{MAD}(x)$ (or the user-supplied `scale`). See [Methodological
enhancements](#methodological-enhancements) for the iteration formula.

**Fallback logic:** When `scale` is unknown and $n < 4$, or when `scale`
is known and $n < 3$, the function returns `median(x)` without
iteration. Providing a known `scale` lowers the minimum sample size from
4 to 3 because the MAD (which is unreliable at $n = 3$) is no longer
needed. The flowchart below shows the full control flow including
fallbacks and the Newton–Raphson loop.

``` r
robLoc(c(1, 2, 3, 5, 7, 8))
robLoc(c(1, 2, 3), scale = 1.5)   # known scale enables n = 3
```

``` mermaid
flowchart TD
    A[Set minobs: 3 if scale known, 4 if unknown] --> B[med = median x]
    B --> C{n < minobs?}
    C -- Yes --> D([Return med])
    C -- No --> E{Scale known?}
    E -- Yes --> F[s = scale]
    E -- No --> G[s = MAD x]
    F --> H{s = 0?}
    G --> H
    H -- Yes --> I([Return med])
    H -- No --> J[t = med]
    J --> K["psi_i = tanh( (x_i - t) / 2s )"]
    K --> L["sum_psi = Sum psi_i, sum_dpsi = Sum (1 - psi_i^2)"]
    L --> M["t += 2s * sum_psi / sum_dpsi"]
    M --> N{"abs(v) <= tol?"}
    N -- No --> K
    N -- Yes --> O([Return t])
```

### `robScale(x, loc = NULL, fallback = c("adm", "na"), implbound = 1e-4, na.rm = FALSE, maxit = 80L, tol = sqrt(.Machine$double.eps))`

M-estimator for scale using multiplicative iteration with the rho
function—the square of the logistic psi (Rousseeuw & Verboven 2002, Eq.
27):

$$S^{(k)} = S^{(k-1)} \cdot \sqrt{2 \cdot \frac{1}{n}\sum \psi_{\log}^2\!\left(\frac{x_i - T}{c \cdot S^{(k-1)}}\right)}$$

where $c = 0.37394112142347236$ and $T = \text{median}(x)$ is held
fixed. Starting value: $S^{(0)} = \text{MAD}(x)$.

**Degenerate input handling:** When the sample size falls below the
minimum for iteration (4 for unknown location, 3 for known), the
function returns the initial MAD-based scale directly if it is nonzero.
When the MAD collapses to zero (i.e. $\text{MAD} \leq$ `implbound`), the
`fallback` argument controls the result:

- `fallback = "adm"` (Default): returns `adm(x)`, maintaining a finite
  robust estimate where standard scale measures fail.
- `fallback = "na"`: returns `NA`, strictly matching the behavioral
  profile of the `revss` package.

Providing a known `loc` centers the data at that value and uses the
median-distance-to-zero ($1.4826 \cdot \text{median}(|x_i - \mu|)$) as
the initial scale, lowering the minimum sample size from 4 to 3. The
flowchart below illustrates the control flow, including the `fallback`
logic and SIMD-accelerated loop.

``` r
robScale(c(1, 2, 3, 5, 7, 8))
robScale(c(5, 5, 5, 5, 6), fallback = "na")   # returns NA (revss compatibility)
```

``` mermaid
flowchart TD
    A{Location known?}
    A -- Yes --> B1[w = x - loc]
    B1 --> B2["s = 1.4826 * median( abs(w) )"]
    B2 --> B3[t = 0, minobs = 3]
    A -- No --> C1[med = median x]
    C1 --> C2[s = MAD x, t = med]
    C2 --> C3[minobs = 4]
    B3 --> D{n < minobs?}
    C3 --> D
    D -- Yes --> E{"s <= implbound?"}
    E -- Yes --> F{fallback?}
    F -- "adm" --> F1([Return ADM x])
    F -- "na" --> F2([Return NA])
    E -- No --> G([Return s])
    D -- No --> H{s <= implbound?}
    H -- Yes --> F
    H -- No --> J["psi_i = tanh( (x_i - t) / (2cs) )"]
    J --> K["sum_rho = Sum psi_i^2"]
    K --> L["v = sqrt( 2 * sum_rho / n ), s *= v"]
    L --> M{"abs(v - 1) <= tol?"}
    M -- No --> J
    M -- Yes --> N([Return s])
```

### `qn(x, constant = 2.2191, finite.corr = TRUE, na.rm = FALSE)`

Computes the $Q_n$ estimator of scale (Rousseeuw & Croux, 1993). Unlike
M-estimators, $Q_n$ requires no location estimate and achieves a 50%
breakdown point. `robscale` implements $Q_n$ with a tiered strategy: a
brute-force exact algorithm for small $n$ (below `qn_exact_threshold`)
and a cache-aware parallelized Johnson-style algorithm for larger
samples.

``` r
qn(c(1, 2, 3, 5, 7, 8))
```

### `sn(x, constant = 1.1926, finite.corr = TRUE, na.rm = FALSE)`

Computes the $S_n$ estimator of scale (Rousseeuw & Croux, 1993). $S_n$
is more statistically efficient than the MAD and maintains a 50%
breakdown point. `robscale` uses optimal sorting networks for $n \le 16$
and a highly optimized parallelized inner-median algorithm for general
samples.

``` r
sn(c(1, 2, 3, 5, 7, 8))
```

### `iqr_scaled(x, constant = 0.741301109252801, na.rm = FALSE)`

Computes the interquartile range scaled by a consistency constant for
asymptotic normality under the Gaussian model:

$$\text{IQR}_s(x) = C \cdot (Q_{0.75} - Q_{0.25})$$

where $Q_p$ denotes the Type 7 quantile (R default) and
$C = 1/(\Phi^{-1}(0.75) - \Phi^{-1}(0.25)) \approx 0.7413$ (Bickel &
Lehmann, 1976). Unlike `stats::IQR()`, which requires a full
$O(n \log n)$ sort, this implementation uses dual $O(n)$ Floyd–Rivest
selection—one call per quartile—providing a substantial speedup for
large datasets. The IQR achieves **37% ARE** with a **25% breakdown
point**.

``` r
iqr_scaled(c(1, 2, 3, 5, 7, 8))
iqr_scaled(c(1, 2, 3, 5, 7, 8), constant = 1)   # raw IQR
```

### `mad_scaled(x, center, constant = 1.482602218505602, na.rm = FALSE)`

Computes the median absolute deviation from the median, scaled by a
consistency constant for asymptotic normality:

$$\text{MAD}_s(x) = C \cdot \text{med}_i\, |x_i - \text{med}(x)|$$

where $C = 1/\Phi^{-1}(3/4) \approx 1.4826$. Unlike `stats::mad()`, this
implementation uses $O(n)$ Floyd–Rivest selection with sorting networks
for $n \le 16$, avoiding a full sort. The MAD achieves **36.8% ARE**
(Rousseeuw & Croux, 1993: “about 37%”) with a **50% breakdown point**.

``` r
mad_scaled(c(1, 2, 3, 5, 7, 8))
mad_scaled(c(1, 2, 3, 5, 7, 8), constant = 1)   # raw MAD
```

### `scale_robust(x, method, auto_switch, threshold, n_boot, na.rm)`

Unified dispatcher for robust scale estimation. Operates in three modes:

1.  **Ensemble** (default, `method = "ensemble"`): variance-weighted
    combination of all 7 scale estimators via bootstrap resampling.
2.  **Auto-switch** (`auto_switch = TRUE`, default): when $n \ge$
    `threshold` (default 20), returns `gmd(x)` directly—the most
    efficient robust estimator at negligible cost.
3.  **Explicit method**: dispatches to a specific estimator by name.

``` r
scale_robust(c(1, 2, 3, 5, 7, 8))           # ensemble (n < 20)
scale_robust(rnorm(50))                       # auto-switches to gmd
scale_robust(rnorm(50), auto_switch = FALSE)  # forces ensemble
scale_robust(rnorm(50), method = "qn")        # explicit Qn
```

``` mermaid
flowchart TD
    A["scale_robust(x, method, auto_switch, threshold, n_boot)"] --> B{n < 2?}
    B -- Yes --> C([Return NA])
    B -- No --> D{auto_switch AND n >= threshold?}
    D -- Yes --> E(["Return gmd(x)"])
    D -- No --> F{method?}
    F -- ensemble --> G["Bootstrap n_boot resamples"]
    G --> H["Compute all 7 estimators per resample"]
    H --> I["Inverse-variance weights from bootstrap variances"]
    I --> J["Weighted sum of estimators on original data"]
    J --> K([Return weighted estimate])
    F -- "gmd/sd/mad/iqr/sn/qn/robScale" --> L(["Return method(x)"])
```

The ensemble combines: `sd_c4`, `gmd`, `mad_scaled`, `iqr_scaled`, `sn`,
`qn`, and `robScale`. Bootstrap resampling uses a deterministic
XorShift32 PRNG seeded by replicate index, ensuring reproducible results
without requiring `set.seed()`.

### Estimator comparison: the robustness–efficiency frontier

The 8 scale estimators in `robscale` occupy different points on the
trade-off between statistical efficiency and robustness to outliers:

<div id="tbl-frontier">

Table 3

| Estimator | ARE (%) | Breakdown (%) | Complexity | Best for |
|:---|:---|:---|:---|:---|
| `sd_c4` | 100 | 0 | $O(n)$ | Efficiency anchor; clean data only |
| `gmd` | 98 | 29.3 | $O(n \log n)$ | General-purpose; large samples |
| `adm` | 88.3 | $1/n$ | $O(n)$ | Fallback when MAD = 0 |
| `qn` | 82.3 | 50 | $O(n \log n)$ | High efficiency + high breakdown |
| `sn` | 58.2 | 50 | $O(n \log n)$ | High breakdown, moderate efficiency |
| `robScale` | 55.0 | 50 | $O(n)$ iters | High breakdown, small samples |
| `iqr_scaled` | 37 | 25 | $O(n)$ | Quick robust estimate; quartile-based |
| `mad_scaled` | 36.8 | 50 | $O(n)$ | Maximum breakdown; classical robust |

</div>

At the top of the frontier, `sd_c4` retains full efficiency but
collapses under a single outlier. At the bottom, `mad_scaled` tolerates
up to 50% contamination but wastes over 60% of the available
information. The `gmd` occupies the sweet spot: 98% efficiency with
29.3% breakdown—sufficient for most practical contamination levels. For
adversarial settings where breakdown must be maximized, `qn` provides
the best combination of high breakdown (50%) and high efficiency
(82.3%).

## Methodological enhancements

`robscale` implements the estimators defined by Rousseeuw & Verboven
(2002) and Rousseeuw & Croux (1993). It produces identical numerical
results to `revss` and `robustbase`, but its computational strategies
significantly cut run times across all sample sizes.

### 1. The tanh identity for the logistic psi function

The logistic psi function central to both M-estimators is:

$$\psi_{\log}(x) = \frac{e^x - 1}{e^x + 1}$$

`revss` evaluates this as `2 * plogis(u) - 1`, which calls R’s `plogis`
(the logistic CDF, $1/(1 + e^{-x})$). The computation requires one call
to `exp()` followed by two arithmetic operations, plus the overhead of
R’s vectorised dispatch, intermediate vector allocation, and
garbage-collection pressure.

The algebraic identity

$$\psi_{\log}(x) = \tanh(x/2)$$

is immediate from the definition of the hyperbolic tangent. `robscale`
exploits this identity to reduce $\psi_{\log}$ to a single `tanh` call.
The identity is not merely a cosmetic rewrite:

- **Branch elimination.** A direct implementation of
  $(e^x - 1)/(e^x + 1)$ overflows for large $|x|$, requiring a
  sign-based branch to keep intermediate values bounded. The `tanh`
  function handles this internally with a single code path.

- **Platform vectorization.** `robscale` uses platform-specific
  libraries to evaluate `tanh` in bulk. The implementation ranks and
  selects the fastest available backend:

  1.  **Apple Accelerate.** On macOS (Darwin), it uses `vvtanh` for
      array-wide SIMD processing.
  2.  **SLEEF.** On Linux (x86_64), it uses the SLEEF library to target
      AVX2 or AVX512 instruction sets. On macOS, this backend is
      explicitly disabled in favor of Apple Accelerate.
  3.  **OpenMP SIMD.** A compiler-guided fallback via
      `#pragma omp simd`.
  4.  **Scalar.** Standard `std::tanh` fallback.

### 2. Newton–Raphson iteration for location

`revss` iterates the location estimator using the scoring fixed-point
iteration (Rousseeuw & Verboven 2002, Eq. 21):

$$T^{(k+1)} = T^{(k)} + S \cdot \frac{\frac{1}{n}\sum \psi_{\log}\!\left(\frac{x_i - T^{(k)}}{S}\right)}{\alpha}$$

where $\alpha = \int \psi_{\log}'(u)\,d\Phi(u) \approx 0.4132$ is the
normalization constant. This is a fixed-point iteration with *linear*
convergence: each step reduces the error by a constant factor.

`robscale` instead applies Newton–Raphson to the estimating equation
$\sum \psi_{\log}((x_i - T)/S) = 0$, yielding:

$$T^{(k+1)} = T^{(k)} + \frac{2\,S\sum \psi\!\left(\frac{x_i - T^{(k)}}{2S}\right)}
{\sum \left[1 - \psi^2\!\left(\frac{x_i - T^{(k)}}{2S}\right)\right]}$$

where $\psi(\cdot) = \tanh(\cdot)$. The efficiency follows from
observing that the derivative of the logistic psi satisfies
$\psi_{\log}'(x) = 1 - \psi_{\log}^2(x)
= 1 - \tanh^2(x/2)$. Since $\tanh$ values have already been computed for
the numerator, the denominator requires only squaring and subtraction—no
additional transcendental function calls.

Newton–Raphson achieves *quadratic* convergence near the solution: the
number of correct digits approximately doubles per iteration. The
practical effect is a reduction from 4–8 iterations (scoring) to 3
iterations (Newton–Raphson) for reaching the same tolerance of
$\sqrt{\epsilon_{\text{mach}}} \approx 1.49
\times 10^{-8}$:

| $n$ | Scoring iterations | Newton–Raphson iterations |
|----:|-------------------:|--------------------------:|
|   4 |                  7 |                         3 |
|   5 |                  8 |                         3 |
|   8 |                  7 |                         3 |
|  20 |                  6 |                         3 |
| 100 |                  5 |                         3 |

At small $n$, the scoring iteration count is higher because the starting
value (the median) can be far from the M-estimate in units of the
auxiliary scale. Newton–Raphson absorbs this gap in fewer steps.

### 3. $O(n)$ median selection

Both the median and the MAD require computing a quantile—the median of
the data, and the median of the absolute deviations. `revss` uses R’s
`median()` and `mad()`, which call `sort()` internally: an $O(n \log n)$
operation.

`robscale` uses a tiered $O(n)$ median selection strategy. For even $n$,
a single linear scan over the upper partition locates the $(k{+}1)$th
element needed for averaging.

For $n \leq 16$—the core target regime—the selection step uses optimal
sorting networks (Knuth, TAOCP Vol. 3, Sec. 5.3.4; Dobbelaere’s verified
optimal networks for $n = 9$–$16$). These are conditional
compare-and-swap sequences—typically compiled to branchless machine code
at `-O2`—with the minimum number of comparisons for each $n$:

| $n$ | Comparators |
|----:|------------:|
|   3 |           3 |
|   4 |           5 |
|   5 |           9 |
|   6 |          12 |
|   7 |          16 |
|   8 |          19 |
|   9 |          25 |
|  10 |          29 |
|  12 |          39 |
|  14 |          51 |
|  16 |          60 |

Cross-platform benchmarking confirmed 2–4$\times$ speedups over
`std::sort` for $n = 9$–$16$ on both ARM64 (Apple Silicon) and x86_64
(AMD Zen 3).

For $17 \leq n < 600$, the code delegates to `std::nth_element`
(introselect), which provides $O(n)$ worst-case selection with
median-of-three pivot selection. For $n \geq 600$, the Floyd–Rivest
algorithm (Floyd & Rivest, 1975) applies a statistical narrowing step
that reduces the active window to $O(n^{2/3})$ elements before
partitioning—a constant-factor improvement that amortizes the overhead
of its `log`/`exp`/`sqrt` computation only at scale.

### 4. Arena allocation on the stack

Each estimator requires working arrays: a copy of the input (for
destructive selection), absolute deviations (for MAD), and a temporary
buffer (for bulk `tanh` arguments). `revss` allocates these as R
vectors, incurring R’s SEXPREC header overhead and adding
garbage-collection pressure.

`robscale` uses a tiered stack-allocated arena: a 128-double
micro-buffer for $n \leq 64$ and a 2,048-double buffer per array for
$n \leq 2{,}048$—which covers the target regime and far beyond—with zero
heap allocation. For $n > 2{,}048$, the code falls back to
`new[]`/`delete[]`.

### 5. Compile-time reciprocal constants

The constants $1/\alpha = 1/0.413241928283814$ and
$1/c = 1/0.37394112142347236$ are declared `constexpr`, allowing the
compiler to replace divisions in the iteration loop with
multiplications. On ARM64, this avoids the ~10-cycle `fdiv` instruction
in favour of a ~3-cycle `fmul`.

### 6. Loop-invariant hoisting

Values that are constant across iterations—`inv_s = 1.0 / s`,
`half_inv_s = 0.5 / s`, `inv_n = 1.0 / n`—are computed once before the
loop. The `revss` implementation recomputes `(x - t) / s` as a fresh R
vector each iteration, traversing the interpreter for every vectorised
operation.

### 7. $Q_n$ and $S_n$ algorithm optimizations

`robustbase` implements $Q_n$ and $S_n$ in R-wrapped Fortran (Maechler
et al.), which incurs R dispatch overhead for every inner-loop call.
`robscale` replaces these paths entirely with a self-contained C++17
implementation.

**$Q_n$ — tiered exact/approximate algorithm.** For small $n$ (below a
compile-time threshold), `robscale` enumerates all $\binom{n}{2}$
pairwise absolute differences, selects the $h$th order statistic with
Floyd–Rivest, and applies the finite-sample correction factor. For
larger $n$, it uses a parallel Johnson-style algorithm:

1.  **Sort** the data in $O(n \log n)$.
2.  **Count and bound** — two parallel workers (`QnCountWorker`,
    `QnRefineWorker`) scan the sorted array to count and bracket the
    number of pairs above/below a trial value, using TBB
    `parallel_reduce` and `parallel_for`.
3.  **Weighted median of inner values** — a weighted-median step
    (`whimed_cpp`) refines the bracket; iterations continue until the
    bracket width falls below $n$.
4.  **Final brute-force** — the residual window is enumerated and
    Floyd–Rivest selects the exact $k$th difference.

This structure avoids materializing all $O(n^2)$ pairs and runs at
$O(n \log n)$ per iteration, with parallelism scaling across all
available cores for $n \geq$ `qn_parallel_threshold`.

**$S_n$ — parallelized inner-median sweep.** The $S_n$ statistic is the
low median of the vector $\{\text{med}_j |x_i - x_j|\}_{i=1}^n$.
`robscale` computes each inner median with an initial binary search
seeding a sliding-window linear scan (exploiting sortedness) in
amortized $O(1)$ per element, then dispatches the outer $n$ iterations
across TBB threads via `SnWorker`. For $n \leq 2048$, a stack-allocated
arena avoids heap allocation entirely.

### 8. Welford’s algorithm for numerical stability

The `sd_c4()` estimator uses Welford’s (1962) one-pass online algorithm
for computing variance. The standard two-pass formula
$s^2 = \frac{1}{n-1}\sum(x_i -
\bar{x})^2$ suffers from catastrophic cancellation when $\sum x_i^2$ and
$n\bar{x}^2$ are close. Welford’s algorithm incrementally updates mean
and sum-of-squares-of-differences, maintaining full precision with a
single pass:

$$\delta_i = x_i - \bar{x}_{i-1}, \quad \bar{x}_i = \bar{x}_{i-1} + \delta_i / i, \quad M_{2,i} = M_{2,i-1} + \delta_i(x_i - \bar{x}_i)$$

### 9. Dual Floyd–Rivest selection for IQR

The `iqr_scaled()` estimator requires two quantiles ($Q_{0.25}$ and
$Q_{0.75}$). Rather than sorting the entire array, `robscale` performs
two independent $O(n)$ Floyd–Rivest selections—one per quartile. Each
selection call uses the Type 7 quantile interpolation (the R default),
locating the floor index via `floyd_rivest_select` and then finding the
next-smallest element via a linear scan for the fractional
interpolation. This gives $O(n)$ total complexity vs the $O(n \log n)$
full sort used by `stats::IQR`.

### 10. Variance-weighted ensemble bootstrap

The `scale_robust()` ensemble operates in C++ (`cpp_scale_ensemble`) to
avoid R-level overhead for the $n_{\text{boot}} \times 7$ estimator
evaluations. For each bootstrap replicate
$r = 1, \ldots, n_{\text{boot}}$:

1.  A deterministic XorShift32 PRNG (Marsaglia, 2003) seeded with
    $r + 12345$ draws $n$ indices with replacement.
2.  All 7 estimators are evaluated on the resampled data, sharing
    pre-allocated work buffers.

After bootstrapping, the inverse-variance weight for each estimator $j$
is:

$$w_j = \frac{1/\hat\sigma_j^2}{\sum_{k=1}^{7} 1/\hat\sigma_k^2}$$

where $\hat\sigma_j^2$ is the sample variance of estimator $j$ across
bootstrap replicates. The final estimate is
$\hat\sigma = \sum_j w_j \cdot \hat\sigma_j(x)$ evaluated on the
original data.

<div id="tbl-qn-bench">

Table 4

|      $n$ | `robustbase::Qn` | `robscale::qn` | Speedup  |
|---------:|:-----------------|:---------------|:---------|
|        8 | 9.2 µs           | 1.6 µs         | **5.8x** |
|       16 | 10.5 µs          | 1.7 µs         | **6.0x** |
|       64 | 14.1 µs          | 7.2 µs         | **2.0x** |
|     1024 | 463.9 µs         | 216.2 µs       | **2.1x** |
|    65536 | 57740.3 µs       | 10676.7 µs     | **5.4x** |
| 10000000 | 9.9 s            | 1.8 s          | **5.4x** |

</div>

<div id="tbl-sn-bench">

Table 5

|      $n$ | `robustbase::Sn` | `robscale::sn` | Speedup  |
|---------:|:-----------------|:---------------|:---------|
|        8 | 4.0 µs           | 1.6 µs         | **2.6x** |
|       16 | 4.5 µs           | 1.6 µs         | **2.8x** |
|       64 | 5.3 µs           | 2.1 µs         | **2.5x** |
|     1024 | 35.2 µs          | 19.4 µs        | **1.8x** |
|    65536 | 6680.1 µs        | 1066.1 µs      | **6.3x** |
| 10000000 | 1.4 s            | 0.2 s          | **5.5x** |

</div>

## Architecture overview

`robscale` uses a tiered dispatch architecture to select the most
efficient algorithm based on sample size and hardware capabilities:

``` mermaid
graph TD
    SR["scale_robust() dispatcher"] --> ENS{Ensemble or single?}
    ENS -- "n < threshold" --> BOOT["Bootstrap ensemble kernel<br/>(7 estimators x n_boot resamples)"]
    ENS -- "n >= threshold" --> GMD_FAST["gmd() direct"]
    ENS -- "explicit method" --> SINGLE["Single estimator dispatch"]

    subgraph "Scale Estimators"
        SD["sd_c4"]
        GMD["gmd"]
        ADM["adm"]
        QN["qn"]
        SN["sn"]
        RS["robScale"]
        IQR["iqr_scaled"]
        MAD["mad_scaled"]
    end

    subgraph "Location Estimators"
        RL["robLoc"]
    end

    BOOT --> SD & GMD & MAD & IQR & SN & QN & RS
    SINGLE --> SD & GMD & ADM & QN & SN & RS & IQR & MAD
    GMD_FAST --> GMD

    subgraph "Algorithm Tiers"
        T1["n <= 16: Sorting networks"]
        T2["17 <= n < L2 threshold: Optimized scalar C++"]
        T3["n >= L2 threshold: Parallel TBB kernels"]
    end

    QN & SN --> T1 & T2 & T3

    subgraph "Hardware Acceleration"
        G[AVX2 / AVX512 / NEON]
        H[Apple Accelerate]
        I[SLEEF Library]
    end

    T1 & T2 & T3 --> G & H & I
```

## Benchmarks

`robscale` targets three distinct performance regimes: very small
samples (the M-estimators `robLoc`, `robScale`, and `adm`, compared
against `revss`), general-size samples (the scale estimators `qn` and
`sn`, compared against `robustbase`), and the new scale estimators
(`gmd`, `iqr_scaled`, `mad_scaled`, compared against existing R
implementations). Figure 1 summarizes all three comparisons on Linux
(Ryzen 9 5900HX, CRAN-compatible build with auto-detected SIMD).

<div id="fig-benchmarks">

![](benchmarks/speedup_fig.png)

Figure 1: Median speedup factor (x) vs. sample size $n$. Panel A
compares `robLoc`, `robScale`, and `adm` against `revss`; Panel B
compares `qn` and `sn` against `robustbase`; Panel C compares `gmd`,
`iqr_scaled`, and `mad_scaled` against existing R implementations.

</div>

**Benchmark environment:**

- **Machine:** AMD Ryzen 9 5900HX with Radeon Graphics
- **CPU Governor:** powersave
- **OS:** Arch Linux
- **R version:** R version 4.5.3 (2026-03-11)
- **Package version:** 0.2.0
- **SLEEF (Optimized Build):** Detected
- **Date:** 2026-03-14

### Small-sample M-estimators vs. `revss` (Panel A)

In the target regime ($n \le 20$), `robscale` outperforms `revss` by
**10.6–31.0x**. Drivers include:

- Transitioning from interpreted R to compiled C++17.
- Achieving quadratic convergence with Newton–Raphson (3 iterations vs
  6–8 for scoring).
- Eliminating heap allocation via stack-allocated memory arenas.
- Deploying optimal sorting networks for $n \le 16$.

Even at $n = 16{,}384$, the gains remain **3.4–5.9x** because the
interpreter overhead of `revss` scales with the number of Newton–Raphson
iterations, not just vector length.

### Scale estimators vs. `robustbase` (Panel B)

For $Q_n$ and $S_n$, the performance story follows two regimes separated
by the parallelism threshold:

**Small to medium samples ($n \le 10^3$).** `robscale` leads by
**1.7–6.8x**. The gain comes primarily from the avoidance of R dispatch
overhead and the use of stack memory. For $Q_n$ at $n = 8$, the
brute-force exact algorithm completes in 1.6 µs vs. 9.2 µs for
`robustbase` — a **5.8x** edge.

**Large samples ($n \ge 10^4$).** The advantage grows to **1.9–6.9x** as
TBB parallelism engages. At $n = 10^7$, `qn` runs in 1.8 s vs. 9.9 s for
`robustbase::Qn` (**5.4x**), and `sn` runs in 0.2 s vs. 1.4 s
(**5.5x**). Parallel efficiency is bounded by Amdahl’s Law and memory
bandwidth; while the multi-threaded kernels provide substantial gains
for massive datasets, speedups do not scale linearly with thread count.

### New scale estimators vs. existing R implementations (Panel C)

<div id="tbl-new-bench">

Table 6

|      $n$ | Comparison          | Speedup   |
|---------:|:--------------------|:----------|
|       64 | gmd vs GiniDistance | **10.7x** |
|       64 | gmd vs Hmisc        | **16.3x** |
|       64 | iqr_scaled vs stats | **27.2x** |
|       64 | mad_scaled vs stats | **19.7x** |
|     1024 | gmd vs GiniDistance | **3.3x**  |
|     1024 | gmd vs Hmisc        | **4.7x**  |
|     1024 | iqr_scaled vs stats | **9.2x**  |
|     1024 | mad_scaled vs stats | **8.9x**  |
|    65536 | gmd vs GiniDistance | **2.9x**  |
|    65536 | gmd vs Hmisc        | **3.8x**  |
|    65536 | iqr_scaled vs stats | **1.4x**  |
|    65536 | mad_scaled vs stats | **2.0x**  |
| 10000000 | gmd vs GiniDistance | **2.7x**  |
| 10000000 | gmd vs Hmisc        | **3.7x**  |
| 10000000 | iqr_scaled vs stats | **1.0x**  |
| 10000000 | mad_scaled vs stats | **2.9x**  |

</div>

**GMD** (`robscale::gmd` vs `Hmisc::GiniMd`): **2.4–17.5x** speedup.
`Hmisc::GiniMd` is a pure R implementation using the same $O(n \log n)$
order-statistics formula. The speedup comes from C++ compilation and
sorting networks for small $n$. Against `GiniDistance::gmd` (an
Rcpp-backed C++ implementation), the comparison is **1.9–11.4x**—a
tighter race since both are compiled, with `robscale`’s advantage coming
from sorting networks and the consistency-constant integration.

**IQR** (`robscale::iqr_scaled` vs `stats::IQR`): **1.0–30.2x** speedup.
`stats::IQR` performs a full $O(n \log n)$ sort via `quantile()`.
`robscale` uses dual $O(n)$ Floyd–Rivest selection, which dominates at
large $n$.

**MAD** (`robscale::mad_scaled` vs `stats::mad`): **2.0–21.4x** speedup.
`stats::mad` also performs a full sort for the median step. `robscale`
uses $O(n)$ selection throughout.

> \[!NOTE\] **Source builds recommended.** Installing from source
> (`install.packages("robscale", type = "source")`) enables the
> `configure` script to detect SIMD capabilities (AVX2/FMA on x86_64,
> NEON on ARM64) and link platform-specific libraries (Apple Accelerate,
> SLEEF). Pre-built CRAN binaries use portable settings and may not
> include these optimizations. Parallelism thresholds are derived from
> the detected per-core L2 cache size at runtime on all platforms.

## Numerical equivalence

The test suite verifies `robscale` against reference implementations:

**Legacy M-estimators** (`tests/testthat/test-cross-check.R`): compares
`robscale` and `revss` across 5,400 randomly generated inputs
($n = 3, 4,
\ldots, 20$; 100 replicates per $n$; all three functions `adm`,
`robLoc`, `robScale`). All 5,400 comparisons pass at tolerance
$\sqrt{\epsilon_{\text{mach}}} \approx 1.49 \times 10^{-8}$.

**New estimators:**

- `gmd`: exact match with the R formula
  $C \cdot 2/(n(n-1)) \sum (2i - n - 1) x_{(i)}$ (`test-gmd.R`)
- `iqr_scaled`: matches `IQR(x) * 0.741301109252801` (`test-iqr.R`)
- `mad_scaled`: matches `stats::mad(x)` (`test-mad-scaled.R`)
- `sd_c4`: matches `sd(x) / c4(n)` (`test-sd-c4.R`)
- `scale_robust` ensemble: deterministic via fixed XorShift32 seeds
  (`test-ensemble.R`, `test-scale-robust.R`)

The Newton–Raphson iteration converges to the same fixed point as the
scoring iteration—it solves the same estimating equation—so results
differ only by rounding at the level of the convergence tolerance.

## Mathematical background

The estimators follow Rousseeuw & Verboven (2002), Rousseeuw & Croux
(1993), Gini (1912), and Nair (1936, 1947). A brief summary of the key
definitions:

**Logistic psi function** (Eq. 23):

$$\psi_{\log}(x) = \frac{e^x - 1}{e^x + 1} = \tanh(x/2)$$

Bounded in $(-1, 1)$, smooth ($C^\infty$), strictly monotone.
Boundedness provides robustness; smoothness avoids the corner artifacts
of Huber’s psi at small $n$.

**Decoupled estimation.** Location and scale are estimated separately
with a fixed auxiliary estimate, breaking the positive-feedback loop of
Huber’s Proposal 2. `robLoc` fixes scale at $\text{MAD}(x)$; `robScale`
fixes location at $\text{median}(x)$.

**Rho function for scale** (Eq. 26):

$$\rho_{\log}(x) = \psi_{\log}^2(x / c)$$

where $c = 0.37394112142347236$ is the constant that yields a 50%
breakdown point.

**$Q_n$ and $S_n$ statistics.**
$Q_n = c_n \cdot d \cdot \{|x_i - x_j|; i < j\}_{(k)}$ where
$k = \binom{h}{2}$, $h = \lfloor n/2 \rfloor + 1$, and $d = 2.2191$
(consistency constant for Gaussian data). $S_n = c_n' \cdot 1.1926 \cdot
\text{lomed}_i \{\text{himed}_j |x_i - x_j|\}$, where $\text{lomed}$ and
$\text{himed}$ denote the low and high medians respectively.

**$c_4(n)$ correction factor.** The expected value of the sample
standard deviation under normality is $\sigma \cdot c_4(n)$ where:

$$c_4(n) = \sqrt{\frac{2}{n{-}1}} \cdot \frac{\Gamma(n/2)}{\Gamma((n{-}1)/2)}$$

Dividing $s$ by $c_4(n)$ yields an unbiased estimator of $\sigma$.

**Gini mean difference.** The order-statistics form:

$$\text{GMD}(x) = C_{\text{GMD}} \cdot \frac{2}{n(n{-}1)}\sum_{i=1}^{n} (2i - n - 1)\, x_{(i)}$$

is algebraically equivalent to the pairwise-difference definition
$\frac{1}{\binom{n}{2}}\sum_{i<j}|x_i - x_j|$ but avoids materializing
$O(n^2)$ pairs.

**Ensemble weighting formula.** Given $J$ estimators with bootstrap
variances $\hat\sigma_j^2$, the inverse-variance weighted estimate is:

$$\hat\sigma = \frac{\sum_{j=1}^{J} \hat\sigma_j(x) / \hat\sigma_j^2}{\sum_{j=1}^{J} 1/\hat\sigma_j^2}$$

**Key constants** (full double precision):

| Symbol | Value | Definition |
|:---|:---|:---|
| $\alpha$ | `0.413241928283814` | $\int \psi_{\log}'(u)\,d\Phi(u)$; scoring normalization constant |
| $c$ | `0.37394112142347236` | Solution to $\int \rho_{\log}(u)\,d\Phi(u) = 0.5$; scale rho constant |
| $C_{\text{ADM}}$ | `1.2533141373155001` | $\sqrt{\pi/2}$; ADM consistency constant |
| $C_{\text{MAD}}$ | `1.482602218505602` | $1/\Phi^{-1}(3/4)$; MAD consistency constant |
| $C_{\text{GMD}}$ | `0.886226925452758` | $\sqrt{\pi}/2$; GMD consistency constant |
| $C_{\text{IQR}}$ | `0.741301109252801` | $1/(\Phi^{-1}(0.75) - \Phi^{-1}(0.25))$; IQR consistency constant |

## Relation to revss and robustbase

This package re-implements the M-estimators from the
[revss](https://CRAN.R-project.org/package=revss) package (Avraham
Adler) and the $Q_n$ and $S_n$ estimators from
[robustbase](https://CRAN.R-project.org/package=robustbase) (Maechler et
al.).

The API for the M-estimators is intentionally identical to `revss`:
`adm()`, `robLoc()`, and `robScale()` accept the same arguments and
return the same values. Code that uses `revss` can switch to `robscale`
by changing only the `library()` call.

For `qn()` and `sn()`, the function signatures match `robustbase::Qn()`
and `robustbase::Sn()` (with lowercase names for consistency).

The v0.2.0 additions—`gmd()`, `iqr_scaled()`, `sd_c4()`, `mad_scaled()`,
`scale_robust()`, and `get_consistency_constant()`—have no direct
counterpart in `revss` or `robustbase`. `mad_scaled()` is API-compatible
with `stats::mad()` when `constant = 1`; similarly,
`iqr_scaled(x, constant = 1)` matches `stats::IQR(x)`.

Users who do not need compiled performance—or who prefer a
dependency-free pure-R package—should use `revss` or `robustbase`
directly. Both are mature, well-tested, and widely available.

## References

Bickel, P.J. and Lehmann, E.L. (1976). Descriptive Statistics for
Nonparametric Models III. Dispersion. *Annals of Statistics*, **4**(6),
1139–1158.
[doi:10.1214/aos/1176343648](https://doi.org/10.1214/aos/1176343648)

Floyd, R.W. and Rivest, R.L. (1975). Expected time bounds for selection.
*Communications of the ACM*, **18**(3), 165–172.
[doi:10.1145/360680.360691](https://doi.org/10.1145/360680.360691)

Gini, C. (1912). *Variabilita e mutabilita*. Bologna: Tipografia di
Paolo Cuppini.

Marsaglia, G. (2003). Xorshift RNGs. *Journal of Statistical Software*,
**8**(14), 1–6.
[doi:10.18637/jss.v008.i14](https://doi.org/10.18637/jss.v008.i14)

Nair, K.R. (1936). On the Mean Deviation. *Biometrika*, **28**(3/4),
428–436. [doi:10.2307/2333958](https://doi.org/10.2307/2333958)

Nair, K.R. (1947). A Note on the Mean Deviation from the Median.
*Biometrika*, **34**(3/4), 360–362.
[doi:10.2307/2332448](https://doi.org/10.2307/2332448)

Rousseeuw, P.J. and Croux, C. (1993). Alternatives to the Median
Absolute Deviation. *Journal of the American Statistical Association*,
**88**, 1273–1283.
[doi:10.1080/01621459.1993.10476408](https://doi.org/10.1080/01621459.1993.10476408)

Rousseeuw, P.J. and Verboven, S. (2002). Robust estimation in very small
samples. *Computational Statistics & Data Analysis*, **40**(4), 741–758.
[doi:10.1016/S0167-9473(02)00078-6](https://doi.org/10.1016/S0167-9473(02)00078-6)

Welford, B.P. (1962). Note on a Method for Calculating Corrected Sums of
Squares and Products. *Technometrics*, **4**(3), 419–420.
[doi:10.1080/00401706.1962.10490022](https://doi.org/10.1080/00401706.1962.10490022)

## Author

Dennis Alexis Valin Dittrich
([ORCID](https://orcid.org/0000-0002-4438-8276))

## License

MIT. Copyright 2026 Dennis Alexis Valin Dittrich.
