# robscale: Faster Robustness for Small and Large Samples

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18828607.svg)](https://doi.org/10.5281/zenodo.18828607)

## Overview

Outliers compromise the reliability of classical estimators even in moderate
samples. Yet the most widely used robust implementations in R—`revss` for
small-sample M-estimation and `robustbase` for the $Q_n$ and $S_n$ scale
estimators—carry significant computational overhead that blocks real-time
genomics pipelines, high-frequency finance, and other time-critical workflows.

`robscale` resolves both bottlenecks. Its C++17 kernels replace interpreted R
loops with SIMD-accelerated transcendental functions, Newton–Raphson iteration,
and parallelized $O(n \log n)$ algorithms. Against `revss`, the package achieves
**11–39×** speedups for the small-sample M-estimators. Against `robustbase`, it
achieves **1.6–6.9×** for $S_n$ and **1.6–10×** for $Q_n$ —with gains peaking near
**10×** at $n = 10^7$ as TBB parallelism reduces the computational bottleneck for massive datasets.

## Installation

```r
install.packages("robscale")

# Development version:
# install.packages("remotes")
# remotes::install_github("davdittrich/robscale")
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

Outliers pull classical estimators away from the bulk of the data. The robust
estimators remain stable.

## API reference

| Function | Purpose | ARE | Breakdown | Speedup vs. `revss` | Speedup vs. `robustbase` |
| :--- | :--- | :--- | :--- | :--- | :--- |
| `qn(x)` | $Q_n$ scale estimator | **82.3%** | 50% | — | **1.6x – 10x** |
| `sn(x)` | $S_n$ scale estimator | **58.2%** | 50% | — | **1.8x – 6.9x** |
| `robLoc(x)` | M-estimate of location | **41.3%** | 50% | **15x – 27x** | — |
| `robScale(x)` | M-estimate of scale | **60.5%** | 50% | **32x – 39x** | — |
| `adm(x)` | Average distance to median | **88.3%** | $1/n$ | **11x – 13x** | — |

All functions accept `na.rm` (default `FALSE`).

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
user-supplied `scale`). See [Methodological enhancements](#methodological-enhancements)
for the iteration formula.

**Fallback logic:** When `scale` is unknown and $n < 4$, or when `scale` is
known and $n < 3$, the function returns `median(x)` without iteration.
Providing a known `scale` lowers the minimum sample size from 4 to 3 because the
MAD (which is unreliable at $n = 3$) is no longer needed. The flowchart below
shows the full control flow including fallbacks and the Newton--Raphson loop.

```r
robLoc(c(1, 2, 3, 5, 7, 8))
robLoc(c(1, 2, 3), scale = 1.5)   # known scale enables n = 3
```

```mermaid
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

M-estimator for scale using multiplicative iteration with the rho function—the
square of the logistic psi (Rousseeuw & Verboven 2002, Eq. 27):

$$S^{(k)} = S^{(k-1)} \cdot \sqrt{2 \cdot \frac{1}{n}\sum \psi_{\log}^2\!\left(\frac{x_i - T}{c \cdot S^{(k-1)}}\right)}$$

where $c = 0.37394112142347236$ and $T = \text{median}(x)$ is held fixed.
Starting value: $S^{(0)} = \text{MAD}(x)$.

**Degenerate input handling:** When the sample size falls below the minimum for
iteration (4 for unknown location, 3 for known) or the MAD collapses to zero,
the function applies the logic selected by the `fallback` argument:

- `fallback = "adm"` (Default): returns `adm(x)` if $\text{MAD} \leq$ `implbound`, maintaining a finite
  robust estimate where standard scale measures fail.
- `fallback = "na"`: returns `NA`, strictly matching the behavioral profile of the `revss` package.

Providing a known `loc` centers the data at that value and uses the
median-distance-to-zero ($1.4826 \cdot \text{median}(|x_i - \mu|)$) as the
initial scale, lowering the minimum sample size from 4 to 3. The flowchart
below illustrates the control flow, including the `fallback` logic and SIMD-accelerated loop.

```r
robScale(c(1, 2, 3, 5, 7, 8))
robScale(c(5, 5, 5, 5, 6), fallback = "na")   # returns NA (revss compatibility)
```

```mermaid
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
M-estimators, $Q_n$ requires no location estimate and achieves a 50% breakdown
point. `robscale` implements $Q_n$ with a tiered strategy: a brute-force exact
algorithm for small $n$ (below `qn_exact_threshold`) and a cache-aware
parallelized Johnson-style algorithm for larger samples.

```r
qn(c(1, 2, 3, 5, 7, 8))
```

### `sn(x, constant = 1.1926, finite.corr = TRUE, na.rm = FALSE)`

Computes the $S_n$ estimator of scale (Rousseeuw & Croux, 1993). $S_n$ is more
statistically efficient than the MAD and maintains a 50% breakdown point.
`robscale` uses optimal sorting networks for $n \le 8$ and a highly
optimized parallelized inner-median algorithm for general samples.

```r
sn(c(1, 2, 3, 5, 7, 8))
```

### Estimator Convergence and Consistency

The three scale estimators—`qn()`, `sn()`, and `robScale()`—all serve as robust
alternatives to the standard deviation. Because `robscale` implements their
canonical consistency constants ($c=2.2191$ for $Q_n$, $c=1.1926$ for $S_n$),
they converge asymptotically to the same parameter $\sigma$ under a Gaussian
distribution.

While numerically similar on Gaussian data, they occupy different points on the
robustness--efficiency frontier:

1. **Precision**: $Q_n$ is the most statistically efficient, retaining **82.3%**
   asymptotic relative efficiency (ARE) compared to the sample standard deviation.
2. **Robustness**: All three maintain a **50% breakdown point**, though $S_n$
   and `robScale` trade higher efficiency for lower computational complexity or
   better behavior in specific small-sample regimes.

Users should expect point estimates to align closely on clean data. On
contaminated data, $Q_n$ provides the most stable performance across the
widest range of outliers.

## Methodological enhancements

`robscale` implements the estimators defined by Rousseeuw & Verboven (2002) and
Rousseeuw & Croux (1993). It produces identical numerical results to `revss`
and `robustbase`, but its computational strategies significantly cut run times
across all sample sizes.

### 1. The tanh identity for the logistic psi function

The logistic psi function central to both M-estimators is:

$$\psi_{\log}(x) = \frac{e^x - 1}{e^x + 1}$$

`revss` evaluates this as `2 * plogis(u) - 1`, which calls R's `plogis`
(the logistic CDF, $1/(1 + e^{-x})$). The computation requires one call to
`exp()` followed by two arithmetic operations, plus the overhead of R's
vectorised dispatch, intermediate vector allocation, and garbage-collection
pressure.

The algebraic identity

$$\psi_{\log}(x) = \tanh(x/2)$$

is immediate from the definition of the hyperbolic tangent. `robscale` exploits
this identity to reduce $\psi_{\log}$ to a single `tanh` call. The identity is not merely a cosmetic rewrite:

- **Branch elimination.** A direct implementation of $(e^x - 1)/(e^x + 1)$
  overflows for large $|x|$, requiring a sign-based branch to keep intermediate
  values bounded. The `tanh` function handles this internally with a single code
  path.

- **Platform vectorization.** `robscale` uses platform-specific libraries to evaluate `tanh` in bulk.
  The implementation ranks and selects the fastest available backend:
  1. **Apple Accelerate.** On macOS (Darwin), it uses `vvtanh` for array-wide SIMD processing.
  2. **SLEEF.** On Linux (x86_64), it uses the SLEEF library to target AVX2 or AVX512 instruction sets. On macOS, this backend is explicitly disabled in favor of Apple Accelerate.
  3. **OpenMP SIMD.** A compiler-guided fallback via `#pragma omp simd`.
  4. **Scalar.** Standard `std::tanh` fallback.

### 2. Newton--Raphson iteration for location

`revss` iterates the location estimator using the scoring fixed-point iteration
(Rousseeuw & Verboven 2002, Eq. 21):

$$T^{(k+1)} = T^{(k)} + S \cdot \frac{\frac{1}{n}\sum \psi_{\log}\!\left(\frac{x_i - T^{(k)}}{S}\right)}{\alpha}$$

where $\alpha = \int \psi_{\log}'(u)\,d\Phi(u) \approx 0.4132$ is the
normalization constant. This is a fixed-point iteration with *linear*
convergence: each step reduces the error by a constant factor.

`robscale` instead applies Newton--Raphson to the estimating equation
$\sum \psi_{\log}((x_i - T)/S) = 0$, yielding:

$$T^{(k+1)} = T^{(k)} + \frac{2\,S\sum \psi\!\left(\frac{x_i - T^{(k)}}{2S}\right)}
{\sum \left[1 - \psi^2\!\left(\frac{x_i - T^{(k)}}{2S}\right)\right]}$$

where $\psi(\cdot) = \tanh(\cdot)$. The efficiency follows from observing that the
derivative of the logistic psi satisfies $\psi_{\log}'(x) = 1 - \psi_{\log}^2(x)
= 1 - \tanh^2(x/2)$.
 Since $\tanh$ values have already been computed for the
numerator, the denominator requires only squaring and subtraction---no additional
transcendental function calls.

Newton--Raphson achieves *quadratic* convergence near the solution: the number
of correct digits approximately doubles per iteration. The practical effect is
a reduction from 4--8 iterations (scoring) to 3 iterations (Newton--Raphson)
for reaching the same tolerance of $\sqrt{\epsilon_{\text{mach}}} \approx 1.49
\times 10^{-8}$:

| $n$ | Scoring iterations | Newton--Raphson iterations |
| ---: | ---: | ---: |
| 4 | 7 | 3 |
| 5 | 8 | 3 |
| 8 | 7 | 3 |
| 20 | 6 | 3 |
| 100 | 5 | 3 |

At small $n$, the scoring iteration count is higher because the starting value
(the median) can be far from the M-estimate in units of the auxiliary scale.
Newton--Raphson absorbs this gap in fewer steps.

### 3. $O(n)$ median selection

Both the median and the MAD require computing a quantile---the median of the
data, and the median of the absolute deviations. `revss` uses R's `median()`
and `mad()`, which call `sort()` internally: an $O(n \log n)$ operation.

`robscale` uses a tiered $O(n)$ median selection strategy. For even $n$, a
single linear scan over the upper partition locates the $(k{+}1)$th element
needed for averaging.

For $n \leq 8$---the core target regime---the selection step uses optimal
sorting networks (Knuth, TAOCP Vol. 3, Sec. 5.3.4). These are conditional
compare-and-swap sequences---typically compiled to branchless machine code at
`-O2`---with the minimum number of comparisons for each $n$:

| $n$ | Comparators |
| ---: | ---: |
| 3 | 3 |
| 4 | 5 |
| 5 | 9 |
| 6 | 12 |
| 7 | 16 |
| 8 | 19 |

For $9 \leq n < 600$, the code delegates to `std::nth_element` (introselect),
which provides $O(n)$ worst-case selection with median-of-three pivot
selection. For $n \geq 600$, the Floyd--Rivest algorithm (Floyd & Rivest,
1975) applies a statistical narrowing step that reduces the active window to
$O(n^{2/3})$ elements before partitioning---a constant-factor improvement that
amortizes the overhead of its `log`/`exp`/`sqrt` computation only at scale.

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

### 7. $Q_n$ and $S_n$ algorithm optimizations

`robustbase` implements $Q_n$ and $S_n$ in R-wrapped Fortran (Maechler et al.),
which incurs R dispatch overhead for every inner-loop call. `robscale` replaces
these paths entirely with a self-contained C++17 implementation.

**$Q_n$ — tiered exact/approximate algorithm.** For small $n$ (below a
compile-time threshold), `robscale` enumerates all $\binom{n}{2}$ pairwise
absolute differences, selects the $h$th order statistic with Floyd--Rivest,
and applies the finite-sample correction factor. For larger $n$, it uses a parallel
Johnson-style algorithm:

1. **Sort** the data in $O(n \log n)$.
2. **Count and bound** — two parallel workers (`QnCountWorker`, `QnRefineWorker`)
   scan the sorted array to count and bracket the number of pairs above/below a
   trial value, using TBB `parallel_reduce` and `parallel_for`.
3. **Weighted median of inner values** — a weighted-median step (`whimed_cpp`)
   refines the bracket; iterations continue until the bracket width falls below
   $n$.
4. **Final brute-force** — the residual window is enumerated and Floyd--Rivest
   selects the exact $k$th difference.

This structure avoids materializing all $O(n^2)$ pairs and runs at $O(n \log n)$
per iteration, with parallelism scaling across all available cores for
$n \geq$ `qn_parallel_threshold`.

**$S_n$ — parallelized inner-median sweep.** The $S_n$ statistic is the low
median of the vector $\{\text{med}_j |x_i - x_j|\}_{i=1}^n$. `robscale`
computes each inner median with a sliding-window binary search (exploiting
sortedness) in O(log n) per element, then dispatches the outer $n$ iterations
across TBB threads via `SnWorker`. For $n \leq 2048$, a
stack-allocated arena avoids heap allocation entirely.

| $n$ | `robustbase::Qn` ($\mu s$) | `robscale::qn` ($\mu s$) | Speedup |
| ---: | ---: | ---: | ---: |
| 8 | 8.9 | 1.3 | **6.8×** |
| 64 | 13.6 | 5.3 | **2.6×** |
| 1 024 | 462.7 | 276.9 | **1.7×** |
| 65 536 | 49 324 | 12 247 | **4.0×** |
| 10 000 000 | 23 374 593 | 2 403 552 | **9.7×** |

| $n$ | `robustbase::Sn` ($\mu s$) | `robscale::sn` ($\mu s$) | Speedup |
| ---: | ---: | ---: | ---: |
| 8 | 4.1 | 1.3 | **3.3×** |
| 64 | 5.3 | 1.8 | **3.0×** |
| 1 024 | 33.0 | 17.4 | **1.9×** |
| 65 536 | 7 265 | 1 315 | **5.5×** |
| 10 000 000 | 1 461 045 | 250 588 | **5.8×** |

## Architecture overview

`robscale` uses a tiered dispatch architecture to select the most efficient algorithm
based on sample size and hardware capabilities:

```mermaid
graph TD
    A[R API: qn, sn, robLoc, robScale] --> B{Dispatch Layer}
    B -- "n <= 8" --> C[Sorting Networks]
    B -- "8 < n < 32768" --> D[Highly optimized scalar C++]
    B -- "n >= 32768" --> E[Parallel multi-threaded kernels]
    subgraph "Hardware Acceleration"
        G[AVX2 / AVX512 / NEON]
        H[Apple Accelerate]
        I[SLEEF Library]
    end
    E & D & C --> G & H & I
```

## Benchmarks

`robscale` targets two distinct performance regimes: very small samples (the
M-estimators `robLoc`, `robScale`, and `adm`, compared against `revss`) and
general-size samples (the scale estimators `qn` and `sn`, compared against
`robustbase`). Figure 1 summarizes both comparisons on Linux (Ryzen 9 5900HX,
CRAN-compatible build with no ROBSCALE_FAST optimizations).

![Two-panel benchmark: M-estimators vs. revss (left) and scale estimators vs. robustbase (right)](benchmarks/cran_mode/full_speedup_panel.png)
*Figure 1: Median speedup factor (×) vs. sample size $n$ on Linux (Ryzen 9 5900HX).
Panel A compares `robLoc`, `robScale`, and `adm` against `revss`; Panel B compares
`qn` and `sn` against `robustbase`. Shaded bands: 95% CI. Dashed line: parity (1×).*

### Small-sample M-estimators vs. `revss` (Panel A)

In the target regime ($n \le 20$), `robscale` outperforms `revss` by
**11×–39×**. Drivers include:

- Transitioning from interpreted R to compiled C++17.
- Achieving quadratic convergence with Newton–Raphson (3 iterations vs 6–8 for scoring).
- Eliminating heap allocation via stack-allocated memory arenas.
- Deploying optimal sorting networks for $n \le 8$.

Even at $n = 16{,}384$, the gains remain **3.5×–4.3×** because the interpreter
overhead of `revss` scales with the number of Newton–Raphson iterations, not
just vector length.

### Scale estimators vs. `robustbase` (Panel B)

For $Q_n$ and $S_n$, the performance story follows two regimes separated by
the parallelism threshold:

**Small to medium samples ($n \le 10^3$).** `robscale` leads by **1.7×–6.8×**.
The gain comes primarily from the avoidance of R dispatch overhead and the
use of stack memory. For $Q_n$ at $n = 8$, the brute-force exact algorithm
completes in 1.3 µs vs. 8.9 µs for `robustbase` — a **6.8×** edge.

**Large samples ($n \ge 10^4$).** The advantage grows to **4×–10×** as TBB
parallelism engages. At $n = 10^7$, `qn` runs in 2.4 s vs. 23.4 s for
`robustbase::Qn` (**9.7×**), and `sn` runs in 0.25 s vs. 1.46 s (**5.8×**).
Parallel efficiency is bounded by Amdahl's Law and memory bandwidth; while the
multi-threaded kernels provide substantial gains for massive datasets,
speedups do not scale linearly with thread count.

> [!IMPORTANT]
> **Peak Performance Note**: To unlock maximum hardware-specific speedups
> (including SIMD auto-discovery and microbenchmark-tuned thresholds), install
> the package from source with the `ROBSCALE_FAST=1` environment variable set:
>
> **From the terminal:**
> ```bash
> ROBSCALE_FAST=1 R CMD INSTALL robscale
> ```
>
> **From an R session (RStudio):**
> ```r
> Sys.setenv(ROBSCALE_FAST = "1")
> # Install from CRAN (compiling from source for peak performance)
> install.packages("robscale", type = "source")
>
> # Or install from GitHub (requires remotes/devtools)
> remotes::install_github("davdittrich/robscale")
> ```
>
> Default CRAN binaries use conservative, portable settings to ensure
> compatibility across all architectures.

## Numerical equivalence

The test suite (`tests/testthat/test-cross-check.R`) compares `robscale` and
`revss` across 5,400 randomly generated inputs:

- Sample sizes $n = 3, 4, \ldots, 20$
- 100 replicates per $n$ (uniform on $[-100, 100]$, seed 42)
- All three functions: `adm`, `robLoc`, `robScale`
- Tolerance: $\sqrt{\epsilon_{\text{mach}}} \approx 1.49 \times 10^{-8}$

All 5,400 comparisons pass. The Newton--Raphson iteration converges to the same
fixed point as the scoring iteration---it solves the same estimating equation---so the results
differ only by rounding at the level of the convergence tolerance.

## Mathematical background

The estimators follow Rousseeuw & Verboven (2002) and Rousseeuw & Croux (1993).
A brief summary of the key definitions:

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

where $c = 0.37394112142347236$ is the constant that yields a 50% breakdown point.

**$Q_n$ and $S_n$ statistics.** $Q_n = c_n \cdot d \cdot \{|x_i - x_j|; i < j\}_{(k)}$
where $k = \binom{h}{2}$, $h = \lfloor n/2 \rfloor + 1$, and $d = 2.2219$
(consistency constant for Gaussian data). $S_n = c_n' \cdot 1.1926 \cdot
\text{lomed}_i \{\text{himed}_j |x_i - x_j|\}$, where $\text{lomed}$ and
$\text{himed}$ denote the low and high medians respectively.

**Key constants** (full double precision):

| Symbol | Value | Definition |
| :--- | :--- | :--- |
| $\alpha$ | `0.413241928283814` | $\int \psi_{\log}'(u)\,d\Phi(u)$; scoring normalization constant |
| $c$ | `0.37394112142347236` | Solution to $\int \rho_{\log}(u)\,d\Phi(u) = 0.5$; scale rho constant |
| $C_{\text{ADM}}$ | `1.2533141373155001` | $\sqrt{\pi/2}$; ADM consistency constant |
| $C_{\text{MAD}}$ | `1.4826` | $1/\Phi^{-1}(3/4)$; MAD consistency constant |

## Relation to revss and robustbase

This package re-implements the M-estimators from the
[revss](https://CRAN.R-project.org/package=revss) package (Avraham Adler) and
the $Q_n$ and $S_n$ estimators from
[robustbase](https://CRAN.R-project.org/package=robustbase) (Maechler et al.).

The API for the M-estimators is intentionally identical to `revss`: `adm()`,
`robLoc()`, and `robScale()` accept the same arguments and return the same values.
Code that uses `revss` can switch to `robscale` by changing only the `library()` call.

For `qn()` and `sn()`, the function signatures match `robustbase::Qn()` and
`robustbase::Sn()` (with lowercase names for consistency).

Users who do not need compiled performance—or who prefer a dependency-free
pure-R package—should use `revss` or `robustbase` directly. Both are mature,
well-tested, and widely available.

## References

Rousseeuw, P.J. and Verboven, S. (2002). Robust estimation in very small
samples. *Computational Statistics & Data Analysis*, **40**(4), 741--758.
[doi:10.1016/S0167-9473(02)00078-6](https://doi.org/10.1016/S0167-9473(02)00078-6)

Rousseeuw, P. J., and Croux, C. (1993). Alternatives to the Median Absolute Deviation.
*Journal of the American Statistical Association*, **88**, 1273--1283.
[doi:10.1080/01621459.1993.10476408](https://doi.org/10.1080/01621459.1993.10476408)

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
