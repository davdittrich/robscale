# robscale 0.1.0

Initial release.

* Three functions: `adm()`, `robLoc()`, `robScale()` implementing the robust
  location and scale M-estimators of Rousseeuw & Verboven (2002) for very small
  samples.
* API-compatible drop-in replacement for the `revss` package.
* C++17 implementation via Rcpp with:
    - Newton--Raphson iteration for location (quadratic convergence).
    - Algebraic `tanh(x/2)` identity for the logistic psi function with
      platform-vectorized bulk evaluation (Apple Accelerate on macOS,
      OpenMP SIMD hints on Linux).
    - Floyd--Rivest O(n) selection for median computation.
    - Optimal sorting networks for n <= 8.
    - Stack-allocated arena buffers (zero heap allocation for n <= 512).
* Numerical equivalence with `revss` verified across 5,400 systematic
  comparisons (tolerance: sqrt(.Machine$double.eps)).
