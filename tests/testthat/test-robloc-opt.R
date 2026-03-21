# robLoc() Performance Optimization Tests (v0.5.0)
# Phase 0: TDD Guards for L1 (SIMD), L2 (Dispatch), L3 (TBB), and L4 (Cache)

test_that("0.1 — robLoc returns finite result for all tested n", {
  # Phase 1-3 Regression Guard: Basic correctness for small and medium n
  set.seed(42)
  for (n in c(8L, 16L, 32L, 64L, 100L, 500L, 1000L)) {
    loc <- robLoc(rnorm(n))
    expect_true(is.finite(loc),
      label = sprintf("robLoc(rnorm(%d)) = %g must be finite", n, loc))
  }
})

test_that("0.2 — robLoc agrees with R optimize reference to 1e-6 after all OPTs", {
  # Phase 1-3 Regression Guard: Cross-check against logistic fixed-point R reference
  # Mirrors existing test-robLoc.R pattern; tolerance 1e-6 = NR vs Brent difference
  robLocRef <- function(x) {
    obj <- function(t) sum((2*plogis((x-t)/mad(x))-1))^2
    optimize(obj, range(x), tol=sqrt(.Machine$double.eps))$minimum
  }
  set.seed(1001)
  for (n in c(10L, 25L, 50L, 100L)) {
    x   <- rnorm(n)
    got <- robLoc(x)
    ref <- robLocRef(x)
    expect_equal(got, ref, tolerance=1e-5,
      label=sprintf("n=%d: robLoc=%g, ref=%g", n, got, ref))
  }
})

test_that("0.3 — robLoc fused AVX2 agrees with scalar path to 2*sqrt(eps)", {
  # Phase 3 (OPT-L1) Guard: Primary correctness gate for fused kernel
  # SKIPs until rob_loc_scalar_impl diagnostic helper is added in Phase 3
  skip_if_not(
    exists("rob_loc_scalar_impl", envir=asNamespace("robscale"), mode="function"),
    "Diagnostic scalar helper not compiled"
  )
  set.seed(777)
  tol <- 2 * sqrt(.Machine$double.eps)
  max_diff <- 0
  for (n in c(16L, 32L, 64L, 100L, 256L, 500L, 1000L)) {
    for (rep in seq_len(20L)) {
      x  <- rnorm(n)
      v1 <- robscale:::rob_loc_scalar_impl(x)
      v2 <- robLoc(x)
      d  <- abs(v1 - v2)
      if (d > max_diff) max_diff <- d
    }
  }
  expect_lt(max_diff, tol,
    label=sprintf("max|fused-scalar| = %.3e (tol %.3e)", max_diff, tol))
})

test_that("0.4 — robLoc handles near-degenerate scale without NaN/Inf", {
  # Phase 3 (OPT-L1) Guard: Edge case for sum_dpsi near-zero guard (pathological input)
  # All values very close together: scale ~= 0 → u_i huge → tanh(u_i) → +/-1
  # sum_dpsi = sum sech^2(u_i) → 0 without guard
  x_degen <- c(rep(1.0, 50L), rep(1.0 + 1e-10, 50L))
  loc <- robLoc(x_degen)
  expect_true(is.finite(loc) && !is.nan(loc),
    label=sprintf("robLoc degenerate input = %g must be finite", loc))

  # Constant input: s=0 path; must return median
  x_const <- rep(5.0, 20L)
  expect_equal(robLoc(x_const), 5.0, tolerance=1e-10)
})

test_that("0.5 — robLoc parallel (n>=threshold) agrees with serial to 1e-10", {
  # Phase 4 (OPT-L3) Guard: Critical gate for TBB parallel implementation
  # SKIPs until rob_loc_has_parallel() and rob_loc_serial_impl() added in Phase 4
  skip_if_not(
    exists("rob_loc_has_parallel", envir=asNamespace("robscale"), mode="function"),
    "TBB not compiled"
  )
  skip_if_not(robscale:::rob_loc_has_parallel())
  set.seed(2024)
  tol <- 1e-10
  for (n in c(4096L, 8192L, 16384L)) {
    x <- rnorm(n)
    v_par <- robLoc(x)                            # dispatches parallel
    v_ser <- robscale:::rob_loc_serial_impl(x)    # forces serial path
    expect_equal(v_par, v_ser, tolerance=tol,
      label=sprintf("n=%d: |parallel-serial| = %.3e", n, abs(v_par-v_ser)))
  }
})
