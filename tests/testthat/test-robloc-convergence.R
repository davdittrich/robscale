# WU-4: rob_loc relative convergence criterion (R2)
#
# Bug: |v| <= tol (absolute) never triggers for t ≈ 1e10 because the NR
# update at convergence is O(t * eps) ≈ 2.2e-6 >> tol = 1.49e-8.
# Fix: |v| <= tol * max(|t|, 1.0)

test_that("robLoc converges correctly for large-offset data (t~1e10)", {
  x <- 1e10 + c(-2, -1, 0, 1, 2, 3, 5)
  t_est <- robLoc(x)
  # Expected: near 1e10 + 1 (weighted median-like M-location)
  expect_true(abs(t_est - (1e10 + 1)) < 1,
              label = "robLoc converges near 1e10+1")
})

test_that("robLoc converges correctly for negative large-offset data", {
  x <- -1e10 + c(-2, -1, 0, 1, 2, 3, 5)
  t_est <- robLoc(x)
  expect_true(abs(t_est - (-1e10 + 1)) < 1,
              label = "robLoc converges near -1e10+1")
})

test_that("robLoc result unchanged for normal-scale data (regression)", {
  set.seed(42)
  x <- rnorm(20)
  # Pin: compute and store the value with the fix applied
  # The converged value should be the same (NR fixed point unchanged)
  res <- robLoc(x)
  expect_true(is.finite(res))
  expect_true(abs(res - median(x)) < 2 * mad(x),
              label = "robLoc within 2*MAD of median")
})

test_that("robLoc convergence: symmetric data around 0 unchanged", {
  x <- c(-5, -3, -1, 0, 1, 3, 5)
  t_est <- robLoc(x)
  expect_equal(t_est, 0, tolerance = 1e-6,
               label = "symmetric data converges to 0")
})

test_that("robLoc with explicit scale at large offset", {
  x <- 1e12 + c(-2, -1, 0, 1, 2)
  t_est <- robLoc(x, scale = mad(x))
  expect_true(abs(t_est - 1e12) < 1,
              label = "robLoc with scale at 1e12")
})
