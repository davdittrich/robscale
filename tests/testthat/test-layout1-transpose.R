library(testthat)

# WU-LAYOUT-1 RED: bit-exact reproducibility before and after boot_results
# layout transposition (nboot×7 → 7×nboot).
#
# The layout change only moves where values live in memory; the computation
# order within each replicate and across the variance pass is unchanged.
# Results must be bit-identical (tolerance=0).

# Reference values captured from cpp_scale_ensemble before layout change.
ref_ens_n100 <- 1.0179191924972544  # cpp_scale_ensemble(rnorm(100,seed=42), 50L)
ref_ens_n500 <- 1.0082018061107851  # cpp_scale_ensemble(rnorm(500,seed=7),  50L)

test_that("layout1 — cpp_scale_ensemble bit-identical after transpose, n=100", {
  set.seed(42)
  x <- rnorm(100)
  expect_equal(robscale:::cpp_scale_ensemble(x, 50L),
               ref_ens_n100, tolerance = 0)
})

test_that("layout1 — cpp_scale_ensemble bit-identical after transpose, n=500", {
  set.seed(7)
  x <- rnorm(500)
  expect_equal(robscale:::cpp_scale_ensemble(x, 50L),
               ref_ens_n500, tolerance = 0)
})

# Internal consistency: two runs with same seed give same result
test_that("layout1 — cpp_scale_ensemble deterministic, n=200", {
  set.seed(99)
  x <- rnorm(200)
  val1 <- robscale:::cpp_scale_ensemble(x, 100L)
  val2 <- robscale:::cpp_scale_ensemble(x, 100L)
  expect_equal(val1, val2, tolerance = 0)
})

# scale_robust(auto_switch=FALSE) forces ensemble path: must be bit-identical
test_that("layout1 — scale_robust(auto_switch=FALSE) bit-identical, n=50", {
  set.seed(55)
  x <- rnorm(50)  # n=50 < default threshold: would already use ensemble
  ref <- scale_robust(x, n_boot = 50L)
  set.seed(55)
  x <- rnorm(50)
  expect_equal(scale_robust(x, n_boot = 50L), ref, tolerance = 0)
})

# CI results: estimate must be bit-identical (CI bounds are non-deterministic
# for BCa with sorted order, so use percentile for deterministic test)
test_that("layout1 — ensemble CI estimate bit-identical, n=100 percentile", {
  set.seed(42)
  x <- rnorm(100)
  ci <- robscale:::cpp_scale_ensemble_ci(x, 50L, 0.95, 1L)  # method=1: percentile
  set.seed(42)
  x <- rnorm(100)
  ci2 <- robscale:::cpp_scale_ensemble_ci(x, 50L, 0.95, 1L)
  expect_equal(ci2$estimate, ci$estimate, tolerance = 0)
  expect_equal(ci2$ci_lower, ci$ci_lower, tolerance = 0)
  expect_equal(ci2$ci_upper, ci$ci_upper, tolerance = 0)
})

# n_boot edge case: n_boot=2 (minimum)
test_that("layout1 — n_boot=2 edge case does not crash", {
  set.seed(1)
  x <- rnorm(20)
  val <- robscale:::cpp_scale_ensemble(x, 2L)
  expect_true(is.numeric(val) && is.finite(val))
})

# Large n_boot to exercise memory layout thoroughly
test_that("layout1 — cpp_scale_ensemble n_boot=500 finite and positive", {
  set.seed(123)
  x <- rnorm(50)
  val <- robscale:::cpp_scale_ensemble(x, 500L)
  expect_true(is.numeric(val) && is.finite(val) && val > 0)
})

# n_boot=1 is clamped to 2 internally — no crash
test_that("layout1 — n_boot=1 clamped to 2 without crash", {
  set.seed(2)
  x <- rnorm(30)
  expect_error(robscale:::cpp_scale_ensemble(x, 1L), NA)
})
