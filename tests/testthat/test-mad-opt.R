# tests/testthat/test-mad-opt.R
#
# TDD gate tests for OPT-M1..M7 (mad_scaled() performance optimizations).
# Tests 1.1-1.3, 1.5 are RED before Phase 0 pin is captured (1.1, 1.5) or
# because mad_abs_diff_n4 does not exist (1.4).
# All must be GREEN after Phase 8.

K_MAD <- 1.482602218505602

# Test 1.1: Regression guard — pin-exact result before and after optimization
test_that("mad_scaled pin unchanged after OPT-M1..M7", {
  tryCatch(devtools::load_all(quiet = TRUE), error = function(e) invisible(NULL))
  pin <- readRDS(testthat::test_path("fixtures/mad_baseline_pin.rds"))

  set.seed(99)
  expect_equal(mad_scaled(rnorm(500)), pin$mad_auto_500,
               tolerance = 8 * .Machine$double.eps)

  set.seed(42)
  expect_equal(mad_scaled(rnorm(100), center = 0), pin$mad_center_100,
               tolerance = 8 * .Machine$double.eps)
})

# Test 1.2: NOINLINE frame boundary — n=64 and n=65 both correct (OPT-M1)
test_that("NOINLINE frame boundary: n=64 and n=65 give correct results", {
  set.seed(7)
  x64 <- rnorm(64)
  x65 <- c(x64, rnorm(1))

  expect_equal(mad_scaled(x64),
               stats::mad(x64, constant = 1) * K_MAD,
               tolerance = sqrt(.Machine$double.eps),
               label = "n=64")
  expect_equal(mad_scaled(x65),
               stats::mad(x65, constant = 1) * K_MAD,
               tolerance = sqrt(.Machine$double.eps),
               label = "n=65")
})

# Test 1.3: Known-center out-of-place RESTRICT path correctness (OPT-M3)
test_that("mad_impl_center RESTRICT path: known exact result and stats::mad agreement", {
  # |c(1,2,5,10) - 4| = c(3,2,1,6), sorted = c(1,2,3,6), median = (2+3)/2 = 2.5
  expect_equal(mad_scaled(c(1, 2, 5, 10), center = 4, constant = 1.0),
               2.5)

  set.seed(42)
  x <- rnorm(20)
  expect_equal(mad_scaled(x, center = 0),
               stats::mad(x, center = 0, constant = K_MAD),
               tolerance = sqrt(.Machine$double.eps))
})

# Test 1.4: bulk_abs_diff kernel — verified via direct call during Phase 4b and
# confirmed GREEN before Phase 8. Diagnostic export removed in Phase 8.
# Kernel correctness is now tested transitively through tests 1.1–1.3.

# Test 1.5: Ensemble pin — bit-identical result after MAD block optimization (OPT-M6/M7)
# n=10 forces the ensemble path (below the auto_switch threshold of 20).
# Bit-identical: XorShift32 per replicate; output is deterministic regardless of TBB scheduling.
test_that("ensemble scale_robust pin unchanged after ensemble MAD optimization", {
  tryCatch(devtools::load_all(quiet = TRUE), error = function(e) invisible(NULL))
  pin <- readRDS(testthat::test_path("fixtures/mad_baseline_pin.rds"))
  set.seed(77)
  result <- scale_robust(rnorm(10), n_boot = 50)
  expect_equal(result, pin$ens_n10, tolerance = 8 * .Machine$double.eps)
})
