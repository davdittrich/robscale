# test-robScale-nr-production.R — TDD production RED/GREEN tests for NR deployment.
#
# WU-PROD-1: Tests 1, 2, 5 are RED before WU-PROD-2 (robScale uses Aitken;
# Aitken and NR differ by up to ~6×10⁻⁸, so 1e-9 tolerance catches the gap).
# Tests 3, 4 may be GREEN before WU-PROD-2.
#
# After WU-PROD-2: all 5 tests GREEN (robScale uses NR; oracle = NR impl).
# After WU-PROD-4: Tests 1, 2, 5 oracle replaced with pinned numeric values.
#
# Worst-case input (diff ≈ 5.97×10⁻⁸ Aitken vs NR):
#   c(2.048861, -0.954442, 0.28295, 0.20733)
# Empirically found: n=4 inputs with 2 sign changes maximise Aitken/NR divergence.

test_that("Test 1: robScale tight match to NR on worst-case 4-element input", {
  x <- c(2.048861, -0.954442, 0.28295, 0.20733)
  nr  <- robscale:::rob_scale_nr_impl(x, maxit = 80L, tol = sqrt(.Machine$double.eps))
  ref <- robScale(x)
  # RED before WU-PROD-2: diff ≈ 5.97×10⁻⁸; GREEN after (ref IS nr)
  expect_equal(ref, nr, tolerance = 1e-9)
})

test_that("Test 2: robScale tight match to NR on random seed-837 4-element input", {
  set.seed(837); x <- rnorm(4)
  nr  <- robscale:::rob_scale_nr_impl(x, maxit = 80L, tol = sqrt(.Machine$double.eps))
  ref <- robScale(x)
  # RED before WU-PROD-2: diff ≈ same order as Test 1; GREEN after
  expect_equal(ref, nr, tolerance = 1e-9)
})

test_that("Test 3: robScale matches NR on seed-42 n=20 input (regression guard)", {
  set.seed(42); x <- rnorm(20)
  nr  <- robscale:::rob_scale_nr_impl(x, maxit = 80L, tol = sqrt(.Machine$double.eps))
  ref <- robScale(x)
  # Diff ≈ 7×10⁻¹⁵ for this input — may already pass before WU-PROD-2;
  # retained as regression guard after NR deployment.
  expect_equal(ref, nr, tolerance = 1e-9)
})

test_that("Test 4: robScale fallback (adm_core) still fires for near-degenerate input", {
  # MAD of c(5,5,5,5,6) is 0 → s_init = 0 → adm_core fallback invoked.
  # This path bypasses nr_scale_compute entirely; behaviour must be unchanged.
  result <- robScale(c(5, 5, 5, 5, 6))
  expect_true(length(result) == 1L && is.finite(result) && result > 0)
})

test_that("Test 5: robScale value pin — seed-42 n=20 against NR oracle", {
  set.seed(42); x <- rnorm(20)
  nr  <- robscale:::rob_scale_nr_impl(x, maxit = 80L, tol = sqrt(.Machine$double.eps))
  ref <- robScale(x)
  # RED before WU-PROD-2 only if Aitken/NR differ at 1e-9 for this seed
  # (diff ≈ 7×10⁻¹⁵; this test may already GREEN — see Test 3).
  # After WU-PROD-4: oracle replaced with pinned numeric value from WU-PROD-2.
  expect_equal(ref, nr, tolerance = 1e-9)
})
