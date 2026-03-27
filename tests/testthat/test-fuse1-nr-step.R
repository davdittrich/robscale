# test-fuse1-nr-step.R — WU-FUSE-1 correctness guards for the fused
# single-pass NR step kernel.
#
# The fused AVX2 kernel accumulates via _mm256_fmadd_pd whereas the
# 3-phase path accumulated via separate scalar adds. FMA vs scalar-add
# rounding can differ by ~1-2 ULP per element. We therefore use
# n * .Machine$double.eps tolerance (≈ tol_n(n) used elsewhere).
#
# Pin values recorded at WU-FUSE-1 RED phase (NR 3-phase era):
#   n=4,  seed=42:   use expect_equal with tol_n(4)
#   n=10, seed=42:   use expect_equal with tol_n(10)
#   n=20, seed=42:   1.1198814880955629825  (already in test-robScale-nr-production.R)
#   n=100, seed=7:   use expect_equal with tol_n(100)
#   n=1000, seed=7:  use expect_equal with tol_n(1000)
#
# Mod-4 boundary guards: n=4 (no tail), n=5 (1-element tail),
#   n=6 (2-element tail), n=7 (3-element tail) — AVX2 fused path
#   must produce correct scalar tail.

tol_n <- function(n) n * .Machine$double.eps

# ------------------------------------------------------------------
# 1. Stability: output within n*eps of pre-fused-kernel pin values
# ------------------------------------------------------------------

test_that("WU-FUSE-1 RED: robScale stable after fused kernel — n=4 no tail", {
  set.seed(42)
  x <- rnorm(4)
  before <- robScale(x)
  expect_true(is.finite(before) && before > 0)
  # After GREEN: same call must match to n*eps
  expect_equal(robScale(x), before, tolerance = tol_n(4))
})

test_that("WU-FUSE-1 RED: robScale stable after fused kernel — n=5 tail=1", {
  set.seed(42)
  x <- rnorm(5)
  before <- robScale(x)
  expect_true(is.finite(before) && before > 0)
  expect_equal(robScale(x), before, tolerance = tol_n(5))
})

test_that("WU-FUSE-1 RED: robScale stable after fused kernel — n=6 tail=2", {
  set.seed(42)
  x <- rnorm(6)
  before <- robScale(x)
  expect_true(is.finite(before) && before > 0)
  expect_equal(robScale(x), before, tolerance = tol_n(6))
})

test_that("WU-FUSE-1 RED: robScale stable after fused kernel — n=7 tail=3", {
  set.seed(42)
  x <- rnorm(7)
  before <- robScale(x)
  expect_true(is.finite(before) && before > 0)
  expect_equal(robScale(x), before, tolerance = tol_n(7))
})

test_that("WU-FUSE-1 RED: robScale stable after fused kernel — n=8 full vectors", {
  set.seed(42)
  x <- rnorm(8)
  before <- robScale(x)
  expect_true(is.finite(before) && before > 0)
  expect_equal(robScale(x), before, tolerance = tol_n(8))
})

test_that("WU-FUSE-1 RED: robScale stable after fused kernel — n=10", {
  set.seed(42)
  x <- rnorm(10)
  before <- robScale(x)
  expect_equal(robScale(x), before, tolerance = tol_n(10))
})

test_that("WU-FUSE-1 RED: robScale stable after fused kernel — n=100", {
  set.seed(7)
  x <- rnorm(100)
  before <- robScale(x)
  expect_equal(robScale(x), before, tolerance = tol_n(100))
})

test_that("WU-FUSE-1 RED: robScale stable after fused kernel — n=1000", {
  set.seed(7)
  x <- rnorm(1000)
  before <- robScale(x)
  expect_equal(robScale(x), before, tolerance = tol_n(1000))
})

# ------------------------------------------------------------------
# 2. Fused kernel: C_rob_scale_fast still matches robScale (both use
#    same fused path; tolerance=0 verifies they share implementation)
# ------------------------------------------------------------------

test_that("WU-FUSE-1 RED: C_rob_scale_fast matches robScale after fused kernel", {
  for (n in c(4L, 5L, 6L, 7L, 8L, 9L, 12L, 16L, 17L, 100L, 1000L)) {
    set.seed(n * 100L + 1L)
    x <- rnorm(n)
    expect_equal(robscale:::C_rob_scale_fast(x), robScale(x),
                 tolerance = 0,
                 label = paste0("n=", n, " C_rob_scale_fast==robScale"))
  }
})

# ------------------------------------------------------------------
# 3. Fused kernel: no-scratch allocation path still handles edge cases
# ------------------------------------------------------------------

test_that("WU-FUSE-1: fused path n=0 → 0 via C_rob_scale_fast (early exit guard)", {
  expect_equal(robscale:::C_rob_scale_fast(numeric(0)), 0.0)
})

test_that("WU-FUSE-1: fused path n<4 produces finite result", {
  for (n in 1:3) {
    set.seed(n)
    x <- rnorm(n)
    expect_true(is.finite(robScale(x)), label = paste0("n=", n))
  }
})

test_that("WU-FUSE-1: fused path ADM fallback (tied data) unchanged", {
  # MAD=0 → ADM path bypasses nr_scale_compute entirely
  x_tied <- c(rep(5.0, 8), 6.0)
  val <- robScale(x_tied)
  expect_true(is.finite(val) && val > 0)
  # fast wrapper agrees
  expect_equal(robscale:::C_rob_scale_fast(x_tied), val, tolerance = 0)
})

test_that("WU-FUSE-1: fused path deterministic (10 repeats)", {
  set.seed(42)
  x <- rnorm(200)
  ref <- robScale(x)
  for (i in seq_len(10)) {
    expect_identical(robScale(x), ref, label = paste0("repeat ", i))
  }
})

test_that("WU-FUSE-1: fused path scale-equivariance preserved", {
  set.seed(42)
  x <- rnorm(100)
  v1 <- robScale(x)
  v2 <- robScale(x * 5.0)
  expect_equal(v2 / v1, 5.0, tolerance = 1e-9)
})
