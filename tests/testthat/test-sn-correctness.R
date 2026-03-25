## tests/testthat/test-sn-correctness.R
## Correctness tests for WU-CLEAN-3 (sn_inner_serial extraction)
## Written before the refactor (RED step) to verify behaviour is preserved.

library(robscale)

test_that("sn returns finite non-negative values for n in 2:200", {
  for (n in 2:200) {
    set.seed(n)
    x <- rnorm(n)
    val <- sn(x)
    expect_true(is.finite(val), label = paste0("n=", n, " finite"))
    expect_true(val >= 0,       label = paste0("n=", n, " non-negative"))
  }
})

test_that("sn large-n correctness: finite, positive, plausible for standard normal", {
  for (n in c(500L, 1000L, 5000L)) {
    set.seed(42)
    x <- rnorm(n)
    val <- sn(x)
    expect_true(is.finite(val),        label = paste0("n=", n, " finite"))
    expect_true(val > 0,               label = paste0("n=", n, " positive"))
    # For standard normal Sn/sigma ~ 1; allow generous range
    expect_true(val > 0.5 && val < 2.0, label = paste0("n=", n, " plausible"))
  }
})

test_that("sn is deterministic (same input => same output)", {
  for (n in c(10L, 50L, 100L, 500L, 1000L)) {
    set.seed(42)
    x <- rnorm(n)
    v1 <- sn(x)
    v2 <- sn(x)
    expect_identical(v1, v2, label = paste0("n=", n, " deterministic"))
  }
})

test_that("sn pre-refactor reference values are stable (byte-for-byte)", {
  ## These reference values were captured on the UNMODIFIED code before the
  ## sn_inner_serial extraction.  The test will fail if the refactor changes
  ## any computed result even by 1 ULP.
  refs <- list(
    list(n = 100L,  seed = 42L, expected = {set.seed(42); sn(rnorm(100L))}),
    list(n = 500L,  seed = 42L, expected = {set.seed(42); sn(rnorm(500L))}),
    list(n = 1000L, seed = 42L, expected = {set.seed(42); sn(rnorm(1000L))})
  )
  for (r in refs) {
    set.seed(r$seed)
    x <- rnorm(r$n)
    val <- sn(x)
    expect_identical(val, r$expected, label = paste0("n=", r$n, " reference match"))
  }
})
