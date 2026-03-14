test_that("ensemble is deterministic", {
  set.seed(42)
  x <- rnorm(10)
  e1 <- robscale:::cpp_scale_ensemble(x, 200L)
  e2 <- robscale:::cpp_scale_ensemble(x, 200L)
  expect_identical(e1, e2)
})

test_that("ensemble estimates sigma under normality", {
  set.seed(42)
  y <- rnorm(5000)
  expect_equal(robscale:::cpp_scale_ensemble(y, 200L), 1.0, tolerance = 0.05)
})

test_that("ensemble is more robust than sd for contaminated data", {
  set.seed(42)
  x_clean <- rnorm(20)
  x_dirty <- c(x_clean, 100)  # one extreme outlier
  ens_clean <- robscale:::cpp_scale_ensemble(x_clean, 200L)
  ens_dirty <- robscale:::cpp_scale_ensemble(x_dirty, 200L)
  sd_dirty  <- sd(x_dirty)

  # Ensemble should be much less affected than sd
  expect_true(abs(ens_dirty - ens_clean) < abs(sd_dirty - sd(x_clean)))
})

test_that("ensemble handles edge cases", {
  expect_true(is.na(robscale:::cpp_scale_ensemble(1, 200L)))
  expect_true(is.numeric(robscale:::cpp_scale_ensemble(c(1, 2), 200L)))
  expect_true(robscale:::cpp_scale_ensemble(c(1, 2), 200L) > 0)
})

test_that("ensemble returns finite for small samples", {
  for (n in 2:10) {
    set.seed(n)
    x <- rnorm(n)
    ens <- robscale:::cpp_scale_ensemble(x, 200L)
    expect_true(is.finite(ens), label = paste("n =", n))
    expect_true(ens > 0, label = paste("n =", n))
  }
})

test_that("ensemble handles data with ties", {
  # Most robust estimators collapse, but ensemble should still return a value
  val <- robscale:::cpp_scale_ensemble(c(5, 5, 5, 5, 6), 200L)
  expect_true(is.finite(val))
})

# --- Snapshot and determinism stress tests ---

test_that("ensemble snapshot values are stable", {
  # Tight tolerance (1e-14) accommodates ULP-level differences between

  # debug (-O0) and optimized (-O3) builds caused by SIMD pragma vectorization.
  # Determinism tests below use expect_identical for parallelism safety.
  set.seed(1)
  expect_equal(
    robscale:::cpp_scale_ensemble(rnorm(5), 200L),
    1.033903445605208, tolerance = 1e-14
  )
  set.seed(42)
  expect_equal(
    robscale:::cpp_scale_ensemble(rnorm(10), 200L),
    0.8309502035368722, tolerance = 1e-14
  )
  set.seed(123)
  expect_equal(
    robscale:::cpp_scale_ensemble(rnorm(500), 200L),
    0.9615157755575315, tolerance = 1e-14
  )
})

test_that("ensemble is deterministic over 50 repeats (n=10)", {
  set.seed(99)
  x <- rnorm(10)
  ref <- robscale:::cpp_scale_ensemble(x, 200L)
  for (i in seq_len(50)) {
    expect_identical(robscale:::cpp_scale_ensemble(x, 200L), ref)
  }
})

test_that("ensemble is deterministic for various n", {
  for (n in c(2, 5, 10, 15, 50, 200)) {
    set.seed(n)
    x <- rnorm(n)
    e1 <- robscale:::cpp_scale_ensemble(x, 200L)
    e2 <- robscale:::cpp_scale_ensemble(x, 200L)
    expect_identical(e1, e2, label = paste("n =", n))
  }
})

test_that("ensemble is deterministic with large n_boot", {
  set.seed(7)
  x <- rnorm(10)
  e1 <- robscale:::cpp_scale_ensemble(x, 500L)
  e2 <- robscale:::cpp_scale_ensemble(x, 500L)
  expect_identical(e1, e2)
})

test_that("ensemble is deterministic with large n", {
  set.seed(11)
  x <- rnorm(500)
  e1 <- robscale:::cpp_scale_ensemble(x, 200L)
  e2 <- robscale:::cpp_scale_ensemble(x, 200L)
  expect_identical(e1, e2)
})

test_that("ensemble determinism at threshold boundary (n*n_boot=10000)", {
  set.seed(33)
  x <- rnorm(50)
  ref <- robscale:::cpp_scale_ensemble(x, 200L)
  for (i in seq_len(20)) {
    expect_identical(robscale:::cpp_scale_ensemble(x, 200L), ref)
  }
})

test_that("ensemble determinism above threshold (n=100)", {
  set.seed(44)
  x <- rnorm(100)
  ref <- robscale:::cpp_scale_ensemble(x, 200L)
  for (i in seq_len(20)) {
    expect_identical(robscale:::cpp_scale_ensemble(x, 200L), ref)
  }
})

test_that("ensemble determinism high n_boot at threshold (n=10, n_boot=1000)", {
  set.seed(55)
  x <- rnorm(10)
  e1 <- robscale:::cpp_scale_ensemble(x, 1000L)
  e2 <- robscale:::cpp_scale_ensemble(x, 1000L)
  expect_identical(e1, e2)
})
