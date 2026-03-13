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
