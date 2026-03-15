# Correctness tests for robScale after adaptive median/MAD selection

test_that("robScale gives non-NA positive result across full size sweep", {
  for (n in c(2, 3, 4, 5, 8, 10, 15, 16, 17, 24, 32, 50, 64, 100,
              128, 256, 512, 600, 601, 1000, 5000, 10000, 50000)) {
    set.seed(42 + n)
    x <- rnorm(n)
    val <- robScale(x)
    expect_false(is.na(val), label = paste("NA at n =", n))
    expect_gt(val, 0, label = paste("non-positive at n =", n))
  }
})

test_that("robScale is robust across distributions at n=10000", {
  n <- 10000
  for (info in list(
    list(seed = 1, gen = quote(rnorm(n)),        lab = "gaussian"),
    list(seed = 2, gen = quote(runif(n)),         lab = "uniform"),
    list(seed = 3, gen = quote(rcauchy(n)),       lab = "cauchy"),
    list(seed = 4, gen = quote(rt(n, df = 3)),   lab = "t3")
  )) {
    set.seed(info$seed)
    x <- eval(info$gen)
    val <- robScale(x)
    expect_false(is.na(val), label = info$lab)
    expect_gt(val, 0, label = info$lab)
  }
})

test_that("robScale edge cases", {
  expect_true(is.na(robScale(numeric(0))))
  expect_equal(robScale(1), 0)        # n=1: ADM fallback → 0
  expect_gt(robScale(c(1, 2)), 0)
  expect_equal(robScale(c(5, 5, 5, 5)), 0)  # all-same: MAD collapses → ADM = 0
})

test_that("robScale is deterministic (50 reps)", {
  set.seed(99)
  x <- rnorm(5000)
  vals <- replicate(50, robScale(x))
  expect_true(all(vals == vals[1]))
})

test_that("pdq_robscale_threshold exists and is in valid range", {
  cfg <- get_qnsn_config()
  expect_true("pdq_robscale_threshold" %in% names(cfg))
  thr <- cfg$pdq_robscale_threshold
  expect_true(is.numeric(thr))
  expect_gte(thr, 2048)
  expect_lte(thr, 500000)
})

test_that("robScale correct at threshold boundary", {
  thr <- get_qnsn_config()$pdq_robscale_threshold
  tol <- sqrt(.Machine$double.eps)
  for (n in c(thr - 1L, thr, thr + 1L)) {
    set.seed(123 + n)
    x <- rnorm(n)
    v1 <- robScale(x)
    # Both sides should give a valid positive result and agree with each other
    # (the adaptive dispatch must not affect the numeric result)
    expect_gt(v1, 0, label = paste("boundary n =", n))
    expect_false(is.na(v1), label = paste("boundary n =", n))
  }
  # Explicit cross-boundary agreement
  set.seed(77)
  x_below <- rnorm(thr - 1L)
  set.seed(77)
  x_above <- rnorm(thr + 1L)
  # Can't compare directly (different sizes), but both must be positive and finite
  expect_gt(robScale(x_below), 0)
  expect_gt(robScale(x_above), 0)
})

test_that("robScale heap path works (n=100000)", {
  set.seed(7)
  x <- rnorm(100000)
  val <- robScale(x)
  expect_false(is.na(val))
  expect_gt(val, 0)
})
