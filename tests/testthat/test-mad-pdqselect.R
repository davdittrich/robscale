# Extended correctness tests for MAD after pdqselect migration

K_MAD <- 1.482602218505602

test_that("mad_scaled matches stats::mad across full size sweep", {
  for (n in c(2, 3, 4, 5, 8, 10, 15, 16, 17, 24, 32, 50, 64, 100,
              128, 256, 512, 600, 601, 1000, 5000, 10000, 50000)) {
    set.seed(42 + n)
    x <- rnorm(n)
    expected <- stats::mad(x, constant = 1) * K_MAD
    actual <- mad_scaled(x)
    expect_equal(actual, expected, tolerance = sqrt(.Machine$double.eps),
                 label = paste("n =", n))
  }
})

test_that("mad_scaled is robust across distributions", {
  n <- 10000
  tol <- sqrt(.Machine$double.eps)
  for (info in list(
    list(seed = 1, gen = quote(rnorm(n)),           lab = "gaussian"),
    list(seed = 2, gen = quote(rt(n, df = 3)),      lab = "t3"),
    list(seed = 3, gen = quote(runif(n)),            lab = "uniform"),
    list(seed = 4, gen = quote(rexp(n)),             lab = "exponential")
  )) {
    set.seed(info$seed)
    x <- eval(info$gen)
    expect_equal(mad_scaled(x, constant = 1),
                 stats::mad(x, constant = 1),
                 tolerance = tol, label = info$lab)
  }
})

test_that("mad_scaled edge cases", {
  expect_equal(mad_scaled(numeric(0)), NA_real_)
  expect_equal(mad_scaled(1), 0)
  expect_equal(mad_scaled(c(5, 5, 5, 5)), 0)
  expect_equal(mad_scaled(c(1, 2), constant = 1),
               stats::mad(c(1, 2), constant = 1),
               tolerance = sqrt(.Machine$double.eps))
  expect_true(mad_scaled(1:10) > 0)
})

test_that("mad_scaled with center parameter across sizes", {
  for (n in c(5, 10, 50, 100, 1000, 10000)) {
    set.seed(42 + n)
    x <- rnorm(n)
    center <- 0.5
    expect_equal(mad_scaled(x, center = center),
                 stats::mad(x, center = center, constant = 1) * K_MAD,
                 tolerance = sqrt(.Machine$double.eps),
                 label = paste("center, n =", n))
  }
})

test_that("mad_scaled is deterministic (50 reps)", {
  set.seed(99)
  x <- rnorm(10000)
  vals <- replicate(50, mad_scaled(x))
  expect_true(all(vals == vals[1]))
})

test_that("pdq_median_threshold exists and is reasonable", {
  cfg <- get_qnsn_config()
  expect_true("pdq_median_threshold" %in% names(cfg))
  thr <- cfg$pdq_median_threshold
  expect_true(is.numeric(thr))
  expect_true(thr >= 2048)
  expect_true(thr <= 200000)
})

test_that("mad_scaled correct at threshold boundary", {
  thr <- get_qnsn_config()$pdq_median_threshold
  tol <- sqrt(.Machine$double.eps)
  for (n in c(thr - 1L, thr, thr + 1L)) {
    set.seed(123 + n)
    x <- rnorm(n)
    expect_equal(mad_scaled(x, constant = 1),
                 stats::mad(x, constant = 1),
                 tolerance = tol,
                 label = paste("boundary n =", n))
  }
})

test_that("mad_scaled correct at crossover sizes", {
  tol <- sqrt(.Machine$double.eps)
  for (n in c(3000, 5000, 10000, 15000, 20000)) {
    set.seed(42 + n)
    x <- rnorm(n)
    expect_equal(mad_scaled(x, constant = 1),
                 stats::mad(x, constant = 1),
                 tolerance = tol,
                 label = paste("crossover n =", n))
  }
})

test_that("mad_scaled heap path works (n=100000)", {
  set.seed(7)
  x <- rnorm(100000)
  expected <- stats::mad(x, constant = 1) * K_MAD
  actual <- mad_scaled(x)
  expect_equal(actual, expected, tolerance = sqrt(.Machine$double.eps))
})
