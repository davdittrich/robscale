# Extended correctness tests for IQR after pdqselect migration

K_IQR <- 0.741301109252801

test_that("iqr_scaled matches stats::IQR across full size sweep", {
  for (n in c(2, 3, 4, 5, 8, 10, 15, 16, 17, 24, 32, 50, 64, 100,
              128, 256, 512, 600, 601, 1000, 5000, 10000, 50000)) {
    set.seed(42 + n)
    x <- rnorm(n)
    expected <- IQR(x, type = 7) * K_IQR
    actual <- iqr_scaled(x)
    expect_equal(actual, expected, tolerance = sqrt(.Machine$double.eps),
                 label = paste("n =", n))
  }
})

test_that("iqr_scaled is robust across distributions", {
  n <- 10000
  tol <- sqrt(.Machine$double.eps)
  for (info in list(
    list(seed = 1, gen = quote(rnorm(n)),           lab = "gaussian"),
    list(seed = 2, gen = quote(rt(n, df = 3)),      lab = "t3"),
    list(seed = 3, gen = quote(runif(n)),            lab = "uniform"),
    list(seed = 4, gen = quote(rexp(n))  ,           lab = "exponential")
  )) {
    set.seed(info$seed)
    x <- eval(info$gen)
    expect_equal(iqr_scaled(x, constant = 1), IQR(x, type = 7),
                 tolerance = tol, label = info$lab)
  }
})

test_that("iqr_scaled edge cases", {
  expect_equal(iqr_scaled(numeric(0)), NA_real_)
  expect_equal(iqr_scaled(1), 0)
  expect_equal(iqr_scaled(c(5, 5, 5, 5)), 0)
  expect_equal(iqr_scaled(c(1, 2), constant = 1), IQR(c(1, 2), type = 7),
               tolerance = sqrt(.Machine$double.eps))
  expect_true(iqr_scaled(1:10) > 0)
})

test_that("iqr_scaled is deterministic (50 reps)", {
  set.seed(99)
  x <- rnorm(10000)
  vals <- replicate(50, iqr_scaled(x))
  expect_true(all(vals == vals[1]))
})

test_that("iqr_scaled heap path works (n=100000)", {
  set.seed(7)
  x <- rnorm(100000)
  expected <- IQR(x, type = 7) * K_IQR
  actual <- iqr_scaled(x)
  expect_equal(actual, expected, tolerance = sqrt(.Machine$double.eps))
})
