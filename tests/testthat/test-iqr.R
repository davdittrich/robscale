test_that("iqr_scaled() matches base R IQR", {
  x <- c(1, 2, 4, 8, 16)
  K_IQR <- 0.741301109252801
  expect_equal(iqr_scaled(x, constant = 1), IQR(x, type = 7))
  expect_equal(iqr_scaled(x), IQR(x, type = 7) * K_IQR)
})

test_that("iqr_scaled() estimates sigma under normality", {
  set.seed(42)
  y <- rnorm(5000)
  expect_equal(iqr_scaled(y), 1.0, tolerance = 0.05)
})

test_that("iqr_scaled() handles edge cases", {
  expect_equal(iqr_scaled(numeric(0)), NA_real_)
  expect_equal(iqr_scaled(1), 0)
  expect_equal(iqr_scaled(c(3, 3, 3)), 0)
  expect_true(iqr_scaled(c(1, 2)) > 0)
})

test_that("iqr_scaled() handles NA correctly", {
  x <- c(1, 2, 4, NA, 16)
  expect_error(iqr_scaled(x), "NAs")
  expect_equal(iqr_scaled(x, na.rm = TRUE), iqr_scaled(x[!is.na(x)]))
})

test_that("iqr_scaled() custom constant works", {
  x <- c(1, 2, 3, 5, 7)
  expect_equal(iqr_scaled(x, constant = 1) * 0.741301109252801, iqr_scaled(x))
})

test_that("iqr_scaled() works for integer input", {
  expect_no_error(iqr_scaled(1:10))
  expect_true(iqr_scaled(1:10) > 0)
})

test_that("iqr_scaled() matches base IQR across sample sizes", {
  set.seed(99)
  for (n in c(5, 10, 20, 50, 100)) {
    x <- rnorm(n)
    expect_equal(iqr_scaled(x, constant = 1), IQR(x, type = 7),
                 tolerance = 1e-10,
                 label = paste("n =", n))
  }
})

test_that("iqr_scaled matches stats::IQR * constant across all regimes", {
  K_IQR <- 0.741301109252801
  for (n in c(2, 3, 4, 5, 8, 10, 15, 16, 17, 32, 50, 64, 100,
              128, 256, 512, 599, 600, 601, 1000, 5000, 10000)) {
    set.seed(42 + n)
    x <- rnorm(n)
    expected <- IQR(x, type = 7) * K_IQR
    actual <- iqr_scaled(x)
    expect_equal(actual, expected, tolerance = sqrt(.Machine$double.eps),
                 label = paste("n =", n))
  }
})

test_that("iqr_scaled handles ties and degenerate data", {
  K_IQR <- 0.741301109252801
  expect_equal(iqr_scaled(c(1, 1, 1, 1, 1)), 0)
  expect_equal(iqr_scaled(c(1, 2)), IQR(c(1, 2), type = 7) * K_IQR,
               tolerance = sqrt(.Machine$double.eps))
  expect_equal(iqr_scaled(c(1, 1, 1, 2)),
               IQR(c(1, 1, 1, 2), type = 7) * K_IQR,
               tolerance = sqrt(.Machine$double.eps))
  # All identical except one outlier
  expect_equal(iqr_scaled(c(rep(5, 99), 100)),
               IQR(c(rep(5, 99), 100), type = 7) * K_IQR,
               tolerance = sqrt(.Machine$double.eps))
})

test_that("iqr_scaled works at even and odd n near transitions", {
  K_IQR <- 0.741301109252801
  for (n in c(2, 3, 4, 5, 6, 7, 15, 16, 17, 18)) {
    set.seed(123 + n)
    x <- rnorm(n)
    expect_equal(iqr_scaled(x, constant = 1), IQR(x, type = 7),
                 tolerance = sqrt(.Machine$double.eps),
                 label = paste("even/odd n =", n))
  }
})
