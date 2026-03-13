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
