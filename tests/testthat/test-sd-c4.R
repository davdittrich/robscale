test_that("sd_c4() corrects for small-sample bias", {
  x <- c(1, 2, 4, 8, 16)
  n <- length(x)
  c4n <- exp(0.5 * log(2 / (n - 1)) + lgamma(n / 2) - lgamma((n - 1) / 2))
  expect_equal(sd_c4(x), sd(x) / c4n)
})

test_that("sd_c4() estimates sigma under normality", {
  set.seed(42)
  y <- rnorm(5000)
  expect_equal(sd_c4(y), 1.0, tolerance = 0.05)
})

test_that("sd_c4() handles edge cases", {
  expect_equal(sd_c4(numeric(0)), NA_real_)
  expect_true(is.na(sd_c4(1)))
  expect_equal(sd_c4(c(3, 3, 3)), 0)
  expect_true(sd_c4(c(1, 2)) > 0)
})

test_that("sd_c4() handles NA correctly", {
  x <- c(1, 2, 4, NA, 16)
  expect_error(sd_c4(x), "NAs")
  expect_equal(sd_c4(x, na.rm = TRUE), sd_c4(x[!is.na(x)]))
})

test_that("sd_c4() works for integer input", {
  expect_no_error(sd_c4(1:10))
  expect_true(sd_c4(1:10) > 0)
})

test_that("sd_c4() converges to sd() for large n", {
  set.seed(7)
  x <- rnorm(10000)
  # c4(10000) is very close to 1, so sd_c4 ≈ sd
  expect_equal(sd_c4(x), sd(x), tolerance = 1e-4)
})

test_that("sd_c4() is numerically stable for large offsets", {
  # Welford should handle this correctly; naive formula would fail
  x <- c(1e12, 1e12 + 1, 1e12 + 2)
  expect_true(sd_c4(x) > 0)
  expect_equal(sd_c4(x), sd(x) / exp(0.5 * log(2/2) + lgamma(1.5) - lgamma(1)), tolerance = 1e-10)
})
