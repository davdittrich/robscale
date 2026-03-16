test_that("mad_scaled() matches base mad() closely", {
  set.seed(42)
  for (n in c(5, 10, 20, 50, 100)) {
    x <- rnorm(n)
    # Selection-based vs sort-based median may differ at FP level
    expect_equal(mad_scaled(x), mad(x), tolerance = 1e-5,
                 label = paste("n =", n))
  }
})

test_that("mad_scaled() estimates sigma under normality", {
  set.seed(42)
  y <- rnorm(5000)
  expect_equal(mad_scaled(y), 1.0, tolerance = 0.05)
})

test_that("mad_scaled() handles edge cases", {
  expect_equal(mad_scaled(numeric(0)), NA_real_)
  expect_equal(mad_scaled(1), 0)
  expect_equal(mad_scaled(c(3, 3, 3)), 0)
  expect_true(mad_scaled(c(1, 2)) > 0)
})

test_that("mad_scaled() handles NA correctly", {
  x <- c(1, 2, 4, NA, 16)
  expect_error(mad_scaled(x), "NAs")
  expect_equal(mad_scaled(x, na.rm = TRUE), mad_scaled(x[!is.na(x)]))
})

test_that("mad_scaled() accepts custom center", {
  x <- c(1, 2, 3, 5, 7)
  m <- mad_scaled(x, center = 3.0, constant = 1)
  # Should be median(|x - 3|) = median(2, 1, 0, 2, 4) = 2
  expect_equal(m, 2.0)
})

test_that("mad_scaled() custom constant works", {
  x <- c(1, 2, 3, 5, 7)
  expect_equal(mad_scaled(x, constant = 1) * 1.482602218505602, mad_scaled(x),
               tolerance = 1e-10)
})

test_that("mad_scaled() works for integer input", {
  expect_no_error(mad_scaled(1:10))
  expect_true(mad_scaled(1:10) > 0)
})

test_that("mad_scaled matches stats::mad across all regimes", {
  K_MAD <- 1.482602218505602
  for (n in c(2, 3, 4, 5, 8, 10, 15, 16, 17, 32, 50, 64, 100,
              128, 256, 512, 599, 600, 601, 1000, 5000, 10000)) {
    set.seed(42 + n)
    x <- rnorm(n)
    # Compare with raw MAD (constant=1) to avoid stats::mad's truncated 1.4826
    expected <- stats::mad(x, constant = 1) * K_MAD
    actual <- mad_scaled(x)
    expect_equal(actual, expected, tolerance = sqrt(.Machine$double.eps),
                 label = paste("n =", n))
  }
})

test_that("fused MAD handles ties and degenerate inputs", {
  K_MAD <- 1.482602218505602
  expect_equal(mad_scaled(c(1, 1, 1, 1, 1)), 0)
  expect_equal(mad_scaled(c(1, 2)),
               stats::mad(c(1, 2), constant = 1) * K_MAD,
               tolerance = sqrt(.Machine$double.eps))
  # Even n where median is average of two elements
  expect_equal(mad_scaled(c(1, 2, 3, 4)),
               stats::mad(c(1, 2, 3, 4), constant = 1) * K_MAD,
               tolerance = sqrt(.Machine$double.eps))
  # All identical except one outlier
  expect_equal(mad_scaled(c(rep(5, 99), 100)),
               stats::mad(c(rep(5, 99), 100), constant = 1) * K_MAD,
               tolerance = sqrt(.Machine$double.eps))
})

test_that("mad_scaled with supplied center still works", {
  K_MAD <- 1.482602218505602
  for (n in c(5, 10, 50, 100, 1000)) {
    set.seed(42 + n)
    x <- rnorm(n)
    center <- 0.5
    expect_equal(mad_scaled(x, center = center),
                 stats::mad(x, center = center, constant = 1) * K_MAD,
                 tolerance = sqrt(.Machine$double.eps),
                 label = paste("center, n =", n))
  }
})
