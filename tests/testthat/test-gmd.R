test_that("gmd() returns correct values for basic cases", {
  x <- c(1, 2, 4, 8, 16)
  expect_no_error(g <- gmd(x))
  expect_true(is.numeric(g))
  expect_true(g > 0)

  # Manual: sorted = 1,2,4,8,16, n=5
  # sum = (-4)*1 + (-2)*2 + 0*4 + 2*8 + 4*16 = -4 -4 +0 +16 +64 = 72
  # raw_gmd = 2 * 72 / (5 * 4) = 7.2
  # scaled = 7.2 * sqrt(pi)/2
  K_GMD <- 0.886226925452758
  expect_equal(gmd(x, constant = 1), 7.2)
  expect_equal(gmd(x), 7.2 * K_GMD)
})

test_that("gmd() estimates sigma under normality", {
  set.seed(42)
  y <- rnorm(5000)
  expect_equal(gmd(y), 1.0, tolerance = 0.05)
})

test_that("gmd() handles edge cases", {
  expect_equal(gmd(numeric(0)), NA_real_)
  expect_equal(gmd(1), 0)
  expect_equal(gmd(c(3, 3, 3)), 0)
  expect_true(gmd(c(1, 2)) > 0)
})

test_that("gmd() handles NA correctly", {
  x <- c(1, 2, 4, NA, 16)
  expect_error(gmd(x), "NAs")
  expect_equal(gmd(x, na.rm = TRUE), gmd(x[!is.na(x)]))
})

test_that("gmd() custom constant works", {
  x <- c(1, 2, 3, 5, 7)
  expect_equal(gmd(x, constant = 1) * 0.886226925452758, gmd(x))
})

test_that("gmd() works for integer input", {
  expect_no_error(gmd(1:10))
  expect_true(gmd(1:10) > 0)
})
