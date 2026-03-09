test_that("robScale computes correctly", {
  tol <- sqrt(.Machine$double.eps)
  y <- c(9, 2, 14, 4)
  expect_equal(robScale(y), 5.8798343299206977, tolerance = tol)
  expect_equal(robScale(y[1:3]), mad(y[1:3]), tolerance = tol)
})

test_that("robScale handles implosion and fallback", {
  tol <- sqrt(.Machine$double.eps)
  expect_equal(robScale(c(0.00001, 0, 4)), adm(c(0.00001, 0, 4)), tolerance = tol)
  expect_equal(robScale(c(1e-4, 0, 0, 4)), 0.000101530115510382, tolerance = 1e-7)
})

test_that("robScale handles known location", {
  tol <- sqrt(.Machine$double.eps)
  y <- c(9, 2, 14, 4)
  
  robScaleLocTest <- function(x, loc) {
    x <- x - loc
    s <- 1.4826 * median(abs(x))
    converged <- FALSE
    k <- 0
    while (!converged && k < 80) {
      k <- k + 1
      v <- sqrt(2 * mean((2 * plogis(x / (s * 0.37394112142347236)) - 1)^2))
      converged <- abs(v - 1) <= sqrt(.Machine$double.eps)
      s <- s * v
    }
    s
  }
  
  expect_equal(robScale(y, loc = 7), robScaleLocTest(y, loc = 7), tolerance = tol)
})

test_that("robScale handles NAs and edge cases", {
  expect_error(robScale(c(1, 2, NA)), "There are NAs in the data yet na.rm is FALSE")
  expect_true(is.na(robScale(numeric(0))))
  expect_equal(robScale(5), 0)
  expect_equal(robScale(c(5, 5, 5, 5)), 0)
})
