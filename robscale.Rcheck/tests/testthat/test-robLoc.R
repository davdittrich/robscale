test_that("robLoc computes correctly", {
  tol <- sqrt(.Machine$double.eps)
  set.seed(54321)
  x5 <- runif(5, 0, 100)
  
  robLocTest <- function(x, na.rm = FALSE, tol = sqrt(.Machine$double.eps)) {
    if (na.rm) x <- x[!is.na(x)]
    if (length(x) <= 3) return(median(x))
    obj <- function(t, data) {
      sum((2 * plogis((data - t) / mad(data)) - 1))^2
    }
    fit <- optimize(f = obj, interval = range(x), data = x, tol = tol)
    fit$minimum
  }
  
  expect_equal(robLoc(x5), robLocTest(x5), tolerance = tol)
  expect_equal(robLoc(c(1, 9, 7)), median(c(1, 9, 7)), tolerance = tol)
})

test_that("robLoc handles known scale correctly", {
  tol <- sqrt(.Machine$double.eps)
  y <- c(9, 2, 14, 4)
  
  robLocScaleTest <- function(x, scale, na.rm = FALSE,
                              tol = sqrt(.Machine$double.eps)) {
    if (na.rm) x <- x[!is.na(x)]
    if (length(x) <= 2) return(median(x))
    obj <- function(t, data) {
      sum((2 * plogis((data - t) / scale) - 1))^2
    }
    fit <- optimize(f = obj, interval = range(x), data = x, tol = tol)
    fit$minimum
  }
  
  expect_equal(robLoc(y, scale = 5), robLocScaleTest(y, scale = 5), tolerance = tol)
  expect_equal(robLoc(c(1, 8, 12), scale = 5),
               robLocScaleTest(c(1, 8, 12), scale = 5), tolerance = tol)
  expect_equal(robLoc(c(1, 8), scale = 5), median(c(1, 8)), tolerance = tol)
})

test_that("robLoc handles NAs and edge cases", {
  tol <- sqrt(.Machine$double.eps)
  expect_error(robLoc(c(1, 2, NA)), "There are NAs in the data yet na.rm is FALSE")
  expect_true(is.na(robLoc(numeric(0))))
  expect_equal(robLoc(5), 5)
  expect_equal(robLoc(c(5, 5, 5, 5)), 5)
})
