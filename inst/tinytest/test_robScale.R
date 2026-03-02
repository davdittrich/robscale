tol <- sqrt(.Machine$double.eps)

y <- c(9, 2, 14, 4)

## Key reference value
expect_equal(robScale(y), 5.8798343299206977, tolerance = tol)

## n < minobs fallback to mad
expect_equal(robScale(y[1:3]), mad(y[1:3]), tolerance = tol)
expect_equal(robScale(1:3), mad(1:3), tolerance = tol)

## Implosion: mad <= implbound → fallback to adm
expect_equal(robScale(c(0.00001, 0, 4)), adm(c(0.00001, 0, 4)), tolerance = tol)
expect_equal(robScale(c(0.0001, 0, 4)), mad(c(0.0001, 0, 4)), tolerance = tol)

## n=4 with near-implosion
expect_equal(robScale(c(1e-4, 0, 0, 4)), 0.000101530115510382, tolerance = 1e-7)

## Known location
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
expect_equal(robScale(1:3, loc = 3), robScaleLocTest(1:3, loc = 3), tolerance = tol)

## Known location differs from unknown location
expect_false(isTRUE(all.equal(robScale(1:3), robScaleLocTest(1:3, loc = 0))))

## NA error trapping
naErr <- "There are NAs in the data yet na.rm is FALSE"
set.seed(99)
x5 <- runif(5, 0, 100)
expect_error(robScale(c(x5, NA)), pattern = naErr)
expect_equal(robScale(c(x5, NA), na.rm = TRUE), robScale(x5), tolerance = tol)

## Edge cases: n=0, n=1, n=2, all-identical
expect_true(is.na(robScale(numeric(0))))
expect_equal(robScale(5), 0)
expect_equal(robScale(c(3, 7)), mad(c(3, 7)), tolerance = tol)
expect_equal(robScale(c(5, 5, 5, 5)), 0)
expect_equal(robScale(c(5, 5, 5, 5, 5)), 0)
expect_equal(robScale(c(5, 5, 5, 5), loc = 5), 0)

## maxit/tol: single iteration differs from converged result
full <- robScale(y)
partial <- robScale(y, maxit = 1L)
expect_false(isTRUE(all.equal(partial, full)))
expect_equal(robScale(y, maxit = 0L), mad(y), tolerance = tol)
