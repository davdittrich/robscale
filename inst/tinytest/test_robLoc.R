tol <- sqrt(.Machine$double.eps)

## Generate test data
set.seed(54321)
x5 <- runif(5, 0, 100)

## Known values via optimize (ground truth)
robLocTest <- function(x, na.rm = FALSE, tol = sqrt(.Machine$double.eps)) {
  if (na.rm) x <- x[!is.na(x)]
  if (length(x) <= 3) return(median(x))
  obj <- function(t, data) {
    sum((2 * plogis((data - t) / mad(data)) - 1))^2
  }
  fit <- optimize(f = obj, interval = range(x), data = x, tol = tol)
  fit$minimum
}

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

## Correctness vs optimization
expect_equal(robLoc(x5), robLocTest(x5), tolerance = tol)

## n < minobs fallback to median
expect_equal(robLoc(c(1, 9, 7)), median(c(1, 9, 7)), tolerance = tol)

## Known scale variant
y <- c(9, 2, 14, 4)
expect_equal(robLoc(y, scale = 5), robLocScaleTest(y, scale = 5), tolerance = tol)
expect_equal(robLoc(c(1, 8, 12), scale = 5),
             robLocScaleTest(c(1, 8, 12), scale = 5), tolerance = tol)

## Known scale with n < 3 falls back to median
expect_equal(robLoc(c(1, 8), scale = 5), median(c(1, 8)), tolerance = tol)

## Known scale actually differs from median
expect_false(isTRUE(all.equal(robLoc(c(1, 8, 12), scale = 5),
                              median(c(1, 8, 12)))))

## NA error trapping
naErr <- "There are NAs in the data yet na.rm is FALSE"
expect_error(robLoc(c(x5, NA)), pattern = naErr)
expect_equal(robLoc(c(x5, NA), na.rm = TRUE), robLoc(x5), tolerance = tol)

## Edge cases: n=0, n=1, n=2, all-identical
expect_true(is.na(robLoc(numeric(0))))
expect_equal(robLoc(5), 5)
expect_equal(robLoc(c(3, 7)), median(c(3, 7)))
expect_equal(robLoc(c(5, 5, 5, 5)), 5)
expect_equal(robLoc(c(5, 5, 5, 5, 5)), 5)
expect_equal(robLoc(c(5, 5, 5, 5), scale = 1), 5)
