tol <- sqrt(.Machine$double.eps)

## Generate test data
set.seed(12345)
x5 <- runif(5, 0, 100)
t5 <- median(x5)
adm5 <- mean(abs(x5 - t5))

## Correctness
expect_equal(adm(x5), adm5 * sqrt(pi / 2), tolerance = tol)

## Custom constant
expect_equal(adm(x5, constant = 1), adm5, tolerance = tol)

## NA handling: na.rm = TRUE strips NAs
expect_equal(adm(c(x5, NA), na.rm = TRUE), adm5 * sqrt(pi / 2), tolerance = tol)
expect_equal(adm(c(x5, NA), constant = 1, na.rm = TRUE), adm5, tolerance = tol)

## NA propagation when na.rm = FALSE (default)
expect_true(is.na(adm(c(x5, NA))))
expect_true(is.na(adm(c(x5, NA), constant = 1)))

## Known values
expect_equal(adm(c(1:9)), sqrt(pi / 2) * mean(abs(1:9 - 5)), tolerance = tol)

## Edge cases: n=0, n=1, n=2, all-identical
expect_true(is.na(adm(numeric(0))))
expect_equal(adm(5), 0)
expect_equal(adm(c(3, 7)), sqrt(pi / 2) * 2, tolerance = tol)
expect_equal(adm(c(5, 5, 5, 5)), 0)
