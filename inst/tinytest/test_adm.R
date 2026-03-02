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

## NA error when na.rm = FALSE (default)
naErr <- "There are NAs in the data yet na.rm is FALSE"
expect_error(adm(c(x5, NA)), pattern = naErr)
expect_error(adm(c(x5, NA), constant = 1), pattern = naErr)

## Known values
expect_equal(adm(c(1:9)), sqrt(pi / 2) * mean(abs(1:9 - 5)), tolerance = tol)

## Explicit center
x6 <- c(1, 2, 3, 5, 7, 8)
expect_equal(adm(x6, center = 4.0),
             sqrt(pi / 2) * mean(abs(x6 - 4.0)), tolerance = tol)
expect_equal(adm(x6, center = 4.0, constant = 1),
             mean(abs(x6 - 4.0)), tolerance = tol)

## Sorting networks: correct median for reverse-sorted input at each n=2..8
for (n in 2:8) {
  x_rev <- as.double(n:1)
  expect_equal(adm(x_rev, constant = 1),
               mean(abs(x_rev - median(x_rev))),
               tolerance = tol,
               info = paste0("sort_net reverse n=", n))
}

## Edge cases: n=0, n=1, n=2, all-identical
expect_true(is.na(adm(numeric(0))))
expect_equal(adm(5), 0)
expect_equal(adm(c(3, 7)), sqrt(pi / 2) * 2, tolerance = tol)
expect_equal(adm(c(5, 5, 5, 5)), 0)
