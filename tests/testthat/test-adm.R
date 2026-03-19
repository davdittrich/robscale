test_that("adm computes correctly", {
  tol <- sqrt(.Machine$double.eps)
  set.seed(12345)
  x5 <- runif(5, 0, 100)
  t5 <- median(x5)
  adm5 <- mean(abs(x5 - t5))
  
  expect_equal(adm(x5), adm5 * sqrt(pi / 2), tolerance = tol)
  expect_equal(adm(x5, constant = 1), adm5, tolerance = tol)
})

test_that("adm handles NAs correctly", {
  tol <- sqrt(.Machine$double.eps)
  set.seed(12345)
  x5 <- runif(5, 0, 100)
  adm5 <- mean(abs(x5 - median(x5)))
  
  expect_equal(adm(c(x5, NA), na.rm = TRUE), adm5 * sqrt(pi / 2), tolerance = tol)
  expect_error(adm(c(x5, NA)), "There are NAs in the data yet na.rm is FALSE")
})

test_that("adm handles explicit center and edge cases", {
  tol <- sqrt(.Machine$double.eps)
  expect_equal(adm(c(1:9)), sqrt(pi / 2) * mean(abs(1:9 - 5)), tolerance = tol)
  
  x6 <- c(1, 2, 3, 5, 7, 8)
  expect_equal(adm(x6, center = 4.0), sqrt(pi / 2) * mean(abs(x6 - 4.0)), tolerance = tol)
  
  expect_true(is.na(adm(numeric(0))))
  expect_equal(adm(5), 0)
  expect_equal(adm(c(5, 5, 5, 5)), 0)
})

test_that("center = NULL is equivalent to omitting center", {
  tol <- sqrt(.Machine$double.eps)
  set.seed(42)
  x <- rnorm(9)
  expect_equal(adm(x, center = NULL), adm(x), tolerance = tol)
})

test_that("sorting networks work for adm", {
  tol <- sqrt(.Machine$double.eps)
  for (n in 2:8) {
    x_rev <- as.double(n:1)
    expect_equal(adm(x_rev, constant = 1),
                 mean(abs(x_rev - median(x_rev))),
                 tolerance = tol,
                 label = paste0("sort_net reverse n=", n))
  }
})

test_that("median_select is correct for n = 17..64 (odd and even)", {
  tol <- sqrt(.Machine$double.eps)
  set.seed(99)
  for (n in 17:64) {
    x <- rnorm(n)
    expect_equal(adm(x, constant = 1),
                 mean(abs(x - median(x))),
                 tolerance = tol,
                 label = paste0("n=", n))
  }
})

test_that("adm is correct for n = 65..200 (AVX2 path)", {
  tol <- sqrt(.Machine$double.eps)
  set.seed(77)
  for (n in c(65, 66, 67, 68, 100, 128, 200)) {
    x <- rnorm(n)
    expect_equal(adm(x, constant = 1),
                 mean(abs(x - median(x))),
                 tolerance = tol,
                 label = paste0("n=", n))
    # Also with explicit center
    ctr <- median(x)
    expect_equal(adm(x, center = ctr, constant = 1),
                 mean(abs(x - ctr)),
                 tolerance = tol,
                 label = paste0("n=", n, " known center"))
  }
})
