# Exhaustive correctness tests for IQR implementation.
# Written BEFORE any refactor (TDD red step — must be GREEN on unmodified code).
# Verifies: iqr_scaled(x, constant=1) == stats::IQR(x, type=7)
# for all n in 2:200, including small-n edge cases (2, 3, 4, 5).

test_that("iqr matches stats::IQR(type=7) for all n in 2:200 (random normal)", {
  for (n in 2:200) {
    set.seed(n)
    x <- rnorm(n)
    expected <- stats::IQR(x, type = 7)
    actual   <- iqr_scaled(x, constant = 1)
    expect_equal(actual, expected,
                 tolerance = 8 * .Machine$double.eps,
                 label = paste0("n=", n))
  }
})

test_that("iqr matches stats::IQR(type=7) for all n in 2:200 (random uniform)", {
  for (n in 2:200) {
    set.seed(n + 1000L)
    x <- runif(n)
    expected <- stats::IQR(x, type = 7)
    actual   <- iqr_scaled(x, constant = 1)
    expect_equal(actual, expected,
                 tolerance = 8 * .Machine$double.eps,
                 label = paste0("n=", n))
  }
})

test_that("iqr matches stats::IQR(type=7) for small-n edge cases with ties", {
  # n=2: both values identical
  x2_tie <- c(3, 3)
  expect_equal(iqr_scaled(x2_tie, constant = 1), stats::IQR(x2_tie, type = 7),
               tolerance = 8 * .Machine$double.eps, label = "n=2 tie")

  # n=2: distinct values
  x2 <- c(1, 5)
  expect_equal(iqr_scaled(x2, constant = 1), stats::IQR(x2, type = 7),
               tolerance = 8 * .Machine$double.eps, label = "n=2 distinct")

  # n=3: all same
  x3_same <- c(7, 7, 7)
  expect_equal(iqr_scaled(x3_same, constant = 1), stats::IQR(x3_same, type = 7),
               tolerance = 8 * .Machine$double.eps, label = "n=3 same")

  # n=4: one outlier
  x4_out <- c(1, 2, 3, 100)
  expect_equal(iqr_scaled(x4_out, constant = 1), stats::IQR(x4_out, type = 7),
               tolerance = 8 * .Machine$double.eps, label = "n=4 outlier")

  # n=5: sorted input (tests frac1=0 branch for certain n)
  x5 <- c(1, 2, 3, 4, 5)
  expect_equal(iqr_scaled(x5, constant = 1), stats::IQR(x5, type = 7),
               tolerance = 8 * .Machine$double.eps, label = "n=5 sorted")
})

test_that("iqr matches stats::IQR(type=7) at path boundaries n=16,17,256,257", {
  # n=16: last value in small_sort path
  # n=17: first value in micro (pdqselect) path
  # n=256: last value in micro path (IQR_INLINE_LIMIT)
  # n=257: first value in iqr_impl_large path
  for (n in c(16L, 17L, 256L, 257L)) {
    set.seed(n + 9000L)
    x <- rnorm(n)
    expected <- stats::IQR(x, type = 7)
    actual   <- iqr_scaled(x, constant = 1)
    expect_equal(actual, expected,
                 tolerance = 8 * .Machine$double.eps,
                 label = paste0("boundary n=", n))
  }
})

test_that("iqr matches stats::IQR(type=7) for integer-valued inputs", {
  for (n in c(2L, 5L, 10L, 20L, 50L)) {
    set.seed(n + 2000L)
    x <- as.double(sample.int(100L, n, replace = TRUE))
    expected <- stats::IQR(x, type = 7)
    actual   <- iqr_scaled(x, constant = 1)
    expect_equal(actual, expected,
                 tolerance = 8 * .Machine$double.eps,
                 label = paste0("integer-valued n=", n))
  }
})
