# Correctness tests for Qn after adaptive final-diff selection

test_that("qn gives non-NA positive result across full size sweep", {
  for (n in c(2, 3, 4, 5, 8, 10, 15, 16, 17, 24, 32, 50, 64, 100,
              128, 256, 512, 600, 601, 1000, 5000, 10000, 50000)) {
    set.seed(42 + n)
    x <- rnorm(n)
    val <- qn(x)
    expect_false(is.na(val), label = paste("NA at n =", n))
    expect_gt(val, 0, label = paste("non-positive at n =", n))
  }
})

test_that("pdq_qn_final_threshold exists and is in valid range", {
  cfg <- get_qnsn_config()
  expect_true("pdq_qn_final_threshold" %in% names(cfg))
  thr <- cfg$pdq_qn_final_threshold
  expect_true(is.numeric(thr))
  expect_gte(thr, 2048)
  expect_lte(thr, 500000)
})

test_that("qn parity at threshold boundary", {
  thr <- get_qnsn_config()$pdq_qn_final_threshold
  for (n in c(thr - 1L, thr, thr + 1L)) {
    set.seed(55 + n)
    x <- rnorm(n)
    val <- qn(x)
    expect_gt(val, 0, label = paste("boundary n =", n))
    expect_false(is.na(val), label = paste("boundary n =", n))
  }
})

test_that("qn known smoke value (n=5)", {
  x <- c(1, 2, 4, 8, 16)
  # From test-qnsn.R: qn(x, finite.corr = FALSE) == 3 * 2.21914446598508
  expect_equal(qn(x, finite.corr = FALSE),
               3 * 2.21914446598508,
               tolerance = sqrt(.Machine$double.eps))
})

test_that("qn is deterministic (50 reps)", {
  set.seed(99)
  x <- rnorm(5000)
  vals <- replicate(50, qn(x))
  expect_true(all(vals == vals[1]))
})

test_that("qn heap path works (n=100000)", {
  set.seed(7)
  x <- rnorm(100000)
  val <- qn(x)
  expect_false(is.na(val))
  expect_gt(val, 0)
})
