# Tests for C_sn_impl_sorted / C_qn_impl_sorted
# These call thin R-exported wrappers (always available).

# --- Equivalence: C_sn_sorted(sort(x)) == sn(x) ---

test_that("sn_sorted matches sn for various n", {

  n_grid <- c(2, 3, 5, 10, 16, 17, 32, 64, 65, 100, 128, 256, 512)
  tol <- sqrt(.Machine$double.eps)
  for (nn in n_grid) {
    set.seed(42 + nn)
    x <- rnorm(nn)
    expected <- sn(x)
    got <- robscale:::C_sn_sorted(sort(x))
    expect_equal(got, expected, tolerance = tol,
                 label = paste("sn_sorted equivalence n =", nn))
  }
})

# --- Equivalence: C_qn_sorted(sort(x)) == qn(x) ---

test_that("qn_sorted matches qn for various n", {

  n_grid <- c(2, 3, 5, 10, 16, 17, 32, 64, 65, 100, 128, 256, 512)
  tol <- sqrt(.Machine$double.eps)
  for (nn in n_grid) {
    set.seed(42 + nn)
    x <- rnorm(nn)
    expected <- qn(x)
    got <- robscale:::C_qn_sorted(sort(x))
    expect_equal(got, expected, tolerance = tol,
                 label = paste("qn_sorted equivalence n =", nn))
  }
})

# --- Known values ---

test_that("sn_sorted returns correct known value", {

  x <- c(1, 2, 4, 8, 16) # already sorted
  expect_equal(robscale:::C_sn_sorted(x), sn(x, finite.corr = TRUE))
})

test_that("qn_sorted returns correct known value", {

  x <- c(1, 2, 4, 8, 16)
  expect_equal(robscale:::C_qn_sorted(x), qn(x, finite.corr = TRUE))
})

# --- Edge cases ---

test_that("sn_sorted edge cases", {

  # n=1 -> NA
  expect_true(is.na(robscale:::C_sn_sorted(1)))
  # n=2 -> valid positive
  val <- robscale:::C_sn_sorted(c(1, 2))
  expect_true(is.finite(val) && val > 0)
  # constant data -> 0
  expect_equal(robscale:::C_sn_sorted(c(5, 5, 5, 5, 5)), 0)
})

test_that("qn_sorted edge cases", {

  # n=1 -> NA
  expect_true(is.na(robscale:::C_qn_sorted(1)))
  # n=2 -> valid positive
  val <- robscale:::C_qn_sorted(c(1, 2))
  expect_true(is.finite(val) && val > 0)
  # constant data -> 0
  expect_equal(robscale:::C_qn_sorted(c(5, 5, 5, 5, 5)), 0)
})

# --- Integer input ---

test_that("sn_sorted handles integer input (as double)", {

  x_int <- c(1L, 2L, 4L, 8L, 16L)
  expect_equal(
    robscale:::C_sn_sorted(sort(as.double(x_int))),
    sn(x_int, finite.corr = TRUE)
  )
})

test_that("qn_sorted handles integer input (as double)", {

  x_int <- c(1L, 2L, 4L, 8L, 16L)
  expect_equal(
    robscale:::C_qn_sorted(sort(as.double(x_int))),
    qn(x_int, finite.corr = TRUE)
  )
})

# --- Determinism ---

test_that("sn_sorted is deterministic (50 repeats, n=100)", {

  set.seed(42)
  x <- sort(rnorm(100))
  ref <- robscale:::C_sn_sorted(x)
  for (i in seq_len(50)) {
    expect_identical(robscale:::C_sn_sorted(x), ref)
  }
})

test_that("qn_sorted is deterministic (50 repeats, n=100)", {

  set.seed(42)
  x <- sort(rnorm(100))
  ref <- robscale:::C_qn_sorted(x)
  for (i in seq_len(50)) {
    expect_identical(robscale:::C_qn_sorted(x), ref)
  }
})
