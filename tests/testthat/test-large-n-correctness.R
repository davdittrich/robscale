# test-large-n-correctness.R
#
# Estimator correctness at large n, exercising the adaptive pdqselect dispatch.
# Each estimator is tested above its selection threshold to confirm the
# pdqselect path produces correct, deterministic results.

library(robscale)
has_robustbase <- requireNamespace("robustbase", quietly = TRUE)

# ---------------------------------------------------------------------------
# robScale
# ---------------------------------------------------------------------------

test_that("robScale returns finite positive value at large n", {
  for (n in c(100000L, 500000L, 1000000L)) {
    set.seed(42)
    x <- rnorm(n)
    val <- robscale::robScale(x)
    expect_true(is.finite(val) && val > 0,
                label = sprintf("robScale n=%d", n))
  }
})

test_that("robScale is deterministic at large n", {
  n <- 200000L
  set.seed(42)
  x <- rnorm(n)
  v1 <- robscale::robScale(x)
  v2 <- robscale::robScale(x)
  expect_equal(v1, v2)
  expect_true(v1 > 0.9 && v1 < 1.1)
})

# ---------------------------------------------------------------------------
# sn
# ---------------------------------------------------------------------------

test_that("sn returns finite positive value at large n", {
  for (n in c(50000L, 100000L, 500000L)) {
    set.seed(42)
    x <- rnorm(n)
    val <- robscale::sn(x)
    expect_true(is.finite(val) && val > 0,
                label = sprintf("sn n=%d", n))
  }
})

test_that("sn is deterministic at large n", {
  n <- 100000L
  set.seed(42); x <- rnorm(n)
  v1 <- robscale::sn(x)
  v2 <- robscale::sn(x)
  expect_equal(v1, v2)
})

test_that("sn matches robustbase::Sn at large n", {
  skip_if_not(has_robustbase, "robustbase not installed")
  n <- 50000L
  set.seed(42); x <- rnorm(n)
  val_rs  <- robscale::sn(x)
  val_rb  <- robustbase::Sn(x)
  expect_equal(val_rs, val_rb, tolerance = 1e-6,
               label = "sn vs robustbase::Sn n=50000")
})

# ---------------------------------------------------------------------------
# qn
# ---------------------------------------------------------------------------

test_that("qn returns finite positive value at large n", {
  for (n in c(10000L, 50000L, 100000L, 500000L)) {
    set.seed(42)
    x <- rnorm(n)
    val <- robscale::qn(x)
    expect_true(is.finite(val) && val > 0,
                label = sprintf("qn n=%d", n))
  }
})

test_that("qn is deterministic at large n (above parallel threshold)", {
  n <- 100000L
  set.seed(42); x <- rnorm(n)
  v1 <- robscale::qn(x)
  v2 <- robscale::qn(x)
  expect_equal(v1, v2)
})

test_that("qn matches robustbase::Qn at large n", {
  skip_if_not(has_robustbase, "robustbase not installed")
  for (n in c(10000L, 50000L)) {
    set.seed(42); x <- rnorm(n)
    val_rs  <- robscale::qn(x)
    val_rb  <- robustbase::Qn(x)
    # ~2e-6 relative difference expected: robscale uses float intermediates
    # in whimed, robustbase uses double throughout.
    expect_equal(val_rs, val_rb, tolerance = 1e-5,
                 label = sprintf("qn vs robustbase::Qn n=%d", n))
  }
})
