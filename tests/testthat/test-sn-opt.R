# tests/testthat/test-sn-opt.R
#
# Correctness gate tests for OPT-S1..S8 (sn() performance optimizations).
# All tests must pass against CURRENT, UNMODIFIED code — this is a baseline
# sanity check, not a TDD red.  They remain green throughout all work units.
#
# Comparison strategy for robustbase::Sn:
#   Use finite.corr = FALSE in both sn() and robustbase::Sn.
#   Our package uses CONST_SN = 1.19259855312321 (high precision); robustbase
#   uses 1.1926 (4-decimal).  The systematic relative difference is ~1.21e-6.
#   Tolerance 1e-5 provides 8× headroom above the constant-rounding bias.

# ---------------------------------------------------------------------------
# S.1 — Known-value test: sn(c(1,2,4,8,16), finite.corr=FALSE) == 3 * CONST_SN
# ---------------------------------------------------------------------------
test_that("S.1 known value: sn(c(1,2,4,8,16), finite.corr=FALSE)", {
  tryCatch(devtools::load_all(quiet = TRUE), error = function(e) invisible(NULL))
  x <- c(1, 2, 4, 8, 16)
  expect_identical(sn(x, finite.corr = FALSE), 3 * 1.19259855312321)
})

# ---------------------------------------------------------------------------
# S.2 — Size sweep n=2..50 + {100,500,1000,5000} vs robustbase::Sn
# ---------------------------------------------------------------------------
test_that("S.2 size sweep n=2..50 + large vs robustbase::Sn", {
  tryCatch(devtools::load_all(quiet = TRUE), error = function(e) invisible(NULL))
  skip_if_not_installed("robustbase")
  tol <- 1e-5  # 8x headroom above the ~1.21e-6 CONST_SN constant-rounding bias

  ns <- c(2L:50L, 100L, 500L, 1000L, 5000L)
  for (n in ns) {
    set.seed(n * 17L + 3L)
    x <- rnorm(n)
    got <- sn(x, finite.corr = FALSE)
    ref <- robustbase::Sn(x, finite.corr = FALSE)
    expect_equal(got, ref, tolerance = tol,
                 label = paste("size sweep n =", n))
  }
})

# ---------------------------------------------------------------------------
# S.3 — Sorted variant equivalence: C_sn_sorted(sort(x)) == sn(x)
# ---------------------------------------------------------------------------
test_that("S.3 sorted variant equivalence: C_sn_sorted(sort(x)) == sn(x)", {
  tryCatch(devtools::load_all(quiet = TRUE), error = function(e) invisible(NULL))
  tol <- sqrt(.Machine$double.eps)
  for (n in c(5L, 50L, 500L, 5000L)) {
    set.seed(n * 31L)
    x <- rnorm(n)
    got <- C_sn_sorted(sort(x))
    ref <- sn(x)
    expect_equal(got, ref, tolerance = tol,
                 label = paste("C_sn_sorted equivalence n =", n))
  }
})

# ---------------------------------------------------------------------------
# S.4 — Ensemble integration: scale_robust(x, n_boot=50) finite, positive, deterministic
# ---------------------------------------------------------------------------
test_that("S.4 ensemble integration: scale_robust is finite, positive, deterministic", {
  tryCatch(devtools::load_all(quiet = TRUE), error = function(e) invisible(NULL))
  set.seed(88)
  x <- rnorm(10)
  v1 <- scale_robust(x, n_boot = 50)
  expect_true(is.finite(v1), label = "ensemble result is finite")
  expect_gt(v1, 0, label = "ensemble result is positive")
  # determinism: repeat 5 times, must be identical
  for (i in seq_len(5L)) {
    v2 <- scale_robust(x, n_boot = 50)
    expect_identical(v1, v2, label = paste("ensemble determinism rep", i))
  }
})

# ---------------------------------------------------------------------------
# S.5 — Edge cases
# ---------------------------------------------------------------------------
test_that("S.5 edge case: n=1 returns NA", {
  tryCatch(devtools::load_all(quiet = TRUE), error = function(e) invisible(NULL))
  expect_true(is.na(sn(1)))
  expect_true(is.na(sn(1.0)))
})

test_that("S.5 edge case: n=2 returns a valid positive value", {
  tryCatch(devtools::load_all(quiet = TRUE), error = function(e) invisible(NULL))
  v <- sn(c(1, 3))
  expect_true(is.finite(v))
  expect_gt(v, 0)
})

test_that("S.5 edge case: constant data returns 0", {
  tryCatch(devtools::load_all(quiet = TRUE), error = function(e) invisible(NULL))
  expect_equal(sn(rep(5.0, 10)), 0)
  expect_equal(sn(rep(5.0, 100)), 0)
})

test_that("S.5 edge case: integer input matches double input", {
  tryCatch(devtools::load_all(quiet = TRUE), error = function(e) invisible(NULL))
  x_int <- c(1L, 2L, 4L, 8L, 16L)
  x_dbl <- as.double(x_int)
  expect_identical(sn(x_int, finite.corr = FALSE),
                   sn(x_dbl, finite.corr = FALSE))
  # Also verify known value via integer path
  expect_identical(sn(x_int, finite.corr = FALSE), 3 * 1.19259855312321)
})

# ---------------------------------------------------------------------------
# S.6 — Determinism: 30 reps at n=500, all identical
# ---------------------------------------------------------------------------
test_that("S.6 determinism: 30 reps at n=500 all identical", {
  tryCatch(devtools::load_all(quiet = TRUE), error = function(e) invisible(NULL))
  set.seed(500L)
  x <- rnorm(500L)
  ref <- sn(x)
  for (i in seq_len(30L)) {
    expect_identical(sn(x), ref,
                     label = paste("determinism rep", i, "n=500"))
  }
})

# ---------------------------------------------------------------------------
# S.7 — Boundary: n=16 and n=17 (sort-network / general path boundary)
# ---------------------------------------------------------------------------
test_that("S.7 boundary: n=16 and n=17 both correct", {
  tryCatch(devtools::load_all(quiet = TRUE), error = function(e) invisible(NULL))
  skip_if_not_installed("robustbase")
  tol <- 1e-5
  set.seed(42L)
  x16 <- rnorm(16L)
  x17 <- c(x16, rnorm(1L))
  expect_equal(sn(x16, finite.corr = FALSE),
               robustbase::Sn(x16, finite.corr = FALSE),
               tolerance = tol, label = "n=16 boundary")
  expect_equal(sn(x17, finite.corr = FALSE),
               robustbase::Sn(x17, finite.corr = FALSE),
               tolerance = tol, label = "n=17 boundary")
})

# ---------------------------------------------------------------------------
# S.8 — Boundary: n=128 and n=129 (micro-buffer / stack-buffer boundary)
# ---------------------------------------------------------------------------
test_that("S.8 boundary: n=128 and n=129 both correct", {
  tryCatch(devtools::load_all(quiet = TRUE), error = function(e) invisible(NULL))
  skip_if_not_installed("robustbase")
  tol <- 1e-5
  set.seed(128L)
  x128 <- rnorm(128L)
  x129 <- c(x128, rnorm(1L))
  expect_equal(sn(x128, finite.corr = FALSE),
               robustbase::Sn(x128, finite.corr = FALSE),
               tolerance = tol, label = "n=128 boundary")
  expect_equal(sn(x129, finite.corr = FALSE),
               robustbase::Sn(x129, finite.corr = FALSE),
               tolerance = tol, label = "n=129 boundary")
})

# ---------------------------------------------------------------------------
# S.9 — Sorted variant C_sn_sorted: diagnostic export present and correct
#        Graceful skip if C_sn_fast_orig exists (OPT-S4 phase — no change to sorted)
# ---------------------------------------------------------------------------
test_that("S.9 C_sn_sorted sorted-variant: matches sn() on sorted input", {
  tryCatch(devtools::load_all(quiet = TRUE), error = function(e) invisible(NULL))
  tol <- sqrt(.Machine$double.eps)
  for (n in c(2L, 5L, 16L, 17L, 100L, 1000L)) {
    set.seed(600L + n)
    x <- sort(rnorm(n))
    expect_equal(C_sn_sorted(x), sn(x),
                 tolerance = tol,
                 label = paste("C_sn_sorted sorted n =", n))
  }
})

# ---------------------------------------------------------------------------
# S.10 — Optional: C_sn_fast_orig diagnostic export (OPT-S4)
#         Skip gracefully when not present (Phases S0..S3)
# ---------------------------------------------------------------------------
test_that("S.10 C_sn_fast_orig diagnostic export (skip if absent)", {
  tryCatch(devtools::load_all(quiet = TRUE), error = function(e) invisible(NULL))
  if (!exists("C_sn_fast_orig", mode = "function")) {
    skip("C_sn_fast_orig diagnostic export not present (added in OPT-S4)")
  }
  tol <- sqrt(.Machine$double.eps)
  for (n in c(5L, 50L, 500L)) {
    set.seed(700L + n)
    x <- rnorm(n)
    expect_equal(C_sn_fast_orig(x), sn(x),
                 tolerance = tol,
                 label = paste("C_sn_fast_orig n =", n))
  }
})
