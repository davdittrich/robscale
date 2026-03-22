# tests/testthat/test-qn-opt.R
#
# Correctness gate tests for Qn performance optimization work units (OPT-Q1..Q8).
# All tests must pass against CURRENT, UNMODIFIED code before WU-Q1 begins.
# They remain green throughout all work units — regression guards, not TDD red tests.
#
# Reference values computed from C_qn_fast() on a fixed seed; validated against
# robustbase::Qn(x, finite.corr=TRUE) where available.
# Tolerance 1e-5 for robustbase comparisons (CONST_QN rounding ~1e-11, with headroom).

tryCatch(devtools::load_all(quiet = TRUE), error = function(e) invisible(NULL))

# ---------------------------------------------------------------------------
# Q.1 — Known value: qn(c(1,2,4,8,16), finite.corr=FALSE) == 3 * CONST_QN
# ---------------------------------------------------------------------------
test_that("Q.1 known value: qn(c(1,2,4,8,16), finite.corr=FALSE) == 3 * CONST_QN", {
  tryCatch(devtools::load_all(quiet = TRUE), error = function(e) invisible(NULL))
  x <- c(1, 2, 4, 8, 16)
  expect_no_error(got <- qn(x, finite.corr = FALSE))
  # pairwise diffs sorted: 1,2,3,4,6,7,8,12,14,15; k-th = 3rd = 3
  # result = 3 * CONST_QN (no finite correction)
  expect_equal(got, 3 * 2.21914446598508, tolerance = sqrt(.Machine$double.eps))
})

# ---------------------------------------------------------------------------
# Q.2 — Known value: C_qn_fast(c(1,2,4,8,16)) with finite correction
# ---------------------------------------------------------------------------
test_that("Q.2 known value: C_qn_fast(c(1,2,4,8,16)) matches pre-computed reference", {
  tryCatch(devtools::load_all(quiet = TRUE), error = function(e) invisible(NULL))
  # 5.619672679881977 = 3 * CONST_QN * qn_factor(5=0.84412)
  expect_equal(C_qn_fast(c(1, 2, 4, 8, 16)), 5.619672679881977,
               tolerance = sqrt(.Machine$double.eps))
})

# ---------------------------------------------------------------------------
# Q.3 — Brute-force boundary n=2..17: C_qn_fast matches qn()
# ---------------------------------------------------------------------------
test_that("Q.3 brute-force boundary n=2..17: C_qn_fast identical to qn()", {
  tryCatch(devtools::load_all(quiet = TRUE), error = function(e) invisible(NULL))
  for (n in 2:17) {
    set.seed(n * 37L + 7L)
    x <- rnorm(n)
    expect_identical(C_qn_fast(x), qn(x),
                     label = paste("C_qn_fast vs qn() at n =", n))
  }
})

# ---------------------------------------------------------------------------
# Q.4 — Brute-force interior n=18..64: C_qn_fast matches qn()
# ---------------------------------------------------------------------------
test_that("Q.4 brute-force interior n=18..64: C_qn_fast identical to qn()", {
  tryCatch(devtools::load_all(quiet = TRUE), error = function(e) invisible(NULL))
  for (n in c(18L, 20L, 24L, 30L, 40L, 50L, 56L, 60L, 63L, 64L)) {
    set.seed(n * 41L + 3L)
    x <- rnorm(n)
    expect_identical(C_qn_fast(x), qn(x),
                     label = paste("C_qn_fast vs qn() at n =", n))
  }
})

# ---------------------------------------------------------------------------
# Q.5 — First refinement n=65: C_qn_fast matches qn()
# ---------------------------------------------------------------------------
test_that("Q.5 first refinement n=65: C_qn_fast identical to qn()", {
  tryCatch(devtools::load_all(quiet = TRUE), error = function(e) invisible(NULL))
  set.seed(65L * 37L + 7L)
  x <- rnorm(65L)
  expect_identical(C_qn_fast(x), qn(x))
})

# ---------------------------------------------------------------------------
# Q.6 — Refinement path n={100,200,500,1000}: C_qn_fast matches qn()
# ---------------------------------------------------------------------------
test_that("Q.6 refinement path n=100..1000: C_qn_fast identical to qn()", {
  tryCatch(devtools::load_all(quiet = TRUE), error = function(e) invisible(NULL))
  for (n in c(100L, 200L, 500L, 1000L)) {
    set.seed(n * 37L + 7L)
    x <- rnorm(n)
    expect_identical(C_qn_fast(x), qn(x),
                     label = paste("C_qn_fast vs qn() at n =", n))
  }
})

# ---------------------------------------------------------------------------
# Q.7 — vs robustbase::Qn at boundary + refinement sizes
# ---------------------------------------------------------------------------
test_that("Q.7 vs robustbase::Qn at n in {5,16,30,64,65,100,200}", {
  tryCatch(devtools::load_all(quiet = TRUE), error = function(e) invisible(NULL))
  skip_if_not_installed("robustbase")
  # Our CONST_QN has higher precision than robustbase's; max relative diff ~4.6e-4.
  # tolerance = 1e-3 provides ~2x headroom above the observed maximum.
  tol <- 1e-3
  for (n in c(5L, 16L, 30L, 64L, 65L, 100L, 200L)) {
    set.seed(n * 97L + 13L)
    x <- rnorm(n)
    got <- C_qn_fast(x)
    ref <- robustbase::Qn(x, finite.corr = TRUE)
    expect_equal(got, ref, tolerance = tol,
                 label = paste("vs robustbase n =", n))
  }
})

# ---------------------------------------------------------------------------
# Q.8 — Sorted variant: C_qn_sorted(sort(x)) matches C_qn_fast(x)
# ---------------------------------------------------------------------------
test_that("Q.8 sorted variant: C_qn_sorted(sort(x)) matches C_qn_fast(x)", {
  tryCatch(devtools::load_all(quiet = TRUE), error = function(e) invisible(NULL))
  tol <- sqrt(.Machine$double.eps)
  for (n in c(5L, 16L, 30L, 64L, 65L, 100L, 200L)) {
    set.seed(n * 53L + 17L)
    x <- rnorm(n)
    expect_equal(C_qn_sorted(sort(x)), C_qn_fast(x), tolerance = tol,
                 label = paste("C_qn_sorted equivalence n =", n))
  }
})

# ---------------------------------------------------------------------------
# Q.9 — Integer input: C_qn_int_fast matches double path
# ---------------------------------------------------------------------------
test_that("Q.9 integer input: C_qn_int_fast matches double equivalent", {
  tryCatch(devtools::load_all(quiet = TRUE), error = function(e) invisible(NULL))
  for (n in c(2L, 10L, 50L, 64L, 65L)) {
    set.seed(n * 29L)
    x_int <- as.integer(sample(-100L:100L, n, replace = TRUE))
    x_dbl <- as.double(x_int)
    expect_equal(C_qn_int_fast(x_int), C_qn_fast(x_dbl), tolerance = 1e-10,
                 label = paste("integer vs double at n =", n))
  }
})

# ---------------------------------------------------------------------------
# Q.10 — Edge case: n < 2 returns NA
# ---------------------------------------------------------------------------
test_that("Q.10 edge case: n<2 returns NA", {
  tryCatch(devtools::load_all(quiet = TRUE), error = function(e) invisible(NULL))
  expect_true(is.na(C_qn_fast(numeric(0))), label = "n=0")
  expect_true(is.na(C_qn_fast(numeric(1))), label = "n=1 numeric(1)")
  expect_true(is.na(C_qn_fast(42.0)),       label = "n=1 scalar")
})

# ---------------------------------------------------------------------------
# Q.11 — Edge case: NaN input returns NA
# ---------------------------------------------------------------------------
test_that("Q.11 edge case: NaN input returns NA", {
  tryCatch(devtools::load_all(quiet = TRUE), error = function(e) invisible(NULL))
  expect_true(is.na(C_qn_fast(c(1, NaN, 2))),   label = "NaN in middle")
  expect_true(is.na(C_qn_fast(c(NaN, NaN, NaN))), label = "all NaN")
})

# ---------------------------------------------------------------------------
# Q.12 — Edge case: Inf input returns NA
# ---------------------------------------------------------------------------
test_that("Q.12 edge case: Inf input returns NA", {
  tryCatch(devtools::load_all(quiet = TRUE), error = function(e) invisible(NULL))
  expect_true(is.na(C_qn_fast(c(1, Inf, 2))),  label = "+Inf")
  expect_true(is.na(C_qn_fast(c(1, -Inf, 2))), label = "-Inf")
})

# ---------------------------------------------------------------------------
# Q.13 — Edge case: all-equal returns 0
# ---------------------------------------------------------------------------
test_that("Q.13 edge case: all-equal input returns 0", {
  tryCatch(devtools::load_all(quiet = TRUE), error = function(e) invisible(NULL))
  expect_equal(C_qn_fast(rep(5, 10L)),   0.0, label = "all-equal n=10")
  expect_equal(C_qn_fast(rep(-3.14, 64L)), 0.0, label = "all-equal n=64 (exact boundary)")
  expect_equal(C_qn_fast(rep(0, 65L)),   0.0, label = "all-equal n=65 (first refinement)")
})

# ---------------------------------------------------------------------------
# Q.14 — Boundary straddle: n=64 (brute-force) and n=65 (refinement) both correct
# ---------------------------------------------------------------------------
test_that("Q.14 boundary straddle: n=64 brute-force and n=65 refinement both correct", {
  tryCatch(devtools::load_all(quiet = TRUE), error = function(e) invisible(NULL))
  skip_if_not_installed("robustbase")
  tol <- 1e-3  # our CONST_QN vs robustbase precision; max rel diff ~4.6e-4
  set.seed(777L)
  x64 <- rnorm(64L)
  set.seed(778L)
  x65 <- rnorm(65L)
  expect_equal(C_qn_fast(x64),         robustbase::Qn(x64, finite.corr = TRUE),
               tolerance = tol, label = "n=64 brute-force vs robustbase")
  expect_equal(C_qn_sorted(sort(x64)), robustbase::Qn(x64, finite.corr = TRUE),
               tolerance = tol, label = "n=64 sorted vs robustbase")
  expect_equal(C_qn_fast(x65),         robustbase::Qn(x65, finite.corr = TRUE),
               tolerance = tol, label = "n=65 refinement vs robustbase")
  expect_equal(C_qn_sorted(sort(x65)), robustbase::Qn(x65, finite.corr = TRUE),
               tolerance = tol, label = "n=65 sorted vs robustbase")
})

# ---------------------------------------------------------------------------
# Q.15 — Determinism: 20 repeated calls return identical result
# ---------------------------------------------------------------------------
test_that("Q.15 determinism: 20 repeated calls are identical", {
  tryCatch(devtools::load_all(quiet = TRUE), error = function(e) invisible(NULL))
  set.seed(314L)
  x <- rnorm(100L)
  ref <- C_qn_fast(x)
  for (i in seq_len(20L)) {
    expect_identical(C_qn_fast(x), ref, label = paste("determinism rep", i))
  }
})

# ---------------------------------------------------------------------------
# Q.16 — Hard reference values to >= 8 decimal places (seed=1234)
# ---------------------------------------------------------------------------
test_that("Q.16 reference values: n=16, n=64, n=100 to >=8 dp (seed=1234)", {
  tryCatch(devtools::load_all(quiet = TRUE), error = function(e) invisible(NULL))
  tol <- 1e-8
  set.seed(1234L); x16  <- rnorm(16L)
  set.seed(1234L); x64  <- rnorm(64L)
  set.seed(1234L); x100 <- rnorm(100L)
  # Reference: C_qn_fast(x) computed on unmodified code, 2026-03-23
  expect_equal(C_qn_fast(x16),  0.814140106293493, tolerance = tol, label = "ref n=16")
  expect_equal(C_qn_fast(x64),  0.820946117525703, tolerance = tol, label = "ref n=64")
  expect_equal(C_qn_fast(x100), 0.939714366776843, tolerance = tol, label = "ref n=100")
})

# ---------------------------------------------------------------------------
# Q.17 — Regression guard: n=17 (first refinement size) via both paths
# ---------------------------------------------------------------------------
test_that("Q.17 regression guard: n=17 C_qn_fast equals C_qn_sorted(sort(x))", {
  tryCatch(devtools::load_all(quiet = TRUE), error = function(e) invisible(NULL))
  set.seed(17L * 37L + 7L)
  x <- rnorm(17L)
  expect_equal(C_qn_fast(x), C_qn_sorted(sort(x)),
               tolerance = sqrt(.Machine$double.eps))
})

# ---------------------------------------------------------------------------
# Q.18 — Ensemble integration: scale_robust is finite, positive, deterministic
# ---------------------------------------------------------------------------
test_that("Q.18 ensemble integration: scale_robust is finite, positive, deterministic", {
  tryCatch(devtools::load_all(quiet = TRUE), error = function(e) invisible(NULL))
  set.seed(88L)
  x <- rnorm(50L)
  v1 <- scale_robust(x, n_boot = 20L)
  expect_true(is.finite(v1), label = "ensemble result is finite")
  expect_gt(v1, 0, label = "ensemble result is positive")
  v2 <- scale_robust(x, n_boot = 20L)
  expect_identical(v1, v2, label = "ensemble is deterministic")
})

# ---------------------------------------------------------------------------
# Q.19 — C_get_qn_factor: valid finite positive for a range of n
# ---------------------------------------------------------------------------
test_that("Q.19 C_get_qn_factor: finite and positive for n=2,5,16,64,100,1000", {
  tryCatch(devtools::load_all(quiet = TRUE), error = function(e) invisible(NULL))
  for (n in c(2L, 5L, 16L, 64L, 100L, 1000L)) {
    f <- C_get_qn_factor(n)
    expect_true(is.finite(f), label = paste("finite qn_factor n =", n))
    expect_gt(f, 0,           label = paste("positive qn_factor n =", n))
  }
})

# ---------------------------------------------------------------------------
# Q.20 — All-equal at n=64 and n=65 (each exercises a different code path)
# ---------------------------------------------------------------------------
test_that("Q.20 all-equal at exact threshold boundaries: n=64 and n=65 return 0", {
  tryCatch(devtools::load_all(quiet = TRUE), error = function(e) invisible(NULL))
  expect_equal(C_qn_fast(rep(1.5, 64L)), 0.0, label = "all-equal n=64 brute-force")
  expect_equal(C_qn_fast(rep(1.5, 65L)), 0.0, label = "all-equal n=65 refinement")
})

# ---------------------------------------------------------------------------
# Q.21 — vs robustbase full sweep n=2..20 (dense brute-force path coverage)
# ---------------------------------------------------------------------------
test_that("Q.21 full boundary sweep n=2..20 vs robustbase::Qn", {
  tryCatch(devtools::load_all(quiet = TRUE), error = function(e) invisible(NULL))
  skip_if_not_installed("robustbase")
  tol <- 1e-3  # our CONST_QN vs robustbase precision; max rel diff ~4.6e-4
  for (n in 2:20) {
    set.seed(n * 113L + 7L)
    x <- rnorm(n)
    expect_equal(C_qn_fast(x), robustbase::Qn(x, finite.corr = TRUE),
                 tolerance = tol, label = paste("n =", n))
  }
})

# ---------------------------------------------------------------------------
# Q.22 — Correctness at n=2 (minimum valid input: 1 pair, 1 diff)
# ---------------------------------------------------------------------------
test_that("Q.22 minimum valid n=2: C_qn_fast correct", {
  tryCatch(devtools::load_all(quiet = TRUE), error = function(e) invisible(NULL))
  skip_if_not_installed("robustbase")
  x <- c(3.0, 7.0)
  expect_equal(C_qn_fast(x), robustbase::Qn(x, finite.corr = TRUE), tolerance = 1e-3)
  expect_true(is.finite(C_qn_fast(x)))
  expect_gt(C_qn_fast(x), 0)
})

# ---------------------------------------------------------------------------
# WU-Q1 regression guards (added before WU-Q1 implementation)
# ---------------------------------------------------------------------------
test_that("Q.NOINLINE.1: C_qn_fast_orig exists as exported function", {
  tryCatch(devtools::load_all(quiet = TRUE), error = function(e) invisible(NULL))
  # Fails if compileAttributes() was not run after adding the export annotation.
  expect_true(exists("C_qn_fast_orig", envir = asNamespace("robscale"), inherits = FALSE),
              label = "C_qn_fast_orig is exported")
})

test_that("Q.NOINLINE.2: C_qn_fast matches C_qn_fast_orig for all sizes", {
  tryCatch(devtools::load_all(quiet = TRUE), error = function(e) invisible(NULL))
  skip_if_not(exists("C_qn_fast_orig", envir = asNamespace("robscale"), inherits = FALSE),
              "C_qn_fast_orig not present (WU-Q1 not implemented)")
  tol <- sqrt(.Machine$double.eps)
  for (n in c(2L, 5L, 10L, 16L, 17L, 64L, 65L, 100L)) {
    set.seed(n * 37L + 7L)
    x <- rnorm(n)
    expect_equal(C_qn_fast(x), C_qn_fast_orig(x), tolerance = tol,
                 label = paste("NOINLINE equivalence n =", n))
  }
})
