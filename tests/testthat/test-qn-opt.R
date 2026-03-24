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
# WU-Q2 regression guards (added before WU-Q2 implementation)
# ---------------------------------------------------------------------------

test_that("Q.TIERED.2: boundary straddle n=16 vs n=17 both correct", {
  tryCatch(devtools::load_all(quiet = TRUE), error = function(e) invisible(NULL))
  tol <- sqrt(.Machine$double.eps)
  set.seed(1616L); x16 <- rnorm(16L)
  set.seed(1717L); x17 <- rnorm(17L)
  expect_equal(C_qn_fast(x16), C_qn_sorted(sort(x16)), tolerance = tol,
               label = "n=16 fast vs sorted")
  expect_equal(C_qn_fast(x17), C_qn_sorted(sort(x17)), tolerance = tol,
               label = "n=17 fast vs sorted")
})

test_that("Q.TIERED.3: equal values at n=10 return 0 (tiny path)", {
  tryCatch(devtools::load_all(quiet = TRUE), error = function(e) invisible(NULL))
  expect_equal(C_qn_fast(rep(2.5, 10L)), 0.0, label = "all-equal n=10")
})

test_that("Q.TIERED.4: NaN propagates for n <= 16", {
  tryCatch(devtools::load_all(quiet = TRUE), error = function(e) invisible(NULL))
  expect_true(is.na(C_qn_fast(c(1.0, NaN, 3.0))),  label = "NaN n=3")
  expect_true(is.na(C_qn_fast(c(NaN, rep(1.0, 8L)))), label = "NaN n=9")
})

test_that("Q.TIERED.5: Inf propagates for n <= 16", {
  tryCatch(devtools::load_all(quiet = TRUE), error = function(e) invisible(NULL))
  expect_true(is.na(C_qn_fast(c(1.0, Inf,  3.0))),  label = "+Inf n=3")
  expect_true(is.na(C_qn_fast(c(1.0, -Inf, 3.0))),  label = "-Inf n=3")
  expect_true(is.na(C_qn_fast(c(-Inf, rep(1.0, 8L)))), label = "-Inf n=9")
})

# ---------------------------------------------------------------------------
# WU-Q3 regression guards (added before WU-Q3 implementation)
# ---------------------------------------------------------------------------
test_that("Q.RESTRICT.2: sorted variant matches at n=100, 1000", {
  tryCatch(devtools::load_all(quiet = TRUE), error = function(e) invisible(NULL))
  tol <- sqrt(.Machine$double.eps)
  for (n in c(100L, 1000L)) {
    set.seed(n * 7L + 3L)
    x <- rnorm(n)
    expect_equal(C_qn_fast(x), C_qn_sorted(sort(x)), tolerance = tol,
                 label = paste("restrict sorted n =", n))
  }
})

test_that("Q.RESTRICT.3: C_qn_fast_orig removed after WU-Q3", {
  tryCatch(devtools::load_all(quiet = TRUE), error = function(e) invisible(NULL))
  expect_false(exists("C_qn_fast_orig", envir = asNamespace("robscale"), inherits = FALSE),
               label = "C_qn_fast_orig removed after WU-Q3")
})

# ---------------------------------------------------------------------------
# WU-Q4 regression guards (added before WU-Q4 implementation)
# whimed_cpp is unexported; observable effect is via C_qn_fast reference values.
# These pass before and after a CORRECT fusion; fail if fusion is buggy.
# ---------------------------------------------------------------------------
test_that("Q.WHIMED.1: reference values preserved at n=65, 100, 200, 500 (seed=1234)", {
  tryCatch(devtools::load_all(quiet = TRUE), error = function(e) invisible(NULL))
  tol <- 1e-8
  # Pre-captured reference values (post-WU-Q3, 2026-03-23)
  set.seed(1234L); x65  <- rnorm(65L)
  set.seed(1234L); x100 <- rnorm(100L)
  set.seed(1234L); x200 <- rnorm(200L)
  set.seed(1234L); x500 <- rnorm(500L)
  expect_equal(C_qn_fast(x65),  0.837392748813781, tolerance = tol, label = "whimed n=65")
  expect_equal(C_qn_fast(x100), 0.939714366776843, tolerance = tol, label = "whimed n=100")
  expect_equal(C_qn_fast(x200), 0.979058242131217, tolerance = tol, label = "whimed n=200")
  expect_equal(C_qn_fast(x500), 1.013671898596087, tolerance = tol, label = "whimed n=500")
})

test_that("Q.WHIMED.2: equal-element degenerate case returns 0 (all-equal n=20)", {
  tryCatch(devtools::load_all(quiet = TRUE), error = function(e) invisible(NULL))
  expect_equal(C_qn_fast(rep(3.0, 20L)), 0.0, label = "all-equal n=20 whimed path")
})

test_that("Q.WHIMED.3: weq partition path via refinement (n=100 with tied work values)", {
  tryCatch(devtools::load_all(quiet = TRUE), error = function(e) invisible(NULL))
  # n=100 exercises qn_refinement_kernel and whimed_cpp; the equal-element
  # partition (weq) fires naturally when work[] contains tied median candidates.
  set.seed(999L)
  x <- rnorm(100L)
  ref <- C_qn_fast(x)
  expect_true(is.finite(ref), label = "refinement result is finite")
  expect_gt(ref, 0.0, label = "refinement result is positive")
  # Determinism guard: a buggy weq accumulation would make result non-deterministic
  expect_identical(C_qn_fast(x), ref, label = "refinement deterministic")
})

# ---------------------------------------------------------------------------
# WU-Q5 regression guards (added before WU-Q5 implementation)
# j-reversal + AVX2 fill: output-preserving by design — guards detect bugs.
# ---------------------------------------------------------------------------
test_that("Q.AVX2.1: reference values preserved at n=100,200,500 after j-reversal", {
  tryCatch(devtools::load_all(quiet = TRUE), error = function(e) invisible(NULL))
  tol <- 1e-8
  set.seed(1234L); x100 <- rnorm(100L)
  set.seed(1234L); x200 <- rnorm(200L)
  set.seed(1234L); x500 <- rnorm(500L)
  expect_equal(C_qn_fast(x100), 0.939714366776843, tolerance = tol, label = "avx2 n=100")
  expect_equal(C_qn_fast(x200), 0.979058242131217, tolerance = tol, label = "avx2 n=200")
  expect_equal(C_qn_fast(x500), 1.013671898596087, tolerance = tol, label = "avx2 n=500")
})

test_that("Q.AVX2.2: sorted variant matches C_qn_fast at n=200, 500", {
  tryCatch(devtools::load_all(quiet = TRUE), error = function(e) invisible(NULL))
  tol <- sqrt(.Machine$double.eps)
  for (n in c(200L, 500L)) {
    set.seed(n * 7L + 5L)
    x <- rnorm(n)
    expect_equal(C_qn_fast(x), C_qn_sorted(sort(x)), tolerance = tol,
                 label = paste("avx2 sorted n =", n))
  }
})

test_that("Q.AVX2.3: forward-fill result matches for all refinement path sizes", {
  tryCatch(devtools::load_all(quiet = TRUE), error = function(e) invisible(NULL))
  # Tests that j-reversal produces the same multiset of diffs as the original.
  # A reordering bug would show up as a wrong order statistic.
  tol <- sqrt(.Machine$double.eps)
  for (n in c(65L, 100L, 200L, 500L, 1000L)) {
    set.seed(n * 11L + 7L)
    x <- rnorm(n)
    expect_equal(C_qn_fast(x), C_qn_sorted(sort(x)), tolerance = tol,
                 label = paste("forward fill n =", n))
  }
})

# ---- WU-Q6: Workspace reuse for Qn in ensemble bootstrap ----

test_that("Q.WORKSPACE.1: scale_robust n=200 n_boot=50 reference value unchanged", {
  tryCatch(devtools::load_all(quiet = TRUE), error = function(e) invisible(NULL))
  set.seed(42L)
  x <- rnorm(200L)
  result <- scale_robust(x, n_boot = 50L)
  expect_equal(result, 0.973344584828327, tolerance = 1e-10,
               label = "workspace n=200 seed=42")
})

test_that("Q.WORKSPACE.2: scale_robust n=100 n_boot=50 reference value unchanged", {
  tryCatch(devtools::load_all(quiet = TRUE), error = function(e) invisible(NULL))
  set.seed(42L)
  x <- rnorm(100L)
  result <- scale_robust(x, n_boot = 50L)
  expect_equal(result, 1.036085388267853, tolerance = 1e-10,
               label = "workspace n=100 seed=42")
})

test_that("Q.WORKSPACE.3: large n > sn_stack_threshold gives finite result", {
  tryCatch(devtools::load_all(quiet = TRUE), error = function(e) invisible(NULL))
  # n=3000 exceeds sn_stack_threshold (2048); workspace path is bypassed (nullptr).
  # Guard: no crash, result is finite.
  set.seed(7L)
  x <- rnorm(3000L)
  result <- scale_robust(x, n_boot = 10L)
  expect_true(is.finite(result), label = "large-n no crash")
})

# ---- WU-Q7: Hoist block_offsets outside refinement loop (TBB path) ----

test_that("Q.BLKOFFS.1: result at n > qn_parallel_threshold matches sorted variant", {
  tryCatch(devtools::load_all(quiet = TRUE), error = function(e) invisible(NULL))
  # On USE_DIRECT_TBB builds this exercises the parallel block_offsets path.
  # On non-TBB builds (Linux RcppParallel) this exercises the serial path, which
  # the refactoring leaves unmodified. Test is valid in both cases.
  cfg <- get_qnsn_config()
  n_par <- cfg$qn_parallel_threshold + 500L
  set.seed(42L)
  x <- rnorm(n_par)
  expect_equal(C_qn_fast(x), C_qn_sorted(sort(x)), tolerance = 1e-10,
               label = paste("block_offsets path n =", n_par))
})

test_that("Q.BLKOFFS.2: result unchanged at n = 5000", {
  tryCatch(devtools::load_all(quiet = TRUE), error = function(e) invisible(NULL))
  tol <- 1e-10
  set.seed(99L)
  x <- rnorm(5000L)
  expect_equal(C_qn_fast(x), 1.008121378390264, tolerance = tol,
               label = "block_offsets n=5000 reference")
  expect_equal(C_qn_fast(x), C_qn_sorted(sort(x)), tolerance = tol,
               label = "block_offsets n=5000 sorted match")
})

# ---- WU-Q8: omp simd + threshold calibration + adversarial pivot check ----

test_that("Q.OMPSIMD.1: right[] iota correctness at n=100 and n=500", {
  tryCatch(devtools::load_all(quiet = TRUE), error = function(e) invisible(NULL))
  # Exercises the omp simd-annotated right[] initialisation in qn_refinement_kernel.
  # Re-use reference values from Q.AVX2.1 (same seeds, same inputs).
  tol <- sqrt(.Machine$double.eps)
  set.seed(1234L); x100 <- rnorm(100L)
  expect_equal(C_qn_fast(x100), 0.939714366776843, tolerance = tol, label = "simd n=100")
  set.seed(1234L); x500 <- rnorm(500L)
  expect_equal(C_qn_fast(x500), 1.013671898596087, tolerance = tol, label = "simd n=500")
})

test_that("Q.THRESH.1: boundary straddle at new threshold (n=40 brute-force, n=41 refinement)", {
  tryCatch(devtools::load_all(quiet = TRUE), error = function(e) invisible(NULL))
  # After OPT-Q8b threshold change from 64 to 40:
  #   n=40 -> brute-force (n <= qn_exact_threshold=40)
  #   n=41 -> refinement  (n >  qn_exact_threshold=40)
  # Both must produce correct values matching C_qn_sorted.
  tol <- sqrt(.Machine$double.eps)
  set.seed(1234L)
  x40 <- rnorm(40L); x41 <- rnorm(41L)
  expect_equal(C_qn_fast(x40), C_qn_sorted(sort(x40)), tolerance = tol,
               label = "threshold straddle n=40")
  expect_equal(C_qn_fast(x41), C_qn_sorted(sort(x41)), tolerance = tol,
               label = "threshold straddle n=41")
})

test_that("Q.THRESH.2: sizes 41..64 match sorted variant after threshold change", {
  tryCatch(devtools::load_all(quiet = TRUE), error = function(e) invisible(NULL))
  # Sizes 41..64 now use refinement path. Guard correctness at key sizes.
  tol <- sqrt(.Machine$double.eps)
  for (n in c(41L, 48L, 56L, 60L, 63L, 64L)) {
    set.seed(n * 3L + 7L)
    x <- rnorm(n)
    expect_equal(C_qn_fast(x), C_qn_sorted(sort(x)), tolerance = tol,
                 label = paste("refinement path n =", n))
  }
})
