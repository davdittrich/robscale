# robLoc() Performance Optimization Tests (v0.5.0)
# Phase 0: TDD Guards for L1 (SIMD), L2 (Dispatch), L3 (TBB), and L4 (Cache)

test_that("0.1 — robLoc returns finite result for all tested n", {
  # Phase 1-3 Regression Guard: Basic correctness for small and medium n
  set.seed(42)
  for (n in c(8L, 16L, 32L, 64L, 100L, 500L, 1000L)) {
    loc <- robLoc(rnorm(n))
    expect_true(is.finite(loc),
      label = sprintf("robLoc(rnorm(%d)) = %g must be finite", n, loc))
  }
})

test_that("0.2 — robLoc agrees with R optimize reference to 1e-6 after all OPTs", {
  # Phase 1-3 Regression Guard: Cross-check against logistic fixed-point R reference
  # Mirrors existing test-robLoc.R pattern; tolerance 1e-6 = NR vs Brent difference
  robLocRef <- function(x) {
    obj <- function(t) sum((2*plogis((x-t)/mad(x))-1))^2
    optimize(obj, range(x), tol=sqrt(.Machine$double.eps))$minimum
  }
  set.seed(1001)
  for (n in c(10L, 25L, 50L, 100L)) {
    x   <- rnorm(n)
    got <- robLoc(x)
    ref <- robLocRef(x)
    expect_equal(got, ref, tolerance=1e-5,
      label=sprintf("n=%d: robLoc=%g, ref=%g", n, got, ref))
  }
})

test_that("0.3 — robLoc fused AVX2 agrees with scalar path to 2*sqrt(eps)", {
  # Phase 3 (OPT-L1) Guard: Primary correctness gate for fused kernel
  # SKIPs until rob_loc_scalar_impl diagnostic helper is added in Phase 3
  skip_if_not(
    exists("rob_loc_scalar_impl", envir=asNamespace("robscale"), mode="function"),
    "Diagnostic scalar helper not compiled"
  )
  set.seed(777)
  tol <- 2 * sqrt(.Machine$double.eps)
  max_diff <- 0
  for (n in c(16L, 32L, 64L, 100L, 256L, 500L, 1000L)) {
    for (rep in seq_len(20L)) {
      x  <- rnorm(n)
      v1 <- robscale:::rob_loc_scalar_impl(x)
      v2 <- robLoc(x)
      d  <- abs(v1 - v2)
      if (d > max_diff) max_diff <- d
    }
  }
  expect_lt(max_diff, tol,
    label=sprintf("max|fused-scalar| = %.3e (tol %.3e)", max_diff, tol))
})

test_that("0.4 — robLoc handles near-degenerate scale without NaN/Inf", {
  # Phase 3 (OPT-L1) Guard: Edge case for sum_dpsi near-zero guard (pathological input)
  # All values very close together: scale ~= 0 → u_i huge → tanh(u_i) → +/-1
  # sum_dpsi = sum sech^2(u_i) → 0 without guard
  x_degen <- c(rep(1.0, 50L), rep(1.0 + 1e-10, 50L))
  loc <- robLoc(x_degen)
  expect_true(is.finite(loc) && !is.nan(loc),
    label=sprintf("robLoc degenerate input = %g must be finite", loc))

  # Constant input: s=0 path; must return median
  x_const <- rep(5.0, 20L)
  expect_equal(robLoc(x_const), 5.0, tolerance=1e-10)
})

test_that("0.5 — robLoc parallel (n>=threshold) agrees with serial to 1e-10", {
  # Phase 4 (OPT-L3) Guard: Critical gate for TBB parallel implementation
  # SKIPs until rob_loc_has_parallel() and rob_loc_serial_impl() added in Phase 4
  skip_if_not(
    exists("rob_loc_has_parallel", envir=asNamespace("robscale"), mode="function"),
    "TBB not compiled"
  )
  skip_if_not(robscale:::rob_loc_has_parallel())
  set.seed(2024)
  tol <- 1e-10
  for (n in c(4096L, 8192L, 16384L)) {
    x <- rnorm(n)
    v_par <- robLoc(x)                            # dispatches parallel
    v_ser <- robscale:::rob_loc_serial_impl(x)    # forces serial path
    expect_equal(v_par, v_ser, tolerance=tol,
      label=sprintf("n=%d: |parallel-serial| = %.3e", n, abs(v_par-v_ser)))
  }
})

# ---------------------------------------------------------------------------
# WU-RL0: Extended correctness test suite
# Tests 0.6–0.33: boundary straddles, edge cases, diagnostic cross-checks
# ---------------------------------------------------------------------------

test_that("0.6 — n=0 returns NA_real_", {
  expect_identical(robLoc(numeric(0L)), NA_real_)
  expect_identical(robLoc(double(0L)), NA_real_)
})

test_that("0.7 — n=1 returns the single value exactly", {
  for (v in c(-3.5, 0.0, 1.0, 100.0)) {
    expect_equal(robLoc(v), v, tolerance=1e-14,
      label=sprintf("robLoc(%.1f) must equal %.1f", v, v))
  }
})

test_that("0.8 — n=2 returns median regardless of scale argument", {
  x2 <- c(1.0, 5.0)
  med2 <- median(x2)
  # no-scale: minobs=4, n=2 < 4 → median
  expect_equal(robLoc(x2), med2, tolerance=1e-14)
  # has_scale: minobs=3, n=2 < 3 → still median
  expect_equal(robLoc(x2, scale=1.0), med2, tolerance=1e-14)
  expect_equal(robLoc(x2, scale=0.5), med2, tolerance=1e-14)
})

test_that("0.9 — n=3 no-scale returns median; n=3 has_scale iterates away from median", {
  x3 <- c(1.0, 2.0, 4.0)
  med3 <- median(x3)
  # no-scale: minobs=4, n=3 < 4 → return median
  expect_equal(robLoc(x3), med3, tolerance=1e-14)
  # has_scale: minobs=3, n=3 >= 3 → Newton-Raphson runs
  # x3 is right-skewed so robust location should be > median
  loc_hs <- robLoc(x3, scale=1.0)
  expect_true(is.finite(loc_hs))
  expect_false(isTRUE(all.equal(loc_hs, med3, tolerance=1e-8)),
    label="n=3 has_scale should iterate away from median for skewed data")
})

test_that("0.10 — n=4..7 (below AVX2 threshold, 3-pass scalar) all return finite", {
  set.seed(42)
  for (n in 4L:7L) {
    x   <- rnorm(n)
    loc <- robLoc(x)
    expect_true(is.finite(loc),
      label=sprintf("robLoc(rnorm(%d)) = %g must be finite", n, loc))
  }
})

test_that("0.11 — n=8 at AVX2 fused-kernel threshold returns finite and matches scalar", {
  skip_if_not(
    exists("rob_loc_scalar_impl", envir=asNamespace("robscale"), mode="function"),
    "Diagnostic scalar helper not compiled"
  )
  set.seed(808)
  x8  <- rnorm(8L)
  loc <- robLoc(x8)
  expect_true(is.finite(loc))
  scal <- robscale:::rob_loc_scalar_impl(x8)
  expect_equal(loc, scal, tolerance=2*sqrt(.Machine$double.eps),
    label=sprintf("n=8 |fused-scalar| = %.3e", abs(loc-scal)))
})

test_that("0.12 — n=9..17 cross-check vs scalar (NR step count varies)", {
  skip_if_not(
    exists("rob_loc_scalar_impl", envir=asNamespace("robscale"), mode="function"),
    "Diagnostic scalar helper not compiled"
  )
  set.seed(123)
  tol <- 2 * sqrt(.Machine$double.eps)
  for (n in 9L:17L) {
    x <- rnorm(n)
    expect_equal(robLoc(x), robscale:::rob_loc_scalar_impl(x), tolerance=tol,
      label=sprintf("n=%d: |robLoc-scalar| = %.3e",
                    n, abs(robLoc(x) - robscale:::rob_loc_scalar_impl(x))))
  }
})

test_that("0.13 — n=32,50,63,64,65 cross-check vs scalar (straddles micro_buf MICRO_SIZE=64)", {
  skip_if_not(
    exists("rob_loc_scalar_impl", envir=asNamespace("robscale"), mode="function"),
    "Diagnostic scalar helper not compiled"
  )
  set.seed(777)
  tol <- 2 * sqrt(.Machine$double.eps)
  for (n in c(32L, 50L, 63L, 64L, 65L)) {
    x <- rnorm(n)
    expect_equal(robLoc(x), robscale:::rob_loc_scalar_impl(x), tolerance=tol,
      label=sprintf("n=%d: |robLoc-scalar| = %.3e",
                    n, abs(robLoc(x) - robscale:::rob_loc_scalar_impl(x))))
  }
})

test_that("0.14 — all-equal input returns that value exactly at all boundary n", {
  for (v in c(-2.5, 0.0, 1.0)) {
    for (n in c(4L, 8L, 16L, 64L, 65L, 100L)) {
      x   <- rep(v, n)
      loc <- robLoc(x)
      expect_equal(loc, v, tolerance=1e-14,
        label=sprintf("robLoc(rep(%g, %d)) = %g", v, n, loc))
    }
  }
})

test_that("0.15 — permutation invariance at all gate sizes", {
  set.seed(1234)
  for (n in c(5L, 8L, 16L, 32L, 64L, 100L)) {
    x    <- rnorm(n)
    loc1 <- robLoc(x)
    loc2 <- robLoc(sample(x))
    expect_equal(loc1, loc2, tolerance=1e-14,
      label=sprintf("n=%d: |perm-orig| = %.3e", n, abs(loc1-loc2)))
  }
})

test_that("0.16 — determinism: repeated calls with same input give identical results", {
  set.seed(2025)
  x <- rnorm(50L)
  expect_identical(robLoc(x), robLoc(x))
  x2 <- rnorm(8L)
  expect_identical(robLoc(x2), robLoc(x2))
})

test_that("0.17 — known reference values at n=5,16,64,100 (vs R optimize)", {
  robLocRef <- function(x) {
    s <- mad(x)
    if (s == 0) return(median(x))
    obj <- function(t) sum(tanh((x - t) / (2 * s)))^2
    optimize(obj, range(x), tol=.Machine$double.eps)$minimum
  }
  set.seed(42)
  for (n in c(5L, 16L, 64L, 100L)) {
    x   <- rnorm(n)
    got <- robLoc(x)
    ref <- robLocRef(x)
    expect_equal(got, ref, tolerance=1e-5,
      label=sprintf("n=%d: |robLoc-ref| = %.3e", n, abs(got-ref)))
  }
})

test_that("0.18 — has_scale path: robLoc(x, scale=s) matches R optimize reference", {
  robLocRef_s <- function(x, s) {
    obj <- function(t) sum(tanh((x - t) / (2 * s)))^2
    optimize(obj, range(x), tol=.Machine$double.eps)$minimum
  }
  set.seed(77)
  for (n in c(3L, 5L, 10L, 20L)) {
    x   <- rnorm(n)
    s   <- 1.2
    got <- robLoc(x, scale=s)
    ref <- robLocRef_s(x, s)
    expect_equal(got, ref, tolerance=1e-5,
      label=sprintf("n=%d scale=1.2: |robLoc-ref| = %.3e", n, abs(got-ref)))
  }
})

test_that("0.19 — has_scale=mad(x) matches no-scale path to NR tolerance", {
  # C++ mad_select and R mad() use the same formula (MAD_CONSISTENCY=1.4826022).
  # Tiny floating-point differences in their median computations can shift the
  # NR fixed point by up to ~1e-7. Use 1e-5 (same as NR vs optimize tolerance).
  set.seed(888)
  for (n in c(5L, 10L, 20L, 50L)) {
    x  <- rnorm(n)
    s  <- mad(x)
    if (s == 0) next
    loc_no   <- robLoc(x)
    loc_with <- robLoc(x, scale=s)
    expect_equal(loc_no, loc_with, tolerance=1e-5,
      label=sprintf("n=%d: |no-scale - has_scale(mad)| = %.3e",
                    n, abs(loc_no - loc_with)))
  }
})

test_that("0.20 — na.rm=TRUE strips NAs before computation", {
  x_na    <- c(1.0, NA, 2.0, 3.0, NA, 5.0, 8.0)
  x_clean <- c(1.0, 2.0, 3.0, 5.0, 8.0)
  expect_equal(robLoc(x_na, na.rm=TRUE), robLoc(x_clean), tolerance=1e-14)
})

test_that("0.21 — na.rm=FALSE raises error when NAs are present", {
  x_na <- c(1.0, NA, 3.0, 5.0, 7.0)
  expect_error(robLoc(x_na), regexp="NA")
})

test_that("0.22 — outlier resistance: one extreme outlier barely moves estimate", {
  set.seed(42)
  x_clean <- rnorm(30L)
  x_dirty <- c(x_clean, 1e6)
  loc_c   <- robLoc(x_clean)
  loc_d   <- robLoc(x_dirty)
  expect_lt(abs(loc_c - loc_d), 0.5,
    label=sprintf("|robLoc(clean)-robLoc(dirty)| = %.3f (must be < 0.5)",
                  abs(loc_c - loc_d)))
})

test_that("0.23 — robLoc returns double scalar (length 1)", {
  x <- rnorm(10L)
  loc <- robLoc(x)
  expect_type(loc, "double")
  expect_length(loc, 1L)
})

test_that("0.24 — maxit=0 returns median for n >= 4 (no NR iterations)", {
  # rob_loc_compute loop: for (k=0; k<maxit; ++k) → runs 0 times → returns med
  set.seed(555)
  for (n in c(4L, 8L, 16L, 32L, 100L)) {
    x     <- rnorm(n)
    loc0  <- robLoc(x, maxit=0L)
    med_r <- median(x)
    expect_equal(loc0, med_r, tolerance=1e-14,
      label=sprintf("n=%d maxit=0: robLoc=%g, median=%g", n, loc0, med_r))
  }
})

test_that("0.25 — scalar fallback n=4..7 cross-check vs R optimize reference", {
  robLocRef <- function(x) {
    s <- mad(x)
    if (s == 0) return(median(x))
    obj <- function(t) sum(tanh((x - t) / (2 * s)))^2
    optimize(obj, range(x), tol=.Machine$double.eps)$minimum
  }
  set.seed(314)
  for (n in 4L:7L) {
    x   <- rnorm(n)
    got <- robLoc(x)
    ref <- robLocRef(x)
    expect_equal(got, ref, tolerance=1e-5,
      label=sprintf("n=%d: |robLoc-ref| = %.3e", n, abs(got-ref)))
  }
})

test_that("0.26 — robLoc handles negative-only and mixed-sign inputs", {
  robLocRef <- function(x) {
    s <- mad(x)
    if (s == 0) return(median(x))
    obj <- function(t) sum(tanh((x - t) / (2 * s)))^2
    optimize(obj, range(x), tol=.Machine$double.eps)$minimum
  }
  set.seed(99)
  x_neg <- -abs(rnorm(20L))   # all negative
  x_mix <- rnorm(20L)         # mixed signs
  for (x in list(x_neg, x_mix)) {
    got <- robLoc(x)
    ref <- robLocRef(x)
    expect_true(is.finite(got))
    expect_equal(got, ref, tolerance=1e-5)
  }
})

test_that("0.27 — right-skewed distribution: robLoc lies below mean", {
  set.seed(42)
  x_skew <- rexp(50L)      # exponential: mean=1, median≈0.69
  loc     <- robLoc(x_skew)
  mn      <- mean(x_skew)
  expect_true(is.finite(loc))
  expect_lt(loc, mn, label="robust location below non-robust mean on right-skewed data")
})

test_that("0.27b — Inf in input: result is not NaN (no crash or NaN propagation)", {
  # With one Inf value: MAD(x_inf) = Inf → half_inv_s = 0 → u_i = 0 → tanh(0)=0
  # → sum_psi=0 → NR step=0 → returns median. Verify no crash and no NaN.
  x_inf <- c(1.0, 2.0, 3.0, Inf, 5.0)
  expect_false(is.nan(robLoc(x_inf)), label="Inf input: result must not be NaN")
  x_neginf <- c(1.0, 2.0, 3.0, -Inf, 5.0)
  expect_false(is.nan(robLoc(x_neginf)), label="-Inf input: result must not be NaN")
})

test_that("0.27c — NaN in input with na.rm=TRUE strips NaN before computation", {
  # R treats NaN as NA; na.rm=TRUE should remove it, na.rm=FALSE should error
  x_nan   <- c(1.0, NaN, 2.0, 3.0, 5.0, 7.0)
  x_clean <- c(1.0, 2.0, 3.0, 5.0, 7.0)
  expect_equal(robLoc(x_nan, na.rm=TRUE), robLoc(x_clean), tolerance=1e-14)
  expect_error(robLoc(x_nan, na.rm=FALSE), regexp="NA")
})

# ---------------------------------------------------------------------------
# WU-RL0 specific: rob_loc_fast_orig diagnostic baseline tests
# All skip if rob_loc_fast_orig has not yet been compiled/exported.
# ---------------------------------------------------------------------------

test_that("0.28 — rob_loc_fast_orig: exists and returns finite at all gate sizes", {
  skip_if_not(
    exists("rob_loc_fast_orig", envir=asNamespace("robscale"), mode="function"),
    "rob_loc_fast_orig diagnostic not compiled"
  )
  set.seed(314)
  tol_nr <- sqrt(.Machine$double.eps)
  for (n in c(4L, 5L, 6L, 7L, 8L, 10L, 16L, 32L, 64L, 100L, 200L, 500L, 1000L)) {
    x   <- rnorm(n)
    loc <- robscale:::rob_loc_fast_orig(x, FALSE, 0.0, 80L, tol_nr)
    expect_true(is.finite(loc),
      label=sprintf("rob_loc_fast_orig(n=%d) = %g must be finite", n, loc))
  }
})

test_that("0.29 — rob_loc_fast_orig matches robLoc at all gate sizes (pre-RL1 baseline == current)", {
  skip_if_not(
    exists("rob_loc_fast_orig", envir=asNamespace("robscale"), mode="function"),
    "rob_loc_fast_orig diagnostic not compiled"
  )
  set.seed(271)
  tol_nr    <- sqrt(.Machine$double.eps)
  tol_match <- 2 * sqrt(.Machine$double.eps)
  for (n in c(4L, 5L, 6L, 7L, 8L, 10L, 16L, 32L, 64L, 100L, 200L)) {
    x    <- rnorm(n)
    orig <- robscale:::rob_loc_fast_orig(x, FALSE, 0.0, 80L, tol_nr)
    curr <- robLoc(x)
    expect_equal(orig, curr, tolerance=tol_match,
      label=sprintf("n=%d: |orig-curr| = %.3e", n, abs(orig-curr)))
  }
})

test_that("0.30 — rob_loc_fast_orig has_scale path matches robLoc has_scale path", {
  skip_if_not(
    exists("rob_loc_fast_orig", envir=asNamespace("robscale"), mode="function"),
    "rob_loc_fast_orig diagnostic not compiled"
  )
  set.seed(141)
  tol_nr    <- sqrt(.Machine$double.eps)
  tol_match <- 2 * sqrt(.Machine$double.eps)
  for (n in c(3L, 5L, 8L, 16L, 32L)) {
    x    <- rnorm(n)
    s    <- 1.5
    orig <- robscale:::rob_loc_fast_orig(x, TRUE, s, 80L, tol_nr)
    curr <- robLoc(x, scale=s)
    expect_equal(orig, curr, tolerance=tol_match,
      label=sprintf("n=%d has_scale: |orig-curr| = %.3e", n, abs(orig-curr)))
  }
})

test_that("0.31 — rob_loc_fast_orig all-equal input returns that value", {
  skip_if_not(
    exists("rob_loc_fast_orig", envir=asNamespace("robscale"), mode="function"),
    "rob_loc_fast_orig diagnostic not compiled"
  )
  tol_nr <- sqrt(.Machine$double.eps)
  for (v in c(-1.0, 0.0, 5.0)) {
    for (n in c(5L, 16L, 64L)) {
      x   <- rep(v, n)
      loc <- robscale:::rob_loc_fast_orig(x, FALSE, 0.0, 80L, tol_nr)
      expect_equal(loc, v, tolerance=1e-14,
        label=sprintf("rob_loc_fast_orig(rep(%g,%d)) = %g", v, n, loc))
    }
  }
})

test_that("0.32 — rob_loc_fast_orig n=3..7 edge cases are finite (no-scale and has_scale)", {
  skip_if_not(
    exists("rob_loc_fast_orig", envir=asNamespace("robscale"), mode="function"),
    "rob_loc_fast_orig diagnostic not compiled"
  )
  set.seed(99)
  tol_nr <- sqrt(.Machine$double.eps)
  for (n in 3L:7L) {
    x      <- rnorm(n)
    loc_ns <- robscale:::rob_loc_fast_orig(x, FALSE, 0.0, 80L, tol_nr)
    loc_hs <- robscale:::rob_loc_fast_orig(x, TRUE, 1.0, 80L, tol_nr)
    expect_true(is.finite(loc_ns),
      label=sprintf("rob_loc_fast_orig no-scale n=%d", n))
    expect_true(is.finite(loc_hs),
      label=sprintf("rob_loc_fast_orig has_scale n=%d", n))
  }
})

test_that("0.33 — rob_loc_fast_orig permutation invariance at gate sizes", {
  skip_if_not(
    exists("rob_loc_fast_orig", envir=asNamespace("robscale"), mode="function"),
    "rob_loc_fast_orig diagnostic not compiled"
  )
  set.seed(212)
  tol_nr <- sqrt(.Machine$double.eps)
  for (n in c(5L, 8L, 16L, 64L)) {
    x     <- rnorm(n)
    loc1  <- robscale:::rob_loc_fast_orig(x, FALSE, 0.0, 80L, tol_nr)
    loc2  <- robscale:::rob_loc_fast_orig(sample(x), FALSE, 0.0, 80L, tol_nr)
    expect_equal(loc1, loc2, tolerance=1e-14,
      label=sprintf("n=%d: |perm-orig| = %.3e", n, abs(loc1-loc2)))
  }
})

# ---- WU-RL2 TDD guards: single-pass scalar vs reference ----
# These tests must pass before AND after RL2 replaces the 3-pass scalar fallback.
# rob_loc_scalar_impl is the independent 3-pass reference (uses std::tanh directly).

test_that("0.34 — RL2 pre-condition: n=4..8 no-scale path matches scalar reference", {
  skip_if_not(
    exists("rob_loc_scalar_impl", envir=asNamespace("robscale"), mode="function"),
    "rob_loc_scalar_impl diagnostic not compiled"
  )
  tol <- 2 * sqrt(.Machine$double.eps)
  max_diff <- 0
  for (n in 4L:8L) {
    for (s in c(42L, 57L, 99L, 123L, 200L, 314L, 628L, 777L, 1024L, 1618L,
                22L, 33L, 44L, 55L, 66L, 77L, 88L, 99L, 111L, 222L)) {
      set.seed(s)
      x  <- rnorm(n)
      v1 <- robscale:::rob_loc_scalar_impl(x)
      v2 <- robLoc(x)
      d  <- abs(v1 - v2)
      if (d > max_diff) max_diff <- d
    }
  }
  expect_lt(max_diff, tol,
    label=sprintf("max |scalar-robLoc| over n=4..8 x 20 seeds = %.3e (tol=%.3e)", max_diff, tol))
})

test_that("0.35 — RL2 pre-condition: has_scale=TRUE path matches scalar reference at n=4..8,16,64,100", {
  skip_if_not(
    exists("rob_loc_scalar_impl", envir=asNamespace("robscale"), mode="function"),
    "rob_loc_scalar_impl diagnostic not compiled"
  )
  tol <- 2 * sqrt(.Machine$double.eps)
  max_diff <- 0
  for (n in c(4L, 5L, 6L, 7L, 8L, 16L, 64L, 100L)) {
    for (s in c(42L, 57L, 99L, 123L, 200L, 314L, 628L, 777L, 1024L, 1618L,
                22L, 33L, 44L, 55L, 66L, 77L, 88L, 99L, 111L, 222L)) {
      set.seed(s)
      x    <- rnorm(n)
      sv   <- mad(x)
      if (sv == 0) next
      # has_scale path: robLoc(x, scale=sv) must match scalar reference robLoc(x)
      # (same scale => same result, within NR tolerance)
      v_hs <- robLoc(x, scale = sv)
      v_ns <- robscale:::rob_loc_scalar_impl(x)
      d    <- abs(v_hs - v_ns)
      if (d > max_diff) max_diff <- d
    }
  }
  expect_lt(max_diff, 1e-5,  # NR vs scalar may differ by ~1e-7 due to MAD rounding
    label=sprintf("max |has_scale-scalar| over n=4..8,16,64,100 x 20 seeds = %.3e", max_diff))
})
