## test-adm-opt.R
## Correctness and regression guards for adm() optimization session.
## WU-ADM0: diagnostic export presence and baseline equivalence.
## WU-ADM1: adm_core kernel correctness (regression guards post-AVX2).
## WU-ADM2: caller-side guards (n=1, buffer boundaries, adaptive median).
## WU-ADM3: adm_core_sorted correctness and ensemble regression.
## WU-ADM4: TBB parallel path correctness, determinism, cleanup guards.

library(robscale)
tol_n <- function(n) n * .Machine$double.eps

# ============================================================
# WU-ADM0: diagnostic export presence + baseline equivalence
# ============================================================

test_that("adm opt: C_adm_orig exists as exported diagnostic function", {
  expect_true(exists("C_adm_orig",
                     envir = asNamespace("robscale"), inherits = FALSE))
})

test_that("adm opt: C_adm_fast exists as exported diagnostic function", {
  expect_true(exists("C_adm_fast",
                     envir = asNamespace("robscale"), inherits = FALSE))
})

test_that("adm opt: C_adm_orig matches adm() at n=8", {
  set.seed(8); x <- rnorm(8)
  expect_equal(robscale:::C_adm_orig(x), adm(x), tolerance = tol_n(8))
})

test_that("adm opt: C_adm_orig matches adm() at n=16", {
  set.seed(16); x <- rnorm(16)
  expect_equal(robscale:::C_adm_orig(x), adm(x), tolerance = tol_n(16))
})

test_that("adm opt: C_adm_orig matches adm() at n=64", {
  set.seed(64); x <- rnorm(64)
  expect_equal(robscale:::C_adm_orig(x), adm(x), tolerance = tol_n(64))
})

test_that("adm opt: C_adm_orig matches adm() at n=128", {
  set.seed(128); x <- rnorm(128)
  expect_equal(robscale:::C_adm_orig(x), adm(x), tolerance = tol_n(128))
})

test_that("adm opt: C_adm_orig matches adm() at n=256", {
  set.seed(256); x <- rnorm(256)
  expect_equal(robscale:::C_adm_orig(x), adm(x), tolerance = tol_n(256))
})

test_that("adm opt: C_adm_orig matches adm() at n=1024", {
  set.seed(1024); x <- rnorm(1024)
  expect_equal(robscale:::C_adm_orig(x), adm(x), tolerance = tol_n(1024))
})

test_that("adm opt: C_adm_fast matches adm() at n=8", {
  set.seed(8); x <- rnorm(8)
  expect_equal(robscale:::C_adm_fast(x), adm(x), tolerance = tol_n(8))
})

test_that("adm opt: C_adm_fast matches adm() at n=16", {
  set.seed(16); x <- rnorm(16)
  expect_equal(robscale:::C_adm_fast(x), adm(x), tolerance = tol_n(16))
})

test_that("adm opt: C_adm_fast matches adm() at n=64", {
  set.seed(64); x <- rnorm(64)
  expect_equal(robscale:::C_adm_fast(x), adm(x), tolerance = tol_n(64))
})

test_that("adm opt: C_adm_fast matches adm() at n=128", {
  set.seed(128); x <- rnorm(128)
  expect_equal(robscale:::C_adm_fast(x), adm(x), tolerance = tol_n(128))
})

test_that("adm opt: C_adm_fast matches adm() at n=256", {
  set.seed(256); x <- rnorm(256)
  expect_equal(robscale:::C_adm_fast(x), adm(x), tolerance = tol_n(256))
})

test_that("adm opt: C_adm_fast matches adm() at n=1024", {
  set.seed(1024); x <- rnorm(1024)
  expect_equal(robscale:::C_adm_fast(x), adm(x), tolerance = tol_n(1024))
})

# ============================================================
# WU-ADM1: adm_core kernel regression guards (AVX2 + 4-wide unroll)
# These tests verify the 4-wide dual-accumulator unroll and (1.0/n)
# multiplication produce the same result as the reference for all
# remainder patterns and edge cases. They must PASS both before AND
# after the kernel change — they are regression guards, not new-feature tests.
# ============================================================

# n=4..20: covers all 4-wide unroll remainder patterns (n mod 4 = 0,1,2,3)
# and the critical scalar-tail path at each boundary.
test_that("adm opt: WU-ADM1 regression guard — C_adm_fast matches adm() for n=4..20", {
  for (n in 4:20) {
    set.seed(100 + n); x <- rnorm(n)
    expect_equal(robscale:::C_adm_fast(x), adm(x), tolerance = tol_n(n),
                 label = paste0("C_adm_fast vs adm() at n=", n))
  }
})

# Tied values: many identical elements → abs deviations are all 0 or identical.
# Tests that accumulator chains don't diverge when subtraction yields 0.
test_that("adm opt: WU-ADM1 regression guard — C_adm_fast handles all-equal input", {
  x <- rep(3.7, 20)
  expect_equal(robscale:::C_adm_fast(x), adm(x), tolerance = tol_n(20))
})

test_that("adm opt: WU-ADM1 regression guard — C_adm_fast handles mostly-tied input", {
  x <- c(rep(5.0, 17), 4.0, 6.0, 7.0)
  expect_equal(robscale:::C_adm_fast(x), adm(x), tolerance = tol_n(20))
})

# Large n: exercises full 4-wide vectorised inner loop many times.
test_that("adm opt: WU-ADM1 regression guard — C_adm_fast matches adm() at n=4096", {
  set.seed(4096); x <- rnorm(4096)
  expect_equal(robscale:::C_adm_fast(x), adm(x), tolerance = tol_n(4096))
})

# Extreme values: ensures no overflow or cancellation in accumulators.
test_that("adm opt: WU-ADM1 regression guard — C_adm_fast handles large-magnitude input", {
  set.seed(7); x <- rnorm(16) * 1e14
  expect_equal(robscale:::C_adm_fast(x), adm(x), tolerance = tol_n(16))
})

# n not divisible by 4: scalar tail path exercised for each remainder.
test_that("adm opt: WU-ADM1 regression guard — C_adm_fast correct at n=5 (rem 1)", {
  set.seed(5); x <- rnorm(5)
  expect_equal(robscale:::C_adm_fast(x), adm(x), tolerance = tol_n(5))
})

test_that("adm opt: WU-ADM1 regression guard — C_adm_fast correct at n=6 (rem 2)", {
  set.seed(6); x <- rnorm(6)
  expect_equal(robscale:::C_adm_fast(x), adm(x), tolerance = tol_n(6))
})

test_that("adm opt: WU-ADM1 regression guard — C_adm_fast correct at n=7 (rem 3)", {
  set.seed(7); x <- rnorm(7)
  expect_equal(robscale:::C_adm_fast(x), adm(x), tolerance = tol_n(7))
})
