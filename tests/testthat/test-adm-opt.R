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

# ============================================================
# WU-ADM2: caller-side guards (n=1, buffer boundaries, adaptive median)
# RED tests: must PASS after the caller changes are applied.
# ============================================================

test_that("adm opt: WU-ADM2 — C_adm_fast(numeric(1)) returns 0.0 (n=1 guard)", {
  expect_equal(robscale:::C_adm_fast(1.0), 0.0)
})

test_that("adm opt: WU-ADM2 — C_adm_fast correct at n=128 (last micro-buffer slot)", {
  set.seed(128); x <- rnorm(128)
  expect_equal(robscale:::C_adm_fast(x), adm(x), tolerance = tol_n(128))
})

test_that("adm opt: WU-ADM2 — C_adm_fast correct at n=129 (first large-n path)", {
  set.seed(129); x <- rnorm(129)
  expect_equal(robscale:::C_adm_fast(x), adm(x), tolerance = tol_n(129))
})

test_that("adm opt: WU-ADM2 — C_adm_fast correct at n=2048 (last stack-buf slot)", {
  set.seed(2048); x <- rnorm(2048)
  expect_equal(robscale:::C_adm_fast(x), adm(x), tolerance = tol_n(2048))
})

test_that("adm opt: WU-ADM2 — C_adm_fast correct at n=2049 (first heap allocation)", {
  set.seed(2049); x <- rnorm(2049)
  expect_equal(robscale:::C_adm_fast(x), adm(x), tolerance = tol_n(2049))
})

# ============================================================
# WU-ADM3: adm_core_sorted correctness (RED tests — will error before
# C_adm_core_sorted is exported; pass after implementation).
# ============================================================

# Helper: compute ADM reference from unsorted x (uses adm_core path)
adm_ref_sorted <- function(x) {
  xs <- sort(x)
  robscale:::C_adm_core_sorted(xs)
}

test_that("adm opt: WU-ADM3 — C_adm_core_sorted matches adm() at n=1..5,8,9,16,17,100", {
  for (n in c(1L, 2L, 3L, 4L, 5L, 8L, 9L, 16L, 17L, 100L)) {
    set.seed(200L + n); x <- rnorm(n)
    expect_equal(adm_ref_sorted(x), adm(x), tolerance = tol_n(n),
                 label = paste0("adm_core_sorted vs adm() at n=", n))
  }
})

test_that("adm opt: WU-ADM3 — adm_core_sorted n=1 returns 0.0", {
  expect_equal(robscale:::C_adm_core_sorted(1.0), 0.0)
})

test_that("adm opt: WU-ADM3 — adm_core_sorted correct for odd n (n=5)", {
  x <- sort(c(1.0, 2.0, 3.0, 4.0, 5.0))
  # med=3, lower=[1,2]->3, upper=[4,5]->9, (9-3)/5=1.2; consistency*1.2
  expect_equal(robscale:::C_adm_core_sorted(x),
               adm(x), tolerance = tol_n(5))
})

test_that("adm opt: WU-ADM3 — adm_core_sorted correct for even n (n=4)", {
  x <- sort(c(1.0, 2.0, 4.0, 5.0))
  # med=3 (avg 2,4), lower=[1,2]->3, upper=[4,5]->9, (9-3)/4=1.5; consistency*1.5
  expect_equal(robscale:::C_adm_core_sorted(x),
               adm(x), tolerance = tol_n(4))
})

test_that("adm opt: WU-ADM3 — adm_core_sorted correct for tied values", {
  x <- sort(c(1.0, 3.0, 3.0, 5.0))
  expect_equal(robscale:::C_adm_core_sorted(x),
               adm(x), tolerance = tol_n(4))
})

test_that("adm opt: WU-ADM3 — adm_core_sorted correct for all-equal input", {
  x <- rep(3.7, 10)
  expect_equal(robscale:::C_adm_core_sorted(x), 0.0, tolerance = tol_n(10))
})
