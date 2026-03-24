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
