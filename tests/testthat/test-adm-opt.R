## test-adm-opt.R
## Correctness and regression guards for adm() optimization session.
## WU-ADM1: adm_core kernel correctness (regression guards post-AVX2 + 4-wide unroll).
## WU-ADM2: caller-side guards (n=1, buffer boundaries, adaptive median).
## WU-ADM3: adm_core_sorted correctness and ensemble regression.
## WU-ADM4: large-n correctness, determinism, diagnostic export cleanup guards.

library(robscale)
tol_n <- function(n) n * .Machine$double.eps

## Pure R reference: ADM = constant * mean(|x - median(x)|)
ADM_CONSISTENCY <- 1.2533141373155001  # sqrt(pi/2)
adm_r_ref <- function(x) ADM_CONSISTENCY * mean(abs(x - median(x)))

## Pure R reference for adm_core_sorted (sorted input formula).
## For sorted x, center = exact median:
##   sum|x_i - center| = upper_sum - lower_sum
##   lower_sum = sum(x[0..k-1]),  upper_sum = sum(x[k+(n&1)..n-1]),  k = n/2
adm_ref_sorted <- function(x) {
  xs <- sort(x)
  n  <- length(xs)
  if (n <= 1L) return(0.0)
  k  <- n %/% 2L
  lower_sum <- if (k > 0L) sum(xs[seq_len(k)]) else 0.0
  upper_start <- k + (n %% 2L) + 1L   # 1-indexed; skip median element for odd n
  upper_sum <- if (upper_start <= n) sum(xs[seq(upper_start, n)]) else 0.0
  ADM_CONSISTENCY * (upper_sum - lower_sum) / n
}

# ============================================================
# WU-ADM1: adm_core kernel regression guards (AVX2 + 4-wide unroll)
# Verify the 4-wide dual-accumulator unroll and (1.0/n) multiplication
# produce the correct result for all remainder patterns and edge cases.
# ============================================================

# n=4..20: covers all 4-wide unroll remainder patterns (n mod 4 = 0,1,2,3)
# and the critical scalar-tail path at each boundary.
test_that("adm opt: WU-ADM1 regression guard — adm() matches R reference for n=4..20", {
  for (n in 4:20) {
    set.seed(100 + n); x <- rnorm(n)
    expect_equal(adm(x), adm_r_ref(x), tolerance = tol_n(n),
                 label = paste0("adm() vs R reference at n=", n))
  }
})

# Tied values: many identical elements → abs deviations are all 0 or identical.
test_that("adm opt: WU-ADM1 regression guard — adm() handles all-equal input", {
  x <- rep(3.7, 20)
  expect_equal(adm(x), adm_r_ref(x), tolerance = tol_n(20))
})

test_that("adm opt: WU-ADM1 regression guard — adm() handles mostly-tied input", {
  x <- c(rep(5.0, 17), 4.0, 6.0, 7.0)
  expect_equal(adm(x), adm_r_ref(x), tolerance = tol_n(20))
})

# Large n: exercises full 4-wide vectorised inner loop many times.
test_that("adm opt: WU-ADM1 regression guard — adm() matches R reference at n=4096", {
  set.seed(4096); x <- rnorm(4096)
  expect_equal(adm(x), adm_r_ref(x), tolerance = tol_n(4096))
})

# Extreme values: ensures no overflow or cancellation in accumulators.
test_that("adm opt: WU-ADM1 regression guard — adm() handles large-magnitude input", {
  set.seed(7); x <- rnorm(16) * 1e14
  expect_equal(adm(x), adm_r_ref(x), tolerance = tol_n(16))
})

# n not divisible by 4: scalar tail path exercised for each remainder.
test_that("adm opt: WU-ADM1 regression guard — adm() correct at n=5 (rem 1)", {
  set.seed(5); x <- rnorm(5)
  expect_equal(adm(x), adm_r_ref(x), tolerance = tol_n(5))
})

test_that("adm opt: WU-ADM1 regression guard — adm() correct at n=6 (rem 2)", {
  set.seed(6); x <- rnorm(6)
  expect_equal(adm(x), adm_r_ref(x), tolerance = tol_n(6))
})

test_that("adm opt: WU-ADM1 regression guard — adm() correct at n=7 (rem 3)", {
  set.seed(7); x <- rnorm(7)
  expect_equal(adm(x), adm_r_ref(x), tolerance = tol_n(7))
})

# ============================================================
# WU-ADM2: caller-side guards (n=1, buffer boundaries, adaptive median)
# ============================================================

test_that("adm opt: WU-ADM2 — adm(numeric(1)) returns 0.0 (n=1 guard)", {
  expect_equal(adm(1.0), 0.0)
})

test_that("adm opt: WU-ADM2 — adm() correct at n=128 (last micro-buffer slot)", {
  set.seed(128); x <- rnorm(128)
  expect_equal(adm(x), adm_r_ref(x), tolerance = tol_n(128))
})

test_that("adm opt: WU-ADM2 — adm() correct at n=129 (first large-n path)", {
  set.seed(129); x <- rnorm(129)
  expect_equal(adm(x), adm_r_ref(x), tolerance = tol_n(129))
})

test_that("adm opt: WU-ADM2 — adm() correct at n=2048 (last stack-buf slot)", {
  set.seed(2048); x <- rnorm(2048)
  expect_equal(adm(x), adm_r_ref(x), tolerance = tol_n(2048))
})

test_that("adm opt: WU-ADM2 — adm() correct at n=2049 (first heap allocation)", {
  set.seed(2049); x <- rnorm(2049)
  expect_equal(adm(x), adm_r_ref(x), tolerance = tol_n(2049))
})

# ============================================================
# WU-ADM3: adm_core_sorted correctness
# Pure R reference (adm_ref_sorted) verifies the sorted-path formula.
# ============================================================

test_that("adm opt: WU-ADM3 — adm_core_sorted matches adm() at n=1..5,8,9,16,17,100", {
  for (n in c(1L, 2L, 3L, 4L, 5L, 8L, 9L, 16L, 17L, 100L)) {
    set.seed(200L + n); x <- rnorm(n)
    expect_equal(adm_ref_sorted(x), adm(x), tolerance = tol_n(n),
                 label = paste0("adm_core_sorted ref vs adm() at n=", n))
  }
})

test_that("adm opt: WU-ADM3 — adm_core_sorted formula: n=1 returns 0.0", {
  expect_equal(adm_ref_sorted(1.0), 0.0)
})

test_that("adm opt: WU-ADM3 — adm_core_sorted formula correct for odd n (n=5)", {
  x <- sort(c(1.0, 2.0, 3.0, 4.0, 5.0))
  # med=3, lower=[1,2]->3, upper=[4,5]->9, (9-3)/5=1.2; consistency*1.2
  expect_equal(adm_ref_sorted(x), adm(x), tolerance = tol_n(5))
})

test_that("adm opt: WU-ADM3 — adm_core_sorted formula correct for even n (n=4)", {
  x <- sort(c(1.0, 2.0, 4.0, 5.0))
  # med=3 (avg 2,4), lower=[1,2]->3, upper=[4,5]->9, (9-3)/4=1.5; consistency*1.5
  expect_equal(adm_ref_sorted(x), adm(x), tolerance = tol_n(4))
})

test_that("adm opt: WU-ADM3 — adm_core_sorted formula correct for tied values", {
  x <- sort(c(1.0, 3.0, 3.0, 5.0))
  expect_equal(adm_ref_sorted(x), adm(x), tolerance = tol_n(4))
})

test_that("adm opt: WU-ADM3 — adm_core_sorted formula correct for all-equal input", {
  x <- rep(3.7, 10)
  expect_equal(adm_ref_sorted(x), 0.0, tolerance = tol_n(10))
})

# ============================================================
# WU-ADM4: large-n correctness, determinism, diagnostic export cleanup
# ============================================================

test_that("adm opt: WU-ADM4 — adm() large-n correctness (n=65536)", {
  set.seed(42L); x <- rnorm(65536L)
  ref <- adm_r_ref(x)
  expect_equal(adm(x), ref, tolerance = 65536L * .Machine$double.eps,
               label = "adm() n=65536 vs R reference")
})

test_that("adm opt: WU-ADM4 — adm() large-n correctness (n=131072)", {
  set.seed(42L); x <- rnorm(131072L)
  ref <- adm_r_ref(x)
  expect_equal(adm(x), ref, tolerance = 131072L * .Machine$double.eps,
               label = "adm() n=131072 vs R reference")
})

test_that("adm opt: WU-ADM4 — adm() large-n determinism (10 replicate calls)", {
  set.seed(42L); x <- rnorm(131072L)
  results <- vapply(seq_len(10L), function(i) adm(x), numeric(1L))
  expect_true(all(results == results[1L]),
              label = "adm() large-n returns identical value across 10 calls")
})

## Removal guards: diagnostic exports must be absent after WU-ADM4 cleanup
test_that("adm opt: WU-ADM4 — C_adm_orig removed from namespace", {
  expect_false(exists("C_adm_orig",
                      envir = asNamespace("robscale"), inherits = FALSE),
               label = "C_adm_orig must be absent after WU-ADM4 cleanup")
})

test_that("adm opt: WU-ADM4 — C_adm_fast removed from namespace", {
  expect_false(exists("C_adm_fast",
                      envir = asNamespace("robscale"), inherits = FALSE),
               label = "C_adm_fast must be absent after WU-ADM4 cleanup")
})

test_that("adm opt: WU-ADM4 — C_adm_core_sorted removed from namespace", {
  expect_false(exists("C_adm_core_sorted",
                      envir = asNamespace("robscale"), inherits = FALSE),
               label = "C_adm_core_sorted must be absent after WU-ADM4 cleanup")
})
