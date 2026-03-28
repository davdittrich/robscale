# WU-2: rob_scale_core Implosion Guard (B2 — REVISED)
#
# Finding B2 claimed that 0 < s_init <= IMPLOSION_BOUND with fallback=0
# is a gap. Verification showed this is NOT a bug:
#
# 1. NR handles any positive s_init correctly (scale is relative)
# 2. ADM does NOT downweight outliers — redirecting to ADM when
#    MAD is small due to outliers (not near-constant data) gives WRONG results
# 3. Example: c(1e-4, 0, 0, 4) has s_init ≈ 7.41e-5 < IMPLOSION_BOUND.
#    NR correctly gives ~1e-4 (bulk scale). ADM gives ~1.25 (includes outlier).
#
# The existing guard is correct:
#   - fallback=1 and s_init <= implbound → NA (user's choice)
#   - s_init == 0.0 exactly → ADM (data is truly constant)
#   - Otherwise → NR (correct M-scale)
#
# These tests document and lock in the correct existing behavior.

test_that("NR runs correctly for small positive s_init (not redirected to ADM)", {
  # Data with outlier: MAD is small because bulk is near 0, not because data is constant.
  # NR should run and give the correct M-scale reflecting the bulk.
  x <- c(1e-4, 0, 0, 4)
  res <- robScale(x)
  expect_equal(res, 0.000101530115510382, tolerance = 1e-7,
               label = "NR runs for small s_init with outlier")
  # ADM would give ~1.25 here — completely wrong (doesn't downweight outlier)
  expect_true(res < 0.01,
              label = "M-scale reflects bulk, not outlier")
})

test_that("exact-zero MAD (truly constant data) returns ADM with fallback=adm", {
  x <- c(5, 5, 5, 5, 5)
  expect_equal(robScale(x), 0, label = "all-identical returns 0 via ADM")
})

test_that("exact-zero MAD returns NA with fallback=na", {
  x <- c(5, 5, 5, 5, 5)
  expect_true(is.na(robScale(x, fallback = "na")))
})

test_that("fallback=na returns NA when s_init in implosion range", {
  # This data has s_init ≈ 7.41e-5 < IMPLOSION_BOUND
  x <- c(1e-4, 0, 0, 4)
  expect_true(is.na(robScale(x, fallback = "na")))
})

test_that("near-constant data with positive MAD: NR converges", {
  # 5 values with tiny spread: MAD ≈ 1.48e-6, NR should converge
  x <- c(1, 1 + 1e-6, 1 + 2e-6, 1 + 3e-6, 1 + 4e-6)
  res <- robScale(x)
  expect_true(is.finite(res), label = "NR converges for near-constant data")
  expect_true(res > 0 && res < 1e-4,
              label = "scale is small but positive")
})

test_that("normal data above IMPLOSION_BOUND uses NR (not ADM)", {
  set.seed(42)
  x <- rnorm(20)
  res <- robScale(x)
  adm_val <- adm(x)
  expect_false(isTRUE(all.equal(res, adm_val, tolerance = 1e-4)),
               label = "normal data uses NR, not ADM fallback")
})
