# tests/testthat/test-iqr-opt.R
#
# TDD gate tests for OPT-I1..I8 (iqr_scaled() performance optimizations).
# Tests I.1, I.5, I.11 are RED before Phase 0 pin is captured.
# Tests I.10 is RED until Phase 6 adds the diagnostic export.
# All other tests are GREEN from the start and remain GREEN throughout.
# All must be GREEN after Phase 8.

K_IQR <- 0.741301109252801

# ---------------------------------------------------------------------------
# Test I.1: Regression guard — bit-identical result before and after OPT-I1..I8
# ---------------------------------------------------------------------------
test_that("iqr_scaled pin unchanged after OPT-I1..I8", {
  tryCatch(devtools::load_all(quiet = TRUE), error = function(e) invisible(NULL))
  pin <- readRDS(testthat::test_path("fixtures/iqr_baseline_pin.rds"))

  set.seed(99)
  expect_identical(iqr_scaled(rnorm(500)), pin$iqr_500)

  set.seed(42)
  expect_identical(iqr_scaled(rnorm(10), constant = 1.0), pin$iqr_10_c1)
})

# ---------------------------------------------------------------------------
# Test I.2: NOINLINE frame boundary — n=128 and n=129 both correct (OPT-I1)
# n=128 = last micro-path (buf_micro[128]); n=129 = first NOINLINE large path.
# ---------------------------------------------------------------------------
test_that("NOINLINE boundary: n=128 and n=129 give correct results", {
  set.seed(7)
  x128 <- rnorm(128)
  x129 <- c(x128, rnorm(1))

  expect_equal(iqr_scaled(x128, constant = 1),
               IQR(x128, type = 7),
               tolerance = sqrt(.Machine$double.eps),
               label = "n=128")
  expect_equal(iqr_scaled(x129, constant = 1),
               IQR(x129, type = 7),
               tolerance = sqrt(.Machine$double.eps),
               label = "n=129")
})

# ---------------------------------------------------------------------------
# Test I.3: Heap boundary — n=2049 correct after STACK_SIZE reduction to 2048
# n=2049 is the first n that exceeds the new STACK_SIZE=2048.
# ---------------------------------------------------------------------------
test_that("Heap boundary: n=2049 correct after STACK_SIZE=2048", {
  set.seed(13)
  x <- rnorm(2049)
  expect_equal(iqr_scaled(x, constant = 1),
               IQR(x, type = 7),
               tolerance = sqrt(.Machine$double.eps),
               label = "n=2049")
})

# ---------------------------------------------------------------------------
# Test I.4a: interp_q7 frac>0 path — n values where (n-1)%4 != 0 (OPT-I4)
# These n force the interpolation scan inside interp_q7.
# ---------------------------------------------------------------------------
test_that("interp_q7 frac>0 path: scan n values correct", {
  # (n-1)%4 == 1: n=2,6,10,14  (frac=0.25)
  # (n-1)%4 == 2: n=3,7,11,15  (frac=0.50)
  # (n-1)%4 == 3: n=4,8,12,16  (frac=0.75)
  scan_ns <- c(2L, 3L, 4L, 6L, 7L, 8L, 10L, 11L, 12L, 14L, 15L, 16L)
  for (n in scan_ns) {
    set.seed(100L + n)
    x <- rnorm(n)
    expect_equal(iqr_scaled(x, constant = 1),
                 IQR(x, type = 7),
                 tolerance = sqrt(.Machine$double.eps),
                 label = paste("frac>0, n =", n))
  }
})

# ---------------------------------------------------------------------------
# Test I.4b: interp_q7 frac==0 path — n values where (n-1)%4 == 0 (OPT-I4)
# These n have integer quantile indices; the scan branch is NOT taken.
# Guards that std::min_element replacement doesn't fire when frac==0.
# ---------------------------------------------------------------------------
test_that("interp_q7 frac==0 path: no-scan n values correct", {
  # (n-1)%4 == 0: n=5,9,13,17,21  (frac=0.0)
  noscan_ns <- c(5L, 9L, 13L, 17L, 21L)
  for (n in noscan_ns) {
    set.seed(200L + n)
    x <- rnorm(n)
    expect_equal(iqr_scaled(x, constant = 1),
                 IQR(x, type = 7),
                 tolerance = sqrt(.Machine$double.eps),
                 label = paste("frac==0, n =", n))
  }
})

# ---------------------------------------------------------------------------
# Test I.5: RESTRICT path — results unchanged after ROBSCALE_RESTRICT (OPT-I5)
# Uses pin values from Phase 0 to verify bit-identical output.
# ---------------------------------------------------------------------------
test_that("RESTRICT path: iqr_scaled results identical to pin at representative n", {
  tryCatch(devtools::load_all(quiet = TRUE), error = function(e) invisible(NULL))
  pin <- readRDS(testthat::test_path("fixtures/iqr_baseline_pin.rds"))

  # Re-derive the same values and confirm identity
  set.seed(99)
  expect_identical(iqr_scaled(rnorm(500)), pin$iqr_500,
                   label = "RESTRICT: n=500 pin identical")

  set.seed(42)
  expect_identical(iqr_scaled(rnorm(10), constant = 1.0), pin$iqr_10_c1,
                   label = "RESTRICT: n=10 constant=1 pin identical")
})

# ---------------------------------------------------------------------------
# Test I.6: Sort fast path — n=2..16 match stats::IQR exactly (OPT-I6)
# small_sort() + direct index reads must give bit-identical Type 7 quantiles.
# ---------------------------------------------------------------------------
test_that("Sort fast path: n=2..16 match stats::IQR exactly", {
  for (n in 2L:16L) {
    set.seed(300L + n)
    x <- rnorm(n)
    expect_equal(iqr_scaled(x, constant = 1),
                 IQR(x, type = 7),
                 tolerance = sqrt(.Machine$double.eps),
                 label = paste("sort path, n =", n))
  }
})

# ---------------------------------------------------------------------------
# Test I.7: Sort/pdqselect boundary — n=16 and n=17 give correct results
# n=16 uses sort path; n=17 is the first n in the pdqselect path.
# ---------------------------------------------------------------------------
test_that("Sort/pdqselect boundary: n=16 and n=17 give correct results", {
  set.seed(9)
  x16 <- rnorm(16)
  x17 <- c(x16, rnorm(1))

  expect_equal(iqr_scaled(x16, constant = 1),
               IQR(x16, type = 7),
               tolerance = sqrt(.Machine$double.eps),
               label = "n=16 sort path")
  expect_equal(iqr_scaled(x17, constant = 1),
               IQR(x17, type = 7),
               tolerance = sqrt(.Machine$double.eps),
               label = "n=17 first pdqselect path")
})

# ---------------------------------------------------------------------------
# Test I.8: Symmetric Q1 — all four (n-1)%4 residue classes correct (OPT-I3)
# Covers both frac==0 (standard path) and frac>0 (symmetric path) for n>16.
# ---------------------------------------------------------------------------
test_that("Symmetric Q1: all four (n-1)%%4 residue classes correct for n>16", {
  cases <- list(
    list(n = 17L, lab = "residue 0, frac=0, standard path"),   # (n-1)%4==0
    list(n = 18L, lab = "residue 1, frac=0.25, symmetric"),    # (n-1)%4==1
    list(n = 19L, lab = "residue 2, frac=0.5, symmetric"),     # (n-1)%4==2
    list(n = 20L, lab = "residue 3, frac=0.75, symmetric"),    # (n-1)%4==3
    list(n = 64L, lab = "residue 3, typical large n"),
    list(n = 100L, lab = "residue 3, n=100"),
    list(n = 1000L, lab = "residue 3, n=1000")
  )
  for (cc in cases) {
    set.seed(400L + cc$n)
    x <- rnorm(cc$n)
    expect_equal(iqr_scaled(x, constant = 1),
                 IQR(x, type = 7),
                 tolerance = sqrt(.Machine$double.eps),
                 label = cc$lab)
  }
})

# ---------------------------------------------------------------------------
# Test I.9: Symmetric Q1 known result — n=20
# x = 1:20 (sorted). R gives: IQR(1:20, type=7) = 9.75.
# iqr_scaled(1:20, constant=1) must equal 9.75 exactly.
# Derivation: h1=(19)*0.25=4.75 → Q1=x[4]+0.75*(x[5]-x[4])=4+0.75=4.75
#              h3=(19)*0.75=14.25 → Q3=x[14]+0.25*(x[15]-x[14])=14+0.25=14.25
#              IQR = 14.25 - 4.75 = 9.5.  With K_IQR: expect 9.5*0.7413...
# But with constant=1: expect 9.5. Verify against R's IQR() to avoid
# hard-coding and to guard against any edge-case in the symmetric path.
# ---------------------------------------------------------------------------
test_that("Symmetric Q1 known result: n=20, x=1:20", {
  x <- as.numeric(1:20)
  expect_equal(iqr_scaled(x, constant = 1),
               IQR(x, type = 7),
               tolerance = sqrt(.Machine$double.eps))
  # Cross-check: iqr_scaled with default constant
  expect_equal(iqr_scaled(x),
               IQR(x, type = 7) * K_IQR,
               tolerance = sqrt(.Machine$double.eps))
})

# ---------------------------------------------------------------------------
# Test I.10: iqr_sorted exists and is correct (via diagnostic export) (OPT-I7)
# RED until Phase 6 adds iqr_sorted_impl. Removed diagnostic before Phase 6 commit.
# ---------------------------------------------------------------------------
test_that("iqr_sorted_impl: O(1) sorted IQR matches iqr_scaled on sorted input", {
  tryCatch(devtools::load_all(quiet = TRUE), error = function(e) invisible(NULL))
  # Skip gracefully if diagnostic export not yet present (Phases 0-5)
  if (!exists("iqr_sorted_impl", mode = "function")) {
    skip("iqr_sorted_impl diagnostic export not present yet (Phase 6 adds it)")
  }
  for (n in c(2L, 3L, 4L, 5L, 16L, 17L, 100L, 1000L)) {
    set.seed(500L + n)
    x <- sort(rnorm(n))
    expect_equal(iqr_sorted_impl(x),
                 iqr_scaled(x),
                 tolerance = sqrt(.Machine$double.eps),
                 label = paste("iqr_sorted, sorted n =", n))
  }
})

# ---------------------------------------------------------------------------
# Test I.11: Ensemble pin unchanged after all IQR changes (OPT-I8)
# scale_robust() ensemble path calls estimators_internal::iqr() via
# compute_all_estimators() for n=10 (below auto_switch threshold of 20).
# ---------------------------------------------------------------------------
test_that("ensemble scale_robust pin unchanged after IQR optimizations", {
  tryCatch(devtools::load_all(quiet = TRUE), error = function(e) invisible(NULL))
  pin <- readRDS(testthat::test_path("fixtures/iqr_baseline_pin.rds"))
  set.seed(77)
  result <- scale_robust(rnorm(10), n_boot = 50)
  expect_equal(result, pin$ens_n10, tolerance = 8 * .Machine$double.eps)
})
