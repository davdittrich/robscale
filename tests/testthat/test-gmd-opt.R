library(testthat)

# ---------------------------------------------------------------------------
# TDD gate for GMD performance optimizations (OPT-G1/G3/G4/G5/G6)
# Phase 0 baseline pins stored in tests/testthat/fixtures/
# ---------------------------------------------------------------------------

test_that("1.1 — gmd() correctness: n = 2..2049 matches Phase 0 pin", {
  skip_if_not(
    file.exists(testthat::test_path("fixtures/gmd_baseline_pin.rds")),
    "Baseline pin missing — run Phase 0"
  )
  tol <- sqrt(.Machine$double.eps)
  pin <- readRDS(testthat::test_path("fixtures/gmd_baseline_pin.rds"))
  set.seed(2025)
  for (i in seq_along(pin$ns)) {
    n   <- pin$ns[i]
    got <- gmd(rnorm(n))
    expect_equal(got, pin$vals[i], tolerance = tol,
      label = sprintf("gmd(rnorm(%d))", n))
  }
})

test_that("1.2 — gmd() equals R-reference sort (sort-algorithm invariance)", {
  # Mathematical invariant: gmd is invariant to the relative ordering of equal
  # elements. For equal value v at positions i and j:
  #   w_i * v + w_j * v = (w_i + w_j) * v — unchanged by swapping i and j.
  # Therefore switching std::sort -> optimized_sort / boost::float_sort cannot
  # change gmd output for the same input multiset.
  tol <- sqrt(.Machine$double.eps)

  gmd_ref <- function(x, k = 0.886226925452758) {
    xs <- sort(x)
    n  <- length(xs)
    w  <- 2.0 * seq_len(n) - n - 1.0
    k * 2.0 * sum(w * xs) / (n * (n - 1.0))
  }

  # No ties
  set.seed(77)
  for (n in c(17L, 33L, 65L, 101L, 513L)) {
    x <- rnorm(n)
    expect_equal(gmd(x), gmd_ref(x), tolerance = tol,
      label = sprintf("no-ties n=%d", n))
  }

  # With ties (bootstrap-style resampling)
  set.seed(88)
  for (n in c(20L, 50L, 100L)) {
    x <- sample(rnorm(10L), n, replace = TRUE)
    expect_equal(gmd(x), gmd_ref(x), tolerance = tol,
      label = sprintf("with-ties n=%d", n))
  }
})

test_that("1.3 — ensemble value unchanged after sort upgrade (G4)", {
  skip_if_not(
    file.exists(testthat::test_path("fixtures/ens_baseline_pin.rds")),
    "Baseline pin missing — run Phase 0"
  )
  # XorShift32 seeded per replicate r: each writes to its own boot_results slot;
  # aggregation is serial post-loop. Output is bit-deterministic regardless of
  # TBB scheduling. expect_identical (strict equality) is correct.
  set.seed(42)
  v_after <- robscale:::cpp_scale_ensemble(rnorm(200L), 500L)
  v_pin   <- readRDS(testthat::test_path("fixtures/ens_baseline_pin.rds"))
  expect_identical(v_after, v_pin,
    label = "ensemble value must be bit-identical after sort algorithm change")
})

test_that("1.4 — G5 scale precompute: |before - after| <= 4 ULP", {
  skip_if_not(
    exists("gmd_scalar_kernel_impl",
           envir = asNamespace("robscale"), mode = "function"),
    "Diagnostic export gmd_scalar_kernel_impl not compiled (added in Phase 4)"
  )
  # constant * 2.0 * sum / (n*(n-1)) vs (constant * 2.0 / (n*(n-1))) * sum
  # differ by at most 1 ULP (one FP rounding step). 4x headroom.
  tol <- 4 * .Machine$double.eps
  set.seed(55)
  for (n in c(2L, 3L, 5L, 10L, 17L, 33L)) {
    x <- rnorm(n)
    v_current <- gmd(x)
    v_precomp <- robscale:::gmd_scalar_kernel_impl(x)
    expect_equal(v_current, v_precomp, tolerance = tol,
      label = sprintf("scale precompute ULP n=%d", n))
  }
})

test_that("1.5 — gmd() frame-split boundary: n=127..2049 correct", {
  tol <- sqrt(.Machine$double.eps)
  gmd_ref <- function(x, k = 0.886226925452758) {
    xs <- sort(x)
    n  <- length(xs)
    k * 2.0 * sum((2.0 * seq_len(n) - n - 1.0) * xs) / (n * (n - 1.0))
  }
  # n < 2 edge case
  expect_equal(gmd(1.0), 0.0)
  # Boundaries: last micro (128), first large-n (129), stack/heap boundary (2048/2049)
  for (n in c(127L, 128L, 129L, 2048L, 2049L)) {
    set.seed(n)
    x <- rnorm(n)
    expect_equal(gmd(x), gmd_ref(x), tolerance = tol,
      label = sprintf("frame-boundary n=%d", n))
  }
})
