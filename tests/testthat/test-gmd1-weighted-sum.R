library(testthat)

# WU-GMD-1 RED: correctness tests for gmd_weighted_sum refactoring.
# These tests guard against regressions when the 3 inline weighted-sum
# loops are replaced by a shared AVX2 FMA kernel in robust_core.h.
# All tests must pass before AND after the kernel is introduced.

GMD_CONSISTENCY <- 0.886226925452758  # sqrt(pi)/2

scalar_gmd_ref <- function(x, constant = GMD_CONSISTENCY) {
  x <- sort(x)
  n <- length(x)
  if (n < 2) return(0.0)
  scale <- constant * 2.0 / (n * (n - 1.0))
  w <- 2.0 * seq_len(n) - n - 1.0  # 2*(i+1) - n - 1
  scale * sum(w * x)
}

# Tolerance for FMA accumulation rounding: n * eps
tol_n <- function(n) n * .Machine$double.eps

# --- Boundary: n < 2 ---

test_that("gmd1 — gmd_impl(numeric(0)) == 0", {
  expect_equal(robscale:::gmd_impl(numeric(0), GMD_CONSISTENCY), 0.0)
})

test_that("gmd1 — gmd_impl(scalar) == 0", {
  expect_equal(robscale:::gmd_impl(42.0, GMD_CONSISTENCY), 0.0)
})

test_that("gmd1 — gmd_impl(n=2) exact", {
  x <- c(1.0, 3.0)
  expect_equal(robscale:::gmd_impl(x, GMD_CONSISTENCY),
               scalar_gmd_ref(x), tolerance = 0)
})

# --- mod-4 boundary sizes (AVX2 tail handling) ---

for (sz in c(4L, 5L, 6L, 7L, 8L, 9L, 10L, 11L, 12L)) {
  local({
    n <- sz
    test_that(paste0("gmd1 — gmd_impl mod-4, n=", n), {
      set.seed(n)
      x <- rnorm(n)
      expect_equal(robscale:::gmd_impl(x, GMD_CONSISTENCY),
                   scalar_gmd_ref(x), tolerance = tol_n(n))
    })
  })
}

# --- sort-network threshold crossing (n=16 vs n=17) ---

test_that("gmd1 — gmd_impl n=16 (sort network)", {
  set.seed(16)
  x <- rnorm(16)
  expect_equal(robscale:::gmd_impl(x, GMD_CONSISTENCY),
               scalar_gmd_ref(x), tolerance = tol_n(16))
})

test_that("gmd1 — gmd_impl n=17 (first optimized_sort)", {
  set.seed(17)
  x <- rnorm(17)
  expect_equal(robscale:::gmd_impl(x, GMD_CONSISTENCY),
               scalar_gmd_ref(x), tolerance = tol_n(17))
})

# --- General sizes ---

for (sz in c(100L, 101L, 102L, 103L, 104L, 1000L, 10000L)) {
  local({
    n <- sz
    test_that(paste0("gmd1 — gmd_impl n=", n), {
      set.seed(n)
      x <- rnorm(n)
      expect_equal(robscale:::gmd_impl(x, GMD_CONSISTENCY),
                   scalar_gmd_ref(x), tolerance = tol_n(n))
    })
  })
}

# --- Equivariance ---

test_that("gmd1 — scale equivariant: gmd(a*x) == a * gmd(x)", {
  set.seed(777)
  x <- rnorm(200)
  a <- 3.14159
  expect_equal(robscale:::gmd_impl(a * x, GMD_CONSISTENCY),
               a * robscale:::gmd_impl(x, GMD_CONSISTENCY),
               tolerance = tol_n(200))
})

test_that("gmd1 — translation invariant: gmd(x+b) == gmd(x)", {
  # GMD is translation-invariant; use modest b to avoid FP cancellation.
  # Large b (e.g., 1e4) causes ~b*n*eps rounding loss in the weighted sum
  # because values of order 1e4 are summed with alternating signs.
  set.seed(888)
  x <- rnorm(200)
  b <- 1.0
  expect_equal(robscale:::gmd_impl(x + b, GMD_CONSISTENCY),
               robscale:::gmd_impl(x, GMD_CONSISTENCY),
               tolerance = tol_n(200))
})

# --- All-equal input (degenerate) ---

test_that("gmd1 — gmd_impl all-equal returns 0", {
  x <- rep(5.0, 50)
  expect_equal(robscale:::gmd_impl(x, GMD_CONSISTENCY), 0.0)
})

# --- Consistency with high-level gmd() R wrapper ---

test_that("gmd1 — gmd_impl consistent with gmd() wrapper, n=100", {
  set.seed(42)
  x <- rnorm(100)
  # gmd() uses GMD_CONSISTENCY internally
  expect_equal(robscale:::gmd_impl(x, GMD_CONSISTENCY),
               gmd(x), tolerance = tol_n(100))
})

test_that("gmd1 — gmd_impl consistent with gmd() wrapper, n=1000", {
  set.seed(43)
  x <- rnorm(1000)
  expect_equal(robscale:::gmd_impl(x, GMD_CONSISTENCY),
               gmd(x), tolerance = tol_n(1000))
})

# --- Zero-sum property: sorted data symmetric around 0 ---

test_that("gmd1 — gmd_impl symmetric data: value positive", {
  x <- c(-3, -1, 1, 3)
  val <- robscale:::gmd_impl(x, GMD_CONSISTENCY)
  ref <- scalar_gmd_ref(x)
  expect_equal(val, ref, tolerance = 0)
  expect_true(val > 0)
})

# --- Negative-valued input ---

test_that("gmd1 — gmd_impl all-negative matches scalar ref", {
  set.seed(555)
  x <- -abs(rnorm(100))
  expect_equal(robscale:::gmd_impl(x, GMD_CONSISTENCY),
               scalar_gmd_ref(x), tolerance = tol_n(100))
})
