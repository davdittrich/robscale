## test-rob-scale-opt.R
## Correctness guards for robScale() performance optimization session.
## These tests must pass before AND after every WU; they detect regressions
## introduced by structural or algorithmic changes.
##
## Coverage:
##   - WU-RS0: diagnostic export existence + equivalence guards + edge cases
##   - WU-RS1 (pre-written): remainder-path and ADM-fallback guards

library(robscale)
tol_n <- function(n) n * .Machine$double.eps

# ============================================================
# 1. Diagnostic export presence
# ============================================================

test_that("rob_scale opt: C_rob_scale_orig exists as exported function", {
  expect_true(exists("C_rob_scale_orig",
                     envir = asNamespace("robscale"), inherits = FALSE))
})

test_that("rob_scale opt: C_rob_scale_fast exists as exported function", {
  expect_true(exists("C_rob_scale_fast",
                     envir = asNamespace("robscale"), inherits = FALSE))
})

# ============================================================
# 2. Equivalence: C_rob_scale_fast == C_rob_scale_orig
#    (Pass before AND after WU-RS1; fail only on numerical error > tol)
# ============================================================

test_that("rob_scale opt: fast==orig for n=4", {
  set.seed(401); x <- rnorm(4)
  expect_equal(robscale:::C_rob_scale_fast(x), robscale:::C_rob_scale_orig(x),
               tolerance = tol_n(4))
})

test_that("rob_scale opt: fast==orig for n=5", {
  set.seed(501); x <- rnorm(5)
  expect_equal(robscale:::C_rob_scale_fast(x), robscale:::C_rob_scale_orig(x),
               tolerance = tol_n(5))
})

test_that("rob_scale opt: fast==orig for n=6", {
  set.seed(601); x <- rnorm(6)
  expect_equal(robscale:::C_rob_scale_fast(x), robscale:::C_rob_scale_orig(x),
               tolerance = tol_n(6))
})

test_that("rob_scale opt: fast==orig for n=7", {
  set.seed(701); x <- rnorm(7)
  expect_equal(robscale:::C_rob_scale_fast(x), robscale:::C_rob_scale_orig(x),
               tolerance = tol_n(7))
})

test_that("rob_scale opt: fast==orig for n=8", {
  set.seed(801); x <- rnorm(8)
  expect_equal(robscale:::C_rob_scale_fast(x), robscale:::C_rob_scale_orig(x),
               tolerance = tol_n(8))
})

test_that("rob_scale opt: fast==orig for n=9", {
  set.seed(901); x <- rnorm(9)
  expect_equal(robscale:::C_rob_scale_fast(x), robscale:::C_rob_scale_orig(x),
               tolerance = tol_n(9))
})

test_that("rob_scale opt: fast==orig for n=10", {
  set.seed(1001); x <- rnorm(10)
  expect_equal(robscale:::C_rob_scale_fast(x), robscale:::C_rob_scale_orig(x),
               tolerance = tol_n(10))
})

test_that("rob_scale opt: fast==orig for n=12", {
  set.seed(1201); x <- rnorm(12)
  expect_equal(robscale:::C_rob_scale_fast(x), robscale:::C_rob_scale_orig(x),
               tolerance = tol_n(12))
})

test_that("rob_scale opt: fast==orig for n=15", {
  set.seed(1501); x <- rnorm(15)
  expect_equal(robscale:::C_rob_scale_fast(x), robscale:::C_rob_scale_orig(x),
               tolerance = tol_n(15))
})

test_that("rob_scale opt: fast==orig for n=16 (is_small boundary)", {
  set.seed(1601); x <- rnorm(16)
  expect_equal(robscale:::C_rob_scale_fast(x), robscale:::C_rob_scale_orig(x),
               tolerance = tol_n(16))
})

test_that("rob_scale opt: fast==orig for n=17 (above is_small boundary)", {
  set.seed(1701); x <- rnorm(17)
  expect_equal(robscale:::C_rob_scale_fast(x), robscale:::C_rob_scale_orig(x),
               tolerance = tol_n(17))
})

test_that("rob_scale opt: fast==orig for n=30", {
  set.seed(3001); x <- rnorm(30)
  expect_equal(robscale:::C_rob_scale_fast(x), robscale:::C_rob_scale_orig(x),
               tolerance = tol_n(30))
})

test_that("rob_scale opt: fast==orig for n=50", {
  set.seed(5001); x <- rnorm(50)
  expect_equal(robscale:::C_rob_scale_fast(x), robscale:::C_rob_scale_orig(x),
               tolerance = tol_n(50))
})

test_that("rob_scale opt: fast==orig for n=64 (small-n arena boundary)", {
  set.seed(6401); x <- rnorm(64)
  expect_equal(robscale:::C_rob_scale_fast(x), robscale:::C_rob_scale_orig(x),
               tolerance = tol_n(64))
})

test_that("rob_scale opt: fast==orig for n=65 (above small-n arena boundary)", {
  set.seed(6501); x <- rnorm(65)
  expect_equal(robscale:::C_rob_scale_fast(x), robscale:::C_rob_scale_orig(x),
               tolerance = tol_n(65))
})

test_that("rob_scale opt: fast==orig for n=100", {
  set.seed(10001); x <- rnorm(100)
  expect_equal(robscale:::C_rob_scale_fast(x), robscale:::C_rob_scale_orig(x),
               tolerance = tol_n(100))
})

test_that("rob_scale opt: fast==orig for n=200", {
  set.seed(20001); x <- rnorm(200)
  expect_equal(robscale:::C_rob_scale_fast(x), robscale:::C_rob_scale_orig(x),
               tolerance = tol_n(200))
})

test_that("rob_scale opt: fast==orig for n=500", {
  set.seed(50001); x <- rnorm(500)
  expect_equal(robscale:::C_rob_scale_fast(x), robscale:::C_rob_scale_orig(x),
               tolerance = tol_n(500))
})

test_that("rob_scale opt: fast==orig for n=1000", {
  set.seed(100001); x <- rnorm(1000)
  expect_equal(robscale:::C_rob_scale_fast(x), robscale:::C_rob_scale_orig(x),
               tolerance = tol_n(1000))
})

# ============================================================
# 3. Stack/heap boundary straddle at SCALE_STACK_SIZE = 2048
# ============================================================

test_that("rob_scale opt: fast==orig for n=2047 (stack path)", {
  set.seed(204701); x <- rnorm(2047)
  expect_equal(robscale:::C_rob_scale_fast(x), robscale:::C_rob_scale_orig(x),
               tolerance = tol_n(2047))
})

test_that("rob_scale opt: fast==orig for n=2048 (stack/heap boundary)", {
  set.seed(204801); x <- rnorm(2048)
  expect_equal(robscale:::C_rob_scale_fast(x), robscale:::C_rob_scale_orig(x),
               tolerance = tol_n(2048))
})

test_that("rob_scale opt: fast==orig for n=2049 (heap path)", {
  set.seed(204901); x <- rnorm(2049)
  expect_equal(robscale:::C_rob_scale_fast(x), robscale:::C_rob_scale_orig(x),
               tolerance = tol_n(2049))
})

# ============================================================
# 4. C_rob_scale_fast matches production robScale()
# ============================================================

test_that("rob_scale opt: C_rob_scale_fast matches robScale for n=4", {
  set.seed(42); x <- rnorm(4)
  expect_equal(robscale:::C_rob_scale_fast(x), robScale(x), tolerance = tol_n(4))
})

test_that("rob_scale opt: C_rob_scale_fast matches robScale for n=20", {
  set.seed(42); x <- rnorm(20)
  expect_equal(robscale:::C_rob_scale_fast(x), robScale(x), tolerance = tol_n(20))
})

test_that("rob_scale opt: C_rob_scale_fast matches robScale for n=100", {
  set.seed(42); x <- rnorm(100)
  expect_equal(robscale:::C_rob_scale_fast(x), robScale(x), tolerance = tol_n(100))
})

# ============================================================
# 5. Edge cases
# ============================================================

test_that("rob_scale opt: n=0 returns 0 via C_rob_scale_fast", {
  expect_equal(robscale:::C_rob_scale_fast(numeric(0)), 0.0)
})

test_that("rob_scale opt: n=0 returns 0 via C_rob_scale_orig", {
  expect_equal(robscale:::C_rob_scale_orig(numeric(0)), 0.0)
})

test_that("rob_scale opt: n=1 constant data returns 0 via C_rob_scale_fast", {
  # MAD=0 → ADM fallback → 0 for single element
  expect_equal(robscale:::C_rob_scale_fast(c(5.0)), 0.0)
})

test_that("rob_scale opt: n=1 constant data returns 0 via C_rob_scale_orig", {
  expect_equal(robscale:::C_rob_scale_orig(c(5.0)), 0.0)
})

test_that("rob_scale opt: n<4 fast==orig (early-exit branch)", {
  # Both must agree on the early-exit path value regardless of what it is
  set.seed(42)
  for (n in c(1L, 2L, 3L)) {
    x <- rnorm(n)
    expect_equal(robscale:::C_rob_scale_fast(x),
                 robscale:::C_rob_scale_orig(x),
                 tolerance = 8 * .Machine$double.eps,
                 label = paste0("n=", n))
  }
})

test_that("rob_scale opt: MAD implosion (all-tied data) triggers ADM fallback via C_rob_scale_fast", {
  # 8 out of 9 elements identical → MAD = 0 → ADM fallback
  x <- c(rep(5.0, 8), 6.0)
  fast_val <- robscale:::C_rob_scale_fast(x)
  orig_val <- robscale:::C_rob_scale_orig(x)
  expect_equal(fast_val, orig_val, tolerance = 8 * .Machine$double.eps)
  expect_gt(fast_val, 0.0)  # ADM returns a positive value
})

test_that("rob_scale opt: exactly constant data returns 0 for n=20", {
  x <- rep(3.14159, 20)
  expect_equal(robscale:::C_rob_scale_fast(x), 0.0)
  expect_equal(robscale:::C_rob_scale_orig(x), 0.0)
})

test_that("rob_scale opt: determinism — 10 identical calls return same value", {
  set.seed(99); x <- rnorm(100)
  vals <- replicate(10, robscale:::C_rob_scale_fast(x))
  expect_true(all(vals == vals[1]))
})

test_that("rob_scale opt: result is positive for typical data", {
  set.seed(42); x <- rnorm(50)
  expect_gt(robscale:::C_rob_scale_fast(x), 0.0)
  expect_gt(robscale:::C_rob_scale_orig(x), 0.0)
})

test_that("rob_scale opt: scaled data gives proportionally scaled result", {
  set.seed(42); x <- rnorm(50)
  fast1 <- robscale:::C_rob_scale_fast(x)
  fast2 <- robscale:::C_rob_scale_fast(x * 10)
  expect_equal(fast2 / fast1, 10.0, tolerance = 1e-6)
})

# ============================================================
# 6. WU-RS1 pre-written guards (written before implementing 8-wide AVX2)
#    These pass before AND after WU-RS1; fail only if the 8-wide
#    accumulator introduces a numerical error > tolerance.
# ============================================================

test_that("rob_scale 8-wide AVX2 pre-guard: n=5 (1 leftover in 4-wide cleanup)", {
  set.seed(5); x <- rnorm(5)
  expect_equal(robscale:::C_rob_scale_fast(x), robscale:::C_rob_scale_orig(x),
               tolerance = tol_n(5), label = "n=5")
})

test_that("rob_scale 8-wide AVX2 pre-guard: n=7 (3 leftover → scalar tail)", {
  set.seed(7); x <- rnorm(7)
  expect_equal(robscale:::C_rob_scale_fast(x), robscale:::C_rob_scale_orig(x),
               tolerance = tol_n(7), label = "n=7")
})

test_that("rob_scale 8-wide AVX2 pre-guard: n=9 (1 leftover after first 8-wide block)", {
  set.seed(9); x <- rnorm(9)
  expect_equal(robscale:::C_rob_scale_fast(x), robscale:::C_rob_scale_orig(x),
               tolerance = tol_n(9), label = "n=9")
})

test_that("rob_scale 8-wide AVX2 pre-guard: n=11 (3 leftover)", {
  set.seed(11); x <- rnorm(11)
  expect_equal(robscale:::C_rob_scale_fast(x), robscale:::C_rob_scale_orig(x),
               tolerance = tol_n(11), label = "n=11")
})

test_that("rob_scale 8-wide AVX2 pre-guard: n=13 (5 leftover → 4-wide + 1 scalar)", {
  set.seed(13); x <- rnorm(13)
  expect_equal(robscale:::C_rob_scale_fast(x), robscale:::C_rob_scale_orig(x),
               tolerance = tol_n(13), label = "n=13")
})

test_that("rob_scale 8-wide AVX2 pre-guard: n=15 (7 leftover → 4-wide + 3 scalar)", {
  set.seed(15); x <- rnorm(15)
  expect_equal(robscale:::C_rob_scale_fast(x), robscale:::C_rob_scale_orig(x),
               tolerance = tol_n(15), label = "n=15")
})

test_that("rob_scale 8-wide AVX2 pre-guard: ADM fallback path unchanged by kernel change", {
  # All-equal data: MAD=0 → ADM. The AVX2 kernel is not reached.
  expect_equal(robscale:::C_rob_scale_fast(rep(3.0, 20)),
               robscale:::C_rob_scale_orig(rep(3.0, 20)),
               tolerance = .Machine$double.eps)
})

# ============================================================
# 7. Removal guard for WU-RS2a (added here as pending test)
#    This test MUST be added AFTER WU-RS1 removes C_rob_scale_orig.
#    It is currently left commented out — uncomment when applying WU-RS2a.
# ============================================================
# test_that("rob_scale opt: C_rob_scale_orig removed after WU-RS1 (WU-RS2a)", {
#   expect_false(exists("C_rob_scale_orig",
#                       envir = asNamespace("robscale"), inherits = FALSE))
# })
