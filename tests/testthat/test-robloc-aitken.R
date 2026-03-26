## Tests 0.36–0.44 — robLoc Aitken Δ² acceleration (WU-RL-A1 correctness baseline)
##
## rob_loc_noaitken_impl: permanent correctness cross-check (like rob_loc_scalar_impl).
## Tests 0.37–0.44 remain GREEN both before and after WU-RL-A2 (same fixed point).
##
## Tests 0.45–0.47 (TDD gate + removal guards) are added at the start of WU-RL-A2.

# ── Test 0.36 — diagnostic callability ───────────────────────────────────────
test_that("0.36 rob_loc_noaitken_impl is callable", {
  skip_if_not(
    exists("rob_loc_noaitken_impl", envir = asNamespace("robscale"), mode = "function"),
    "rob_loc_noaitken_impl diagnostic not compiled"
  )
  x <- rnorm(20, mean = 3, sd = 2)
  result <- robscale:::rob_loc_noaitken_impl(x)
  expect_true(is.finite(result))
})

# ── Tests 0.37–0.44 — correctness contract ───────────────────────────────────

test_that("0.37 Aitken result matches plain NR within tolerance, n=50", {
  skip_if_not(
    exists("rob_loc_noaitken_impl", envir = asNamespace("robscale"), mode = "function"),
    "rob_loc_noaitken_impl diagnostic not compiled"
  )
  set.seed(101)
  x <- rnorm(50, mean = 5, sd = 3)
  ref  <- robscale:::rob_loc_noaitken_impl(x)
  aitk <- robLoc(x)
  expect_equal(aitk, ref, tolerance = 2 * sqrt(.Machine$double.eps))
})

test_that("0.38 Aitken result matches plain NR, n=200", {
  skip_if_not(
    exists("rob_loc_noaitken_impl", envir = asNamespace("robscale"), mode = "function"),
    "rob_loc_noaitken_impl diagnostic not compiled"
  )
  set.seed(202)
  x <- rnorm(200, mean = -2, sd = 5)
  expect_equal(robLoc(x), robscale:::rob_loc_noaitken_impl(x),
               tolerance = 2 * sqrt(.Machine$double.eps))
})

test_that("0.39 Aitken handles negative location, n=30", {
  skip_if_not(
    exists("rob_loc_noaitken_impl", envir = asNamespace("robscale"), mode = "function"),
    "rob_loc_noaitken_impl diagnostic not compiled"
  )
  set.seed(303)
  x <- rnorm(30, mean = -100, sd = 10)
  expect_equal(robLoc(x), robscale:::rob_loc_noaitken_impl(x),
               tolerance = 2 * sqrt(.Machine$double.eps))
})

test_that("0.40 Aitken matches plain NR on heavy-tailed data", {
  skip_if_not(
    exists("rob_loc_noaitken_impl", envir = asNamespace("robscale"), mode = "function"),
    "rob_loc_noaitken_impl diagnostic not compiled"
  )
  set.seed(404)
  x <- c(rt(45, df = 3), 50, -60, 100)
  expect_equal(robLoc(x), robscale:::rob_loc_noaitken_impl(x), tolerance = 1e-6)
})

test_that("0.41 Aitken correct at n=4", {
  skip_if_not(
    exists("rob_loc_noaitken_impl", envir = asNamespace("robscale"), mode = "function"),
    "rob_loc_noaitken_impl diagnostic not compiled"
  )
  x <- c(1.0, 2.0, 5.0, 10.0)
  expect_equal(robLoc(x), robscale:::rob_loc_noaitken_impl(x),
               tolerance = 2 * sqrt(.Machine$double.eps))
})

test_that("0.42 n=3 returns median, Aitken not reached", {
  skip_if_not(
    exists("rob_loc_noaitken_impl", envir = asNamespace("robscale"), mode = "function"),
    "rob_loc_noaitken_impl diagnostic not compiled"
  )
  x <- c(1.0, 5.0, 3.0)
  expect_equal(robLoc(x), 3.0, tolerance = 1e-10)
  expect_equal(robscale:::rob_loc_noaitken_impl(x), 3.0, tolerance = 1e-10)
})

test_that("0.43 degenerate scale: all equal, returns that value", {
  skip_if_not(
    exists("rob_loc_noaitken_impl", envir = asNamespace("robscale"), mode = "function"),
    "rob_loc_noaitken_impl diagnostic not compiled"
  )
  x <- rep(7.0, 8)
  expect_equal(robLoc(x), 7.0, tolerance = 1e-10)
  expect_equal(robscale:::rob_loc_noaitken_impl(x), 7.0, tolerance = 1e-10)
})

test_that("0.44 Aitken matches plain NR across 20 seeds, n=64", {
  skip_if_not(
    exists("rob_loc_noaitken_impl", envir = asNamespace("robscale"), mode = "function"),
    "rob_loc_noaitken_impl diagnostic not compiled"
  )
  for (seed in 1:20) {
    set.seed(seed * 7)
    x <- rnorm(64, mean = seed - 10, sd = seed * 0.5 + 1)
    expect_equal(robLoc(x), robscale:::rob_loc_noaitken_impl(x),
                 tolerance = 2 * sqrt(.Machine$double.eps),
                 info = paste("seed =", seed))
  }
})
