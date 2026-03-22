library(testthat)
library(robscale)

# Helper: high-precision M-scale reference (tol=1e-14, maxit=10000).
# Mirrors rob_scale_core: median → MAD → fixed-point iteration.
rob_scale_ref <- function(x, tol = 1e-14, maxit = 10000L) {
  C     <- 0.37394112142347236
  MAD_C <- 1.482602218505602
  center <- median(x)
  s <- MAD_C * median(abs(x - center))
  if (s == 0) return(0)
  for (k in seq_len(maxit)) {
    u <- (x - center) * (0.5 / (C * s))
    v <- sqrt(2 * mean(tanh(u)^2))
    s <- s * v
    if (abs(v - 1) <= tol) break
  }
  s
}

# Oscillating-case inputs, identified by diagnostic analysis:
#   n=5, seed=47:   rho_evals=50  (Aitken fires only once, then oscillates)
#   n=5, seed=247:  rho_evals=71  (Aitken never fires)
#   n=7, seed=49:   rho_evals=80, converged=FALSE  ← maxit hit
#   n=7, seed=149:  rho_evals=55  (Aitken never fires)

make_x <- function(n, seed) { set.seed(seed); rnorm(n) }

tol2 <- 2 * sqrt(.Machine$double.eps)  # tolerance matching pinning test

# ── Convergence tests ──────────────────────────────────────────────────────

test_that("robScale: n=7 seed=49 converges (was: maxit hit, converged=FALSE)", {
  x <- make_x(7L, 49L)
  diag <- robscale:::rob_scale_diag_impl(x)
  expect_true(diag$converged,
    label = "n=7 seed=49 must reach |v-1|<=tol within maxit iterations")
})

test_that("robScale: n=7 seed=149 converges (was: rho_evals=55, slow)", {
  x <- make_x(7L, 149L)
  diag <- robscale:::rob_scale_diag_impl(x)
  expect_true(diag$converged, label = "n=7 seed=149 must converge")
})

# ── Accuracy tests ─────────────────────────────────────────────────────────

test_that("robScale: n=7 seed=49 result within 2*sqrt(eps) of reference", {
  x     <- make_x(7L, 49L)
  s_imp <- robscale::robScale(x)
  s_ref <- rob_scale_ref(x)
  expect_equal(s_imp, s_ref, tolerance = tol2,
    label = "n=7 seed=49: implementation must match high-precision reference")
})

test_that("robScale: n=5 seed=47 result within 2*sqrt(eps) of reference", {
  x     <- make_x(5L, 47L)
  s_imp <- robscale::robScale(x)
  s_ref <- rob_scale_ref(x)
  expect_equal(s_imp, s_ref, tolerance = tol2,
    label = "n=5 seed=47: implementation must match high-precision reference")
})

test_that("robScale: n=5 seed=247 result within 2*sqrt(eps) of reference", {
  x     <- make_x(5L, 247L)
  s_imp <- robscale::robScale(x)
  s_ref <- rob_scale_ref(x)
  expect_equal(s_imp, s_ref, tolerance = tol2,
    label = "n=5 seed=247: implementation must match high-precision reference")
})

# ── Iteration-count tests (rho_evals) ─────────────────────────────────────
# Thresholds are generous: fast path uses 5-11, oscillation fix should use <30.
# (Current: 50, 71, 80, 55 respectively — all must fall below threshold.)

test_that("robScale: n=5 seed=47 rho_evals reduced (was 50)", {
  x    <- make_x(5L, 47L)
  diag <- robscale:::rob_scale_diag_impl(x)
  expect_lt(diag$rho_evals, 30L,
    label = "n=5 seed=47: geometric-mean step must cut rho_evals below 30")
})

test_that("robScale: n=5 seed=247 rho_evals reduced (was 71)", {
  x    <- make_x(5L, 247L)
  diag <- robscale:::rob_scale_diag_impl(x)
  expect_lt(diag$rho_evals, 30L,
    label = "n=5 seed=247: geometric-mean step must cut rho_evals below 30")
})

test_that("robScale: n=7 seed=49 rho_evals reduced (was 80, non-convergent)", {
  x    <- make_x(7L, 49L)
  diag <- robscale:::rob_scale_diag_impl(x)
  expect_lt(diag$rho_evals, 40L,
    label = "n=7 seed=49: geometric-mean step must cut rho_evals below 40")
})

test_that("robScale: n=7 seed=149 rho_evals reduced (was 55)", {
  x    <- make_x(7L, 149L)
  diag <- robscale:::rob_scale_diag_impl(x)
  expect_lt(diag$rho_evals, 30L,
    label = "n=7 seed=149: geometric-mean step must cut rho_evals below 30")
})

# ── Non-regression: previously-fast cases must stay fast ──────────────────

test_that("robScale: n=5 seed=147 still fast (was rho_evals=7)", {
  x    <- make_x(5L, 147L)
  diag <- robscale:::rob_scale_diag_impl(x)
  expect_lt(diag$rho_evals, 15L,
    label = "n=5 seed=147 was already fast; must not regress")
})

test_that("robScale: n=7 seed=249 still fast (was rho_evals=7)", {
  x    <- make_x(7L, 249L)
  diag <- robscale:::rob_scale_diag_impl(x)
  expect_lt(diag$rho_evals, 15L,
    label = "n=7 seed=249 was already fast; must not regress")
})

test_that("robScale: n=4 typical cases still fast", {
  for (seed in c(46L, 146L, 246L)) {
    x    <- make_x(4L, seed)
    diag <- robscale:::rob_scale_diag_impl(x)
    expect_lt(diag$rho_evals, 20L,
      label = paste0("n=4 seed=", seed, ": must remain fast"))
  }
})
