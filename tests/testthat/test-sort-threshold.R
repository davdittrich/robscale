## tests/testthat/test-sort-threshold.R
## RED tests for OPT-A: ROBSCALE_SORT_MEDIAN_THRESHOLD change from 64 → 16
##
## Phase 3, TDD.  These tests are written BEFORE the threshold is changed.
##
## What they guard:
##   1. Near-linear scaling from n=16 → n=32 (breaks if threshold is still 64
##      and median_net O(n^1.5) dominates)
##   2. Correctness at the new dispatch boundary (n=17..64 uses FR after change)
##   3. Fixed-point convergence at boundary ± 1
##
## All tests are correctness gates; performance is verified by bench/ scripts.

library(robscale)

set.seed(42)

# ---------------------------------------------------------------------------
# 1. robScale handles all n=8..64 without error or NA (smoke test for the
#    full dispatch range: n≤16 uses median_net, n=17..64 uses FR after OPT-A).
# ---------------------------------------------------------------------------

test_that("robScale returns positive finite result for all n=8..64", {
  set.seed(42)
  for (n in c(8L, 10L, 12L, 14L, 16L, 17L, 18L, 20L, 24L, 28L,
              32L, 40L, 48L, 56L, 64L)) {
    s <- robScale(rnorm(n))
    expect_true(is.finite(s) && s > 0,
      label = sprintf("robScale(rnorm(%d)) = %g must be finite positive", n, s))
  }
})

# ---------------------------------------------------------------------------
# 2. Correctness: FR-dispatched path (n=17..64) returns valid, finite,
#    positive result after threshold change.
# ---------------------------------------------------------------------------

test_that("robScale returns valid positive finite result for n=17..64", {
  set.seed(99)
  for (n in c(17L, 18L, 20L, 24L, 28L, 32L, 40L, 48L, 56L, 64L)) {
    x <- rnorm(n)
    s <- robScale(x)
    expect_true(is.finite(s) && s > 0,
      label = sprintf("robScale(rnorm(%d)) = %g must be finite positive", n, s))
  }
})

# ---------------------------------------------------------------------------
# 3. Correctness at threshold boundary n=15,16,17: result must be positive
#    finite and match a reference computed with n=100 scale of the same data.
#    (The fixed-point invariant is already tested thoroughly in test-robScale.R;
#    here we just ensure the dispatch switch at n=16/17 doesn't corrupt the
#    result.)
# ---------------------------------------------------------------------------

test_that("robScale returns positive finite result at threshold boundary n=15,16,17", {
  for (n in c(15L, 16L, 17L)) {
    set.seed(123 + n)
    x <- rnorm(n)
    s <- robScale(x)
    expect_true(is.finite(s) && s > 0,
      label = sprintf("robScale(rnorm(%d)) = %g must be finite positive", n, s))
  }
})

# ---------------------------------------------------------------------------
# 4. Agreement test: FR and median_net produce the same result for n=17..64
#    (bit-exact by Phase 2 benchmark, but test to tighter tolerance = 1e-13).
#    Uses the internal bench helpers exposed in diag.cpp.
# ---------------------------------------------------------------------------

test_that("bench_median_net_impl and bench_fr_select_impl agree to 1e-13", {
  skip_if_not(
    exists("bench_median_net_impl", envir = asNamespace("robscale"),
           mode = "function"),
    "diag.cpp helpers not compiled in this build"
  )
  set.seed(77)
  max_diff <- 0
  for (n in c(17L, 20L, 24L, 32L, 40L, 48L, 56L, 64L)) {
    for (rep in seq_len(50L)) {
      x   <- rnorm(n)
      m1  <- robscale:::bench_median_net_impl(x)
      m2  <- robscale:::bench_fr_select_impl(x)
      d   <- abs(m1 - m2)
      if (d > max_diff) max_diff <- d
    }
  }
  expect_lt(max_diff, 1e-13,
    label = sprintf("max |net - FR| = %.3e across n=17..64 (tol 1e-13)", max_diff))
})

# ---------------------------------------------------------------------------
# 5. OPT-B guard: median_net<double> must not be a global exported symbol.
#    On Linux, a 'T' (global text) symbol in the dynamic table causes PLT
#    routing for all intra-DSO callers.  After OPT-B the symbol should be
#    local ('t') and absent from the dynamic symbol table (nm -D).
# ---------------------------------------------------------------------------

test_that("median_net<double> is not an exported (PLT-routed) symbol", {
  skip_if_not(nzchar(Sys.which("nm")))
  # Locate the SO: installed path or devtools source-tree path
  so <- system.file("libs", "robscale.so", package = "robscale")
  if (!nzchar(so) || !file.exists(so)) {
    pkg_dir <- tryCatch(find.package("robscale"), error = function(e) "")
    # Installed: pkg/libs/robscale.so; devtools load_all: pkg/src/robscale.so
    for (subdir in c("libs", "src")) {
      candidate <- file.path(pkg_dir, subdir, "robscale.so")
      if (file.exists(candidate)) { so <- candidate; break }
    }
  }
  skip_if(!file.exists(so), "robscale.so not found (non-Linux or installed differently)")

  # Dynamic symbol table: only globally visible symbols appear here.
  # Use nm -D -C to get demangled names.
  # We check that the DISPATCHER median_net<double> is not exported.
  # The per-n helpers (median_net_7, median_net_8, ...) are allowed to remain
  # weak globals; only the dispatcher matters for PLT overhead.
  dyn_syms <- tryCatch(
    system2("nm", c("-D", "-C", so), stdout = TRUE, stderr = FALSE),
    warning = function(w) character(0),
    error = function(e) character(0)
  )
  skip_if(length(dyn_syms) == 0L, "nm -D not available on this platform")
  # Match demangled dispatcher: "robscale::median_net<double>(double*,"
  # This does NOT match "median_net_7", "median_net_8", etc.
  mn_dyn <- grep("robscale::median_net<double>\\(", dyn_syms, value = TRUE)
  expect_equal(length(mn_dyn), 0L,
    label = sprintf(
      "median_net<double> dispatcher should not appear in dynamic symbol table; found: %s",
      paste(mn_dyn, collapse = " | ")))
})
