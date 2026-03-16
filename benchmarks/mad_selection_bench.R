#!/usr/bin/env Rscript
# mad_selection_bench.R — Benchmark MAD selection-algorithm variants
#
# Compares 5 C++ variants (in robscale) against stats::mad.
# Outputs a timestamped CSV in the same directory.

suppressPackageStartupMessages({
  library(robscale)
  library(bench)
})

# ---------- constants ---------------------------------------------------------

MAD_CONSISTENCY <- 1.482602218505602
TOLERANCE <- 1e-12

SAMPLE_SIZES <- c(1000, 4096, 10000, 50000, 100000, 500000, 1000000, 5000000, 10000000)

DISTRIBUTIONS <- list(
  gaussian = function(n) rnorm(n),
  t3       = function(n) rt(n, df = 3),
  uniform  = function(n) runif(n)
)

# C++ variant functions (internal, not exported in NAMESPACE)
CPP_VARIANTS <- list(
  current        = robscale:::mad_bench_current,
  pdq            = robscale:::mad_bench_pdq,
  pdq_branchless = robscale:::mad_bench_pdq_branchless,
  nth            = robscale:::mad_bench_nth,
  miniselect_fr  = robscale:::mad_bench_miniselect_fr
)

# Check for collapse::fmad availability
has_collapse_fmad <- requireNamespace("collapse", quietly = TRUE) &&
  exists("fmad", where = asNamespace("collapse"), inherits = FALSE)

if (has_collapse_fmad) {
  cat("collapse::fmad detected — including in benchmark.\n")
} else {
  cat("collapse::fmad not available — skipping.\n")
}

# ---------- correctness gate --------------------------------------------------

cat("=== Correctness gate ===\n")
set.seed(42)
x_test <- rnorm(100000)
ref <- stats::mad(x_test, constant = MAD_CONSISTENCY)

for (nm in names(CPP_VARIANTS)) {
  val <- CPP_VARIANTS[[nm]](x_test)
  err <- abs(val - ref)
  status <- if (err < TOLERANCE) "PASS" else sprintf("FAIL (err = %.2e)", err)
  cat(sprintf("  %-25s %s\n", nm, status))
  if (err >= TOLERANCE) stop(sprintf("Variant '%s' failed correctness check", nm))
}

if (has_collapse_fmad) {
  val_coll <- collapse::fmad(x_test)
  err_coll <- abs(val_coll - ref)
  cat(sprintf("  %-25s %s\n", "collapse_fmad",
              if (err_coll < TOLERANCE) "PASS" else sprintf("FAIL (err = %.2e)", err_coll)))
}

cat("All variants passed.\n\n")

# ---------- benchmark ---------------------------------------------------------

cat("=== Benchmark ===\n")

if (has_collapse_fmad) {
  results <- bench::press(
    n = SAMPLE_SIZES,
    dist = names(DISTRIBUTIONS),
    {
      set.seed(123 + n)
      x <- DISTRIBUTIONS[[dist]](n)
      iters <- if (n <= 10000) 200L else if (n <= 100000) 50L else if (n <= 1000000) 20L else 5L

      bench::mark(
        current        = robscale:::mad_bench_current(x),
        pdq            = robscale:::mad_bench_pdq(x),
        pdq_branchless = robscale:::mad_bench_pdq_branchless(x),
        nth            = robscale:::mad_bench_nth(x),
        miniselect_fr  = robscale:::mad_bench_miniselect_fr(x),
        stats_mad      = stats::mad(x, constant = MAD_CONSISTENCY),
        collapse_fmad  = collapse::fmad(x),
        min_iterations = iters,
        check = FALSE
      )
    }
  )
} else {
  results <- bench::press(
    n = SAMPLE_SIZES,
    dist = names(DISTRIBUTIONS),
    {
      set.seed(123 + n)
      x <- DISTRIBUTIONS[[dist]](n)
      iters <- if (n <= 10000) 200L else if (n <= 100000) 50L else if (n <= 1000000) 20L else 5L

      bench::mark(
        current        = robscale:::mad_bench_current(x),
        pdq            = robscale:::mad_bench_pdq(x),
        pdq_branchless = robscale:::mad_bench_pdq_branchless(x),
        nth            = robscale:::mad_bench_nth(x),
        miniselect_fr  = robscale:::mad_bench_miniselect_fr(x),
        stats_mad      = stats::mad(x, constant = MAD_CONSISTENCY),
        min_iterations = iters,
        check = FALSE
      )
    }
  )
}

# ---------- post-process & save -----------------------------------------------

res_df <- as.data.frame(results)
# expression column from bench is a list — coerce to character
res_df$expression <- as.character(res_df$expression)

# Keep useful columns
keep <- c("expression", "n", "dist", "min", "median", "mem_alloc", "n_itr", "n_gc")
keep <- intersect(keep, names(res_df))
res_df <- res_df[, keep]

# Convert bench_time to numeric milliseconds
for (col in c("min", "median")) {
  if (col %in% names(res_df)) {
    res_df[[col]] <- as.numeric(res_df[[col]], units = "ms")
  }
}

# Compute speedup vs stats::mad within each (n, dist) group
res_df$speedup_vs_stats <- NA_real_
for (nn in unique(res_df$n)) {
  for (dd in unique(res_df$dist)) {
    mask <- res_df$n == nn & res_df$dist == dd
    ref_time <- res_df$median[mask & res_df$expression == "stats_mad"]
    if (length(ref_time) == 1 && !is.na(ref_time) && ref_time > 0) {
      res_df$speedup_vs_stats[mask] <- ref_time / res_df$median[mask]
    }
  }
}

# Platform tag
platform <- if (Sys.info()["sysname"] == "Darwin") "macos" else "linux"
timestamp <- format(Sys.time(), "%Y%m%d")
outfile <- file.path("benchmarks",
                     sprintf("mad_selection_%s_%s.csv", timestamp, platform))

write.csv(res_df, outfile, row.names = FALSE)
cat(sprintf("\nResults saved to: %s\n", outfile))

# ---------- summary table -----------------------------------------------------

cat("\n=== Summary (median ms, speedup vs stats::mad) ===\n\n")

# Print compact summary: one row per (variant, n), averaged over distributions
agg <- aggregate(
  cbind(median, speedup_vs_stats) ~ expression + n,
  data = res_df,
  FUN = function(x) round(mean(x, na.rm = TRUE), 3)
)
agg <- agg[order(agg$n, -agg$speedup_vs_stats), ]
print(agg, row.names = FALSE)
