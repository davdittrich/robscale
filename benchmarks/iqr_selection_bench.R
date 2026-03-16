#!/usr/bin/env Rscript
# iqr_selection_bench.R — Benchmark IQR selection-algorithm variants
#
# Compares 6 C++ variants (in robscale) against stats::IQR and collapse::fquantile.
# Outputs a timestamped CSV in the same directory.

suppressPackageStartupMessages({
  library(robscale)
  library(bench)
  library(collapse)
})

# ---------- constants ---------------------------------------------------------

IQR_K <- 0.741301109252801
TOLERANCE <- 1e-12

SAMPLE_SIZES <- c(1000, 4096, 10000, 50000, 100000, 500000, 1000000, 5000000, 10000000)

DISTRIBUTIONS <- list(
  gaussian = function(n) rnorm(n),
  t3       = function(n) rt(n, df = 3),
  uniform  = function(n) runif(n)
)

# C++ variant functions (internal, not exported in NAMESPACE)
CPP_VARIANTS <- list(
  current            = robscale:::iqr_current,
  fr_incremental     = robscale:::iqr_fr_incremental,
  pdq_dual           = robscale:::iqr_pdq_dual,
  pdq_incremental    = robscale:::iqr_pdq_incremental,
  pdq_branchless_inc = robscale:::iqr_pdq_branchless_inc,
  nth_incremental    = robscale:::iqr_nth_incremental
)

# ---------- correctness gate --------------------------------------------------

cat("=== Correctness gate ===\n")
set.seed(42)
x_test <- rnorm(100000)
ref <- IQR(x_test, type = 7) * IQR_K

for (nm in names(CPP_VARIANTS)) {
  val <- CPP_VARIANTS[[nm]](x_test)
  err <- abs(val - ref)
  status <- if (err < TOLERANCE) "PASS" else sprintf("FAIL (err = %.2e)", err)
  cat(sprintf("  %-25s %s\n", nm, status))
  if (err >= TOLERANCE) stop(sprintf("Variant '%s' failed correctness check", nm))
}

# Also sanity-check collapse
val_coll <- diff(fquantile(x_test, c(0.25, 0.75))) * IQR_K
err_coll <- abs(val_coll - ref)
cat(sprintf("  %-25s %s\n", "collapse_fquantile",
            if (err_coll < TOLERANCE) "PASS" else sprintf("FAIL (err = %.2e)", err_coll)))

cat("All variants passed.\n\n")

# ---------- benchmark ---------------------------------------------------------

cat("=== Benchmark ===\n")

results <- bench::press(
  n = SAMPLE_SIZES,
  dist = names(DISTRIBUTIONS),
  {
    set.seed(123 + n)
    x <- DISTRIBUTIONS[[dist]](n)

    # Adaptive iterations: more for small n, fewer for large n
    iters <- if (n <= 10000) 200L else if (n <= 100000) 50L else if (n <= 1000000) 20L else 5L

    bench::mark(
      current            = robscale:::iqr_current(x),
      fr_incremental     = robscale:::iqr_fr_incremental(x),
      pdq_dual           = robscale:::iqr_pdq_dual(x),
      pdq_incremental    = robscale:::iqr_pdq_incremental(x),
      pdq_branchless_inc = robscale:::iqr_pdq_branchless_inc(x),
      nth_incremental    = robscale:::iqr_nth_incremental(x),
      stats_IQR          = IQR(x, type = 7) * IQR_K,
      collapse_fquantile = diff(fquantile(x, c(0.25, 0.75))) * IQR_K,
      min_iterations = iters,
      check = FALSE   # already verified correctness above
    )
  }
)

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

# Compute speedup vs stats::IQR within each (n, dist) group
res_df$speedup_vs_stats <- NA_real_
for (nn in unique(res_df$n)) {
  for (dd in unique(res_df$dist)) {
    mask <- res_df$n == nn & res_df$dist == dd
    ref_time <- res_df$median[mask & res_df$expression == "stats_IQR"]
    if (length(ref_time) == 1 && !is.na(ref_time) && ref_time > 0) {
      res_df$speedup_vs_stats[mask] <- ref_time / res_df$median[mask]
    }
  }
}

# Platform tag
platform <- if (Sys.info()["sysname"] == "Darwin") "macos" else "linux"
timestamp <- format(Sys.time(), "%Y%m%d")
outfile <- file.path("benchmarks",
                     sprintf("iqr_selection_%s_%s.csv", timestamp, platform))

write.csv(res_df, outfile, row.names = FALSE)
cat(sprintf("\nResults saved to: %s\n", outfile))

# ---------- summary table -----------------------------------------------------

cat("\n=== Summary (median ms, speedup vs stats::IQR) ===\n\n")

# Print compact summary: one row per (variant, n), averaged over distributions
agg <- aggregate(
  cbind(median, speedup_vs_stats) ~ expression + n,
  data = res_df,
  FUN = function(x) round(mean(x, na.rm = TRUE), 3)
)
agg <- agg[order(agg$n, -agg$speedup_vs_stats), ]
print(agg, row.names = FALSE)
