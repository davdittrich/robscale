#!/usr/bin/env Rscript
# select_crossover_sweep.R
#
# Directly times Floyd-Rivest vs pdqselect for each selection type used
# in the package, sweeping n from 1K to 100M. Finds empirical crossover n.
#
# Requires the package to be built with src/select_bench.cpp (exposes
# sel_fr_*, sel_pdq_* functions that copy input to a heap buffer and run
# selection once — no surrounding estimator overhead).
#
# Usage:
#   Rscript --vanilla benchmarks/select_crossover_sweep.R
#
# Output: CSV written to benchmarks/select_crossover_<date>_<platform>.csv

suppressPackageStartupMessages(library(robscale))
suppressPackageStartupMessages(library(bench))

platform <- ifelse(Sys.info()[["sysname"]] == "Darwin", "macos", "linux")
cat(sprintf("Platform: %s\n", platform))

# Three distributions per call site, to check if crossover is distribution-sensitive.
make_x <- function(n, dist) {
  set.seed(42)
  switch(dist,
    gaussian   = rnorm(n),
    uniform    = runif(n),
    contam     = c(rnorm(round(0.9 * n)), rnorm(round(0.1 * n), mean = 10, sd = 0.1))
  )
}

# n grid: fine-grained around expected crossover (~10K–100K), coarser elsewhere.
sizes <- c(
  1000L, 2000L, 5000L,
  10000L, 15000L, 20000L, 30000L, 50000L,
  75000L, 100000L, 200000L, 500000L,
  1000000L, 5000000L, 10000000L, 50000000L, 100000000L
)

dists   <- c("gaussian", "uniform", "contam")
n_iter  <- 30L   # bench::mark iterations per cell

rows <- list()
k    <- 1L

fr_med  <- robscale:::sel_fr_median
pdq_med <- robscale:::sel_pdq_median
fr_low  <- robscale:::sel_fr_lowmedian
pdq_low <- robscale:::sel_pdq_lowmedian
fr_kth  <- robscale:::sel_fr_kth
pdq_kth <- robscale:::sel_pdq_kth

for (dist in dists) {
  cat(sprintf("\nDistribution: %s\n", dist))
  for (n in sizes) {
    x <- make_x(n, dist)
    k_idx <- as.integer(n / 4L)   # 25th percentile — arbitrary k for Qn proxy

    # Median selection
    bm_fr  <- bench::mark(fr_med(x),  min_iterations = n_iter, memory = FALSE, check = FALSE)
    bm_pdq <- bench::mark(pdq_med(x), min_iterations = n_iter, memory = FALSE, check = FALSE)
    rows[[k]] <- data.frame(
      sel_type = "median", dist = dist, n = n,
      fr_us  = as.numeric(bm_fr$median)  * 1e6,
      pdq_us = as.numeric(bm_pdq$median) * 1e6
    )
    k <- k + 1L

    # Low-median selection
    bm_fr  <- bench::mark(fr_low(x),  min_iterations = n_iter, memory = FALSE, check = FALSE)
    bm_pdq <- bench::mark(pdq_low(x), min_iterations = n_iter, memory = FALSE, check = FALSE)
    rows[[k]] <- data.frame(
      sel_type = "lowmedian", dist = dist, n = n,
      fr_us  = as.numeric(bm_fr$median)  * 1e6,
      pdq_us = as.numeric(bm_pdq$median) * 1e6
    )
    k <- k + 1L

    # Arbitrary-k selection (k = n/4)
    bm_fr  <- bench::mark(fr_kth(x, k_idx),  min_iterations = n_iter, memory = FALSE, check = FALSE)
    bm_pdq <- bench::mark(pdq_kth(x, k_idx), min_iterations = n_iter, memory = FALSE, check = FALSE)
    rows[[k]] <- data.frame(
      sel_type = "kth", dist = dist, n = n,
      fr_us  = as.numeric(bm_fr$median)  * 1e6,
      pdq_us = as.numeric(bm_pdq$median) * 1e6
    )
    k <- k + 1L

    cat(sprintf(
      "  n=%9d | med FR=%7.1f pdq=%7.1f (%.2fx) | low FR=%7.1f pdq=%7.1f (%.2fx) | kth FR=%7.1f pdq=%7.1f (%.2fx)\n",
      n,
      rows[[k-3]]$fr_us, rows[[k-3]]$pdq_us, rows[[k-3]]$fr_us / rows[[k-3]]$pdq_us,
      rows[[k-2]]$fr_us, rows[[k-2]]$pdq_us, rows[[k-2]]$fr_us / rows[[k-2]]$pdq_us,
      rows[[k-1]]$fr_us, rows[[k-1]]$pdq_us, rows[[k-1]]$fr_us / rows[[k-1]]$pdq_us
    ))
  }
}

df <- do.call(rbind, rows)
df$speedup <- df$fr_us / df$pdq_us

# ---------------------------------------------------------------------------
# Summary: crossover n per sel_type (first n where pdqselect is faster)
# ---------------------------------------------------------------------------
cat("\n\n=== Crossover summary (first n where speedup > 1.00 for gaussian dist) ===\n")
for (st in c("median", "lowmedian", "kth")) {
  sub <- df[df$sel_type == st & df$dist == "gaussian", ]
  sub <- sub[order(sub$n), ]
  cross <- sub$n[which(sub$speedup > 1.0)[1]]
  if (is.na(cross)) cross <- Inf
  cat(sprintf("  %-10s  crossover n ~ %d\n", st, cross))
}

# ---------------------------------------------------------------------------
# Back-derive divisors from empirical crossover
# ---------------------------------------------------------------------------
cfg <- robscale:::get_qnsn_config()
l2  <- as.numeric(cfg$l2_per_core)
cat(sprintf("\nL2 per core: %d bytes (%.0f KB)\n", as.integer(l2), l2 / 1024))
cat("\n=== Implied divisors (divisor = L2 / (8 * crossover_n)) ===\n")
for (st in c("median", "lowmedian", "kth")) {
  sub <- df[df$sel_type == st & df$dist == "gaussian", ]
  sub <- sub[order(sub$n), ]
  cross <- sub$n[which(sub$speedup > 1.0)[1]]
  if (!is.na(cross) && is.finite(cross)) {
    div <- l2 / (8 * cross)
    cat(sprintf("  %-10s  crossover n=%d  implied divisor=%.2f\n", st, cross, div))
  } else {
    cat(sprintf("  %-10s  no crossover found in range\n", st))
  }
}

# ---------------------------------------------------------------------------
# Current thresholds for comparison
# ---------------------------------------------------------------------------
cat(sprintf("\n=== Current thresholds ===\n"))
cat(sprintf("  pdq_median_threshold     = %d\n", cfg$pdq_median_threshold))
cat(sprintf("  pdq_robscale_threshold   = %d\n", cfg$pdq_robscale_threshold))
cat(sprintf("  pdq_lowmedian_threshold  = %d\n", cfg$pdq_lowmedian_threshold))
cat(sprintf("  pdq_qn_final_threshold   = %d\n", cfg$pdq_qn_final_threshold))

out_file <- sprintf("benchmarks/select_crossover_%s_%s.csv",
                    format(Sys.Date(), "%Y%m%d"), platform)
write.csv(df, out_file, row.names = FALSE)
cat(sprintf("\nSaved to %s\n", out_file))
