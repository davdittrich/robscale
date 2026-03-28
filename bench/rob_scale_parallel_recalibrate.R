# Parallel threshold recalibration: measure per-element NR time across sizes.
#
# The production robScale() uses serial NR below rob_scale_parallel_threshold
# and TBB parallel NR above it. This script measures wall time at sizes
# straddling the threshold to find the actual crossover.
#
# Per-element time should DECREASE at the crossover if parallel is beneficial.
# If per-element time INCREASES at the crossover, the threshold is too low.
#
# Usage:
#   Rscript bench/rob_scale_parallel_recalibrate.R        # auto run number
#   Rscript bench/rob_scale_parallel_recalibrate.R 1      # explicit run

library(robscale)
library(bench)

# Derive the current threshold from RuntimeConfig (not directly exposed to R)
cfg <- robscale:::get_qnsn_config()
l2_per_core <- cfg$l2_per_core
# Formula: max(4096, l2_per_core / (sizeof(double) * 4)) = max(4096, l2_per_core / 32)
threshold <- max(4096L, as.integer(l2_per_core / 32))
cat("L2 per core:", l2_per_core, "bytes\n")
cat("Derived rob_scale_parallel_threshold:", threshold, "\n\n")

# Sweep sizes: below, at, and above the threshold
sizes <- sort(unique(c(
  2048L, 4096L,
  as.integer(threshold / 2),
  as.integer(threshold * 0.75),
  as.integer(threshold - 1L),    # last serial
  as.integer(threshold),          # first parallel
  as.integer(threshold + 1L),
  as.integer(threshold * 1.5),
  as.integer(threshold * 2),
  as.integer(threshold * 4),
  65536L
)))

min_iter_for <- function(n) {
  if (n <= 4096L)  500L
  else if (n <= 16384L) 200L
  else if (n <= 32768L) 100L
  else 50L
}

args <- commandArgs(trailingOnly = TRUE)
run_n <- if (length(args) >= 1L) as.integer(args[1L]) else {
  existing <- sum(file.exists(file.path("bench", paste0("par_recal_run", 1:3, ".rds"))))
  min(existing + 1L, 3L)
}

cat("=== Parallel Threshold Recalibration — Run", run_n, "of 3 ===\n\n")

results <- list()
for (n in sizes) {
  set.seed(42L + n)
  x <- rnorm(n)
  mi <- min_iter_for(n)

  bm <- bench::mark(
    robScale = robScale(x),
    min_iterations = mi,
    max_time = 60,
    check = FALSE
  )

  med_ns <- as.numeric(bm$median[[1L]]) * 1e9
  per_elem_ns <- med_ns / n
  dispatch <- if (n >= threshold) "PARALLEL" else "serial"

  results[[length(results) + 1L]] <- list(
    n = n, med_ns = med_ns, per_elem_ns = per_elem_ns, dispatch = dispatch
  )
}

# Print results
cat(sprintf("%-10s  %10s  %10s  %s\n", "n", "total(ns)", "per_elem", "dispatch"))
cat(strrep("-", 50), "\n")
for (r in results) {
  cat(sprintf("%-10d  %10.1f  %10.2f  %s\n",
      r$n, r$med_ns, r$per_elem_ns, r$dispatch))
}

# Find the crossover: compare per-element time at last serial vs first parallel
serial_results <- Filter(function(r) r$dispatch == "serial", results)
parallel_results <- Filter(function(r) r$dispatch == "PARALLEL", results)

if (length(serial_results) > 0 && length(parallel_results) > 0) {
  last_serial <- serial_results[[length(serial_results)]]
  first_parallel <- parallel_results[[1L]]

  ratio <- first_parallel$per_elem_ns / last_serial$per_elem_ns
  cat(sprintf("\nCrossover analysis:\n"))
  cat(sprintf("  Last serial  (n=%d): %.2f ns/elem\n", last_serial$n, last_serial$per_elem_ns))
  cat(sprintf("  First parallel (n=%d): %.2f ns/elem\n", first_parallel$n, first_parallel$per_elem_ns))
  cat(sprintf("  Ratio (parallel/serial per-elem): %.4f\n", ratio))

  if (ratio < 0.95) {
    cat("  >>> Parallel is FASTER at threshold — threshold is CORRECT or could be LOWERED\n")
  } else if (ratio <= 1.05) {
    cat("  >>> Parallel ≈ serial at threshold — threshold is in the right ballpark\n")
  } else {
    cat("  >>> Parallel is SLOWER at threshold — threshold should be RAISED\n")

    # Find the actual crossover: smallest n where parallel per-elem < serial per-elem
    # Use the last serial per-elem as reference
    ref_per_elem <- last_serial$per_elem_ns
    crossover_n <- NA
    for (r in parallel_results) {
      if (r$per_elem_ns < ref_per_elem * 0.95) {
        crossover_n <- r$n
        break
      }
    }
    if (!is.na(crossover_n)) {
      cat(sprintf("  >>> Actual crossover at n ≈ %d\n", crossover_n))
    } else {
      cat("  >>> No clear crossover found in measured range — serial wins everywhere\n")
    }
  }
}

outfile <- file.path("bench", paste0("par_recal_run", run_n, ".rds"))
saveRDS(results, outfile)
cat("\nSaved:", outfile, "\n")
