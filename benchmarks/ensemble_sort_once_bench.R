#!/usr/bin/env Rscript
# Benchmark: Sort-Once Optimization for Ensemble Bootstrap
#
# Measures the potential gain from pre-sorting each resample once
# instead of letting GMD, Sn, Qn each sort independently.
#
# Requires: package built with -DROBSCALE_BENCH_INTERNALS
# Usage:
#   PKG_CPPFLAGS="-DROBSCALE_BENCH_INTERNALS" Rscript -e 'devtools::load_all()'; \
#     Rscript benchmarks/ensemble_sort_once_bench.R
#   — or —
#   R CMD INSTALL --configure-args="" . && Rscript benchmarks/ensemble_sort_once_bench.R
#   (after adding -DROBSCALE_BENCH_INTERNALS to Makevars)

# --- Load package and extract harness functions ---
# When run via devtools::load_all(), functions are in globalenv.
# When run via R CMD INSTALL, they're in the robscale namespace (unexported).
if (!exists("bench_ensemble_internals", envir = globalenv(), inherits = FALSE)) {
  if (requireNamespace("robscale", quietly = TRUE) &&
      exists("bench_ensemble_internals", where = asNamespace("robscale"),
             inherits = FALSE)) {
    bench_ensemble_internals <- robscale:::bench_ensemble_internals
    bench_sort_cost_us       <- robscale:::bench_sort_cost_us
    bench_estimator_breakdown <- robscale:::bench_estimator_breakdown
  } else {
    stop(
      "Benchmark harness not found.\n",
      "Rebuild with: PKG_CPPFLAGS=\"-DROBSCALE_BENCH_INTERNALS\" ",
      "Rscript -e 'devtools::load_all()'",
      call. = FALSE
    )
  }
}

# --- Platform info ---
platform <- if (Sys.info()[["sysname"]] == "Darwin") {
  arch <- Sys.info()[["machine"]]
  paste0("macOS-", arch)
} else {
  arch <- Sys.info()[["machine"]]
  paste0("Linux-", arch)
}
cat(sprintf("\nEnsemble Sort-Once Benchmark — %s\n", platform))
cat(sprintf("Date: %s\n\n", Sys.Date()))

# --- Dense n grid (55 points) ---
n_grid <- sort(unique(c(
  2:16,                          # every integer in sorting-network range
  seq(18, 32, by = 2),           # transition zone
  seq(36, 64, by = 4),           # Qn brute-force threshold (~64)
  seq(72, 128, by = 8),          # medium
  seq(144, 256, by = 16),        # larger
  seq(288, 512, by = 32)         # target ceiling
)))
cat(sprintf("Grid: %d points, n = %d..%d\n", length(n_grid), min(n_grid), max(n_grid)))

# --- Fixed parameters ---
N_BOOT  <- 200L
REPS_AB <- 101L   # odd for clean median (ensemble A/B)
REPS_SORT <- 1001L  # higher reps for micro-bench (sort cost)
REPS_EST  <- 101L   # estimator breakdown

# --- Preallocate results ---
results <- data.frame(
  platform        = character(length(n_grid)),
  n               = integer(length(n_grid)),
  n_boot          = integer(length(n_grid)),
  sort_us         = numeric(length(n_grid)),
  sd_c4_us        = numeric(length(n_grid)),
  gmd_us          = numeric(length(n_grid)),
  mad_us          = numeric(length(n_grid)),
  iqr_us          = numeric(length(n_grid)),
  sn_us           = numeric(length(n_grid)),
  qn_us           = numeric(length(n_grid)),
  rob_scale_us    = numeric(length(n_grid)),
  current_total_us  = numeric(length(n_grid)),
  presort_total_us  = numeric(length(n_grid)),
  speedup         = numeric(length(n_grid)),
  sort_fraction   = numeric(length(n_grid)),
  pct_in_sn_qn    = numeric(length(n_grid)),
  stringsAsFactors = FALSE
)

# --- Progress bar ---
cat(sprintf("\nRunning %d grid points (n_boot=%d, reps=%d/%d/%d)...\n\n",
            length(n_grid), N_BOOT, REPS_SORT, REPS_EST, REPS_AB))

header_fmt <- "%5s  %7s  %10s  %10s  %7s  %9s  %6s  %s\n"
row_fmt    <- "%5d  %7.2f  %10.1f  %10.1f  %6.2fx  %8.1f%%  %5.1f%%  %s\n"
cat(sprintf(header_fmt, "n", "sort_us", "current_us", "presort_us",
            "speedup", "sort_frac", "sn+qn", "verdict"))
cat(paste(rep("-", 78), collapse = ""), "\n")

for (i in seq_along(n_grid)) {
  n <- n_grid[i]

  # Generate standard normal data
  set.seed(42L)
  x <- rnorm(n)

  # 1. Sort cost isolation
  sort_us <- bench_sort_cost_us(x, REPS_SORT)

  # 2. Estimator breakdown (7-element named vector, median µs each)
  est <- bench_estimator_breakdown(x, REPS_EST)

  # 3. Full ensemble A/B
  ab <- bench_ensemble_internals(x, N_BOOT, REPS_AB)

  # Extract medians from A/B results
  current_rows <- ab[ab$mode == "current", ]
  presort_rows <- ab[ab$mode == "presort", ]
  current_total_us <- median(current_rows$total_us)
  presort_total_us <- median(presort_rows$total_us)

  # Derived metrics
  speedup <- current_total_us / presort_total_us
  sort_fraction <- 3.0 * sort_us * N_BOOT / current_total_us * 100  # percentage
  est_total <- sum(est)
  pct_sn_qn <- if (est_total > 0) (est["sn"] + est["qn"]) / est_total * 100 else 0

  # Verdict
  verdict <- if (speedup < 1.05) {
    "SKIP"
  } else if (speedup < 1.15) {
    "MARGINAL"
  } else {
    "WORTH IT"
  }

  # Store results
  results$platform[i]        <- platform
  results$n[i]               <- n
  results$n_boot[i]          <- N_BOOT
  results$sort_us[i]         <- sort_us
  results$sd_c4_us[i]        <- est["sd_c4"]
  results$gmd_us[i]          <- est["gmd"]
  results$mad_us[i]          <- est["mad"]
  results$iqr_us[i]          <- est["iqr"]
  results$sn_us[i]           <- est["sn"]
  results$qn_us[i]           <- est["qn"]
  results$rob_scale_us[i]    <- est["rob_scale"]
  results$current_total_us[i]  <- current_total_us
  results$presort_total_us[i]  <- presort_total_us
  results$speedup[i]         <- speedup
  results$sort_fraction[i]   <- sort_fraction
  results$pct_in_sn_qn[i]    <- pct_sn_qn

  # Print row
  cat(sprintf(row_fmt, n, sort_us, current_total_us, presort_total_us,
              speedup, sort_fraction, pct_sn_qn, verdict))
}

# --- Save CSV ---
outfile <- sprintf("benchmarks/ensemble_sort_once_%s.csv", format(Sys.Date(), "%Y%m%d"))
# If run from package root, save there; if not, save to /tmp
if (file.exists("benchmarks")) {
  write.csv(results, outfile, row.names = FALSE)
} else {
  outfile <- sprintf("/tmp/ensemble_sort_once_%s.csv", format(Sys.Date(), "%Y%m%d"))
  write.csv(results, outfile, row.names = FALSE)
}
cat(sprintf("\nResults saved to: %s\n", outfile))

# --- Summary ---
cat("\n=== Summary ===\n")
worth_it <- results[results$speedup >= 1.15, ]
marginal <- results[results$speedup >= 1.05 & results$speedup < 1.15, ]
skip     <- results[results$speedup < 1.05, ]

if (nrow(worth_it) > 0) {
  cat(sprintf("WORTH IT (>=1.15x): n >= %d (min speedup %.2fx, max %.2fx)\n",
              min(worth_it$n), min(worth_it$speedup), max(worth_it$speedup)))
} else {
  cat("WORTH IT: none\n")
}
if (nrow(marginal) > 0) {
  cat(sprintf("MARGINAL (1.05-1.15x): n in [%d, %d]\n",
              min(marginal$n), max(marginal$n)))
} else {
  cat("MARGINAL: none\n")
}
cat(sprintf("SKIP (<1.05x): %d grid points\n", nrow(skip)))

# Sn+Qn dominance at crossover
if (nrow(worth_it) > 0) {
  crossover_n <- min(worth_it$n)
  crossover_row <- results[results$n == crossover_n, ]
  cat(sprintf("\nAt crossover n=%d: sn+qn = %.1f%% of single-replicate time\n",
              crossover_n, crossover_row$pct_in_sn_qn))
  if (crossover_row$pct_in_sn_qn > 40) {
    cat("  -> Sn/Qn _sorted variants would yield further gains\n")
  } else {
    cat("  -> Sn/Qn _sorted variants likely NOT needed\n")
  }
}

cat("\nDone.\n")
