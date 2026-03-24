## bench/adm_parallel_threshold.R
## WU-ADM4: calibrate adm_parallel_threshold.
## H2H sweep: serial range (below threshold) vs parallel range (above threshold).
## Decision criterion: threshold = smallest n where parallel ratio < 0.90 (10% gain).
##
## Usage: Rscript bench/adm_parallel_threshold.R
##
## Interpretation:
##   µs/1k column: normalised throughput (lower = faster per element).
##   Serial range: TBB inactive (n < adm_parallel_threshold).
##   Parallel range: TBB active (n >= adm_parallel_threshold).
##   Regression = parallel µs/1k higher than serial µs/1k.

library(robscale)
library(bench)

## Threshold computed from hardware: max(65536, per_core_l2 * 4 / sizeof(double))
## On this machine (512KB L2, x86_64): max(65536, 262144) = 262144
THRESHOLD <- 262144L
WARMUP    <- 5L

cat("=== adm_parallel_threshold calibration ===\n")
cat(sprintf("Threshold on this machine: %d\n\n", THRESHOLD))

sweep <- function(n, min_iter) {
  set.seed(42L); x <- rnorm(n)
  for (i in seq_len(WARMUP)) adm(x)
  bm <- bench::mark(adm(x), min_iterations = min_iter, check = FALSE)
  med_us <- as.numeric(bm$median) * 1e6
  list(n = n, med_us = med_us, us_per_1k = med_us / n * 1000,
       parallel = (n >= THRESHOLD))
}

## Serial range (below threshold)
serial_sizes <- c(32768L, 65536L, 131072L, 196608L, 245760L)
## Parallel range (above threshold)
par_sizes    <- c(262144L, 393216L, 524288L, 786432L, 1048576L, 2097152L)

cat(sprintf("%-10s  %10s  %10s  %s\n", "n", "median(µs)", "µs/1k", "path"))
cat(strrep("-", 50), "\n")

all_results <- list()
for (n in c(serial_sizes, par_sizes)) {
  min_iter <- if (n <= 131072L) 50L else if (n <= 524288L) 20L else 10L
  r <- sweep(n, min_iter)
  all_results[[length(all_results) + 1L]] <- r
  cat(sprintf("%-10d  %10.1f  %10.3f  %s\n",
              r$n, r$med_us, r$us_per_1k,
              if (r$parallel) "PARALLEL (TBB)" else "serial"))
}

## Decision
cat(strrep("-", 50), "\n")

serial_us_1k <- mean(sapply(all_results[!sapply(all_results, `[[`, "parallel")],
                            `[[`, "us_per_1k"))
par_us_1k    <- mean(sapply(all_results[sapply(all_results, `[[`, "parallel")],
                            `[[`, "us_per_1k"))

cat(sprintf("\nSerial avg µs/1k:   %.3f\n", serial_us_1k))
cat(sprintf("Parallel avg µs/1k: %.3f\n", par_us_1k))
cat(sprintf("Parallel / Serial:  %.3f\n", par_us_1k / serial_us_1k))
cat(sprintf("\nDecision criterion: parallel ratio < 0.90 (10%% gain)\n"))

if (par_us_1k / serial_us_1k < 0.90) {
  cat(sprintf("RESULT: TBB provides speedup. Threshold %d is appropriate.\n", THRESHOLD))
} else {
  cat(sprintf("RESULT: TBB does NOT provide speedup (ratio %.3f >= 0.90).\n",
              par_us_1k / serial_us_1k))
  cat("ACTION: adm_parallel_threshold should be set to SIZE_MAX (disable TBB for adm).\n")
}
