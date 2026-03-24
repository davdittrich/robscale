## bench/ensemble_gate_perf1.R
## Timing gate for WU-PERF-1 (hoist TLS singleton from ensemble_one_replicate)
##
## Usage:
##   Rscript bench/ensemble_gate_perf1.R          # CAPTURE baseline (BEFORE change)
##   Rscript bench/ensemble_gate_perf1.R CHECK    # CHECK post-change vs baseline

library(robscale)
library(bench)

RATIO_LIMIT <- 1.05
BASELINE_FILE <- "bench/ensemble_perf1_pre.rds"
MODE <- if (length(commandArgs(trailingOnly = TRUE)) > 0 &&
            toupper(commandArgs(trailingOnly = TRUE)[1]) == "CHECK") "CHECK" else "CAPTURE"

stable_sizes <- c(1000L, 2000L, 5000L)
noisy_sizes  <- c(5L, 10L, 15L, 19L)
N_BOOT <- 500L

cat("=== ensemble_gate_perf1 ===\n")
cat("Mode:", MODE, "\n")

if (MODE == "CAPTURE") {
  all_sizes <- c(noisy_sizes, stable_sizes)
  pre <- lapply(all_sizes, function(n) {
    set.seed(42); x <- rnorm(n)
    bm <- bench::mark(
      robscale:::cpp_scale_ensemble(x, N_BOOT),
      min_iterations = 5L, check = FALSE
    )
    med_ns <- as.numeric(bm$median[1]) * 1e9
    cat(sprintf("  n=%5d: %.0f ns\n", n, med_ns))
    list(n = n, median_ns = med_ns)
  })
  saveRDS(pre, BASELINE_FILE)
  cat("Baseline saved to", BASELINE_FILE, "\n")
} else {
  if (!file.exists(BASELINE_FILE)) { cat("ERROR: baseline not found.\n"); quit(status = 1) }
  pre <- readRDS(BASELINE_FILE)
  stable_fails <- 0L
  for (b in pre) {
    set.seed(42); x <- rnorm(b$n)
    bm <- bench::mark(
      robscale:::cpp_scale_ensemble(x, N_BOOT),
      min_iterations = 5L, check = FALSE
    )
    med_post_ns <- as.numeric(bm$median[1]) * 1e9
    ratio <- med_post_ns / b$median_ns
    is_stable <- b$n > 512L
    pass <- ratio <= RATIO_LIMIT
    if (is_stable && !pass) stable_fails <- stable_fails + 1L
    cat(sprintf("n=%5d: %.4f %s %s\n", b$n, ratio,
                if (pass) "PASS" else "FAIL",
                if (is_stable) "(stable)" else "(noisy)"))
  }
  if (stable_fails > 0L) {
    cat("GATE FAILED\n"); quit(status = 1)
  } else {
    cat("Gate passed.\n")
  }
}
