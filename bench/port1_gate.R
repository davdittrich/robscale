## bench/port1_gate.R
## Regression gate for WU-PORT-1 (portability fixes in qnsn_hardware_info.h)
##
## Usage:
##   Rscript bench/port1_gate.R          # CAPTURE baseline (run BEFORE changes)
##   Rscript bench/port1_gate.R CHECK    # CHECK post-change vs baseline

library(robscale)
library(bench)

RATIO_LIMIT <- 1.05
BASELINE_FILE <- "bench/port1_pre.rds"
MODE <- if (length(commandArgs(trailingOnly = TRUE)) > 0 &&
            toupper(commandArgs(trailingOnly = TRUE)[1]) == "CHECK") "CHECK" else "CAPTURE"

stable_sizes <- c(1000L, 5000L)
N_BOOT <- 200L

cat("=== port1 gate ===\n")
cat("Mode:", MODE, "\n")

if (MODE == "CAPTURE") {
  pre <- lapply(stable_sizes, function(n) {
    set.seed(42); x <- rnorm(n)
    bm <- bench::mark(robscale:::cpp_scale_ensemble(x, N_BOOT), min_iterations = 5L, check = FALSE)
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
    bm <- bench::mark(robscale:::cpp_scale_ensemble(x, N_BOOT), min_iterations = 5L, check = FALSE)
    med_post_ns <- as.numeric(bm$median[1]) * 1e9
    ratio <- med_post_ns / b$median_ns
    pass <- ratio <= RATIO_LIMIT
    if (!pass) stable_fails <- stable_fails + 1L
    cat(sprintf("n=%5d: %.4f %s\n", b$n, ratio, if (pass) "PASS" else "FAIL"))
  }
  if (stable_fails > 0L) {
    cat("GATE FAILED\n"); quit(status = 1)
  } else {
    cat("Gate passed.\n")
  }
}
