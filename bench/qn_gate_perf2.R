## bench/qn_gate_perf2.R
## Timing gate for WU-PERF-2 (AVX2 in TBB Qn diff-fill path)
##
## Usage:
##   Rscript bench/qn_gate_perf2.R          # CAPTURE baseline (BEFORE change)
##   Rscript bench/qn_gate_perf2.R CHECK    # CHECK post-change vs baseline

library(robscale)
library(bench)

## Skip entirely on non-AVX2 machines
if (!isTRUE(robscale:::get_qnsn_config()$has_avx2)) {
  cat("SKIP: no AVX2 support on this machine.\n")
  quit(status = 0)
}

RATIO_LIMIT <- 1.05
BASELINE_FILE <- "bench/qn_perf2_pre.rds"
MODE <- if (length(commandArgs(trailingOnly = TRUE)) > 0 &&
            toupper(commandArgs(trailingOnly = TRUE)[1]) == "CHECK") "CHECK" else "CAPTURE"

stable_sizes <- c(1000L, 5000L, 10000L)

cat("=== qn_gate_perf2 ===\n")
cat("Mode:", MODE, "\n")

if (MODE == "CAPTURE") {
  pre <- lapply(stable_sizes, function(n) {
    set.seed(42); x <- rnorm(n)
    bm <- bench::mark(qn(x), min_iterations = 5L, check = FALSE)
    med_ns <- as.numeric(bm$median[1]) * 1e9
    cat(sprintf("  n=%6d: %.0f ns\n", n, med_ns))
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
    bm <- bench::mark(qn(x), min_iterations = 5L, check = FALSE)
    med_post_ns <- as.numeric(bm$median[1]) * 1e9
    ratio <- med_post_ns / b$median_ns
    pass <- ratio <= RATIO_LIMIT
    if (!pass) stable_fails <- stable_fails + 1L
    cat(sprintf("n=%6d: %.4f %s\n", b$n, ratio, if (pass) "PASS" else "FAIL"))
  }
  if (stable_fails > 0L) {
    cat("GATE FAILED\n"); quit(status = 1)
  } else {
    cat("Gate passed.\n")
  }
}
