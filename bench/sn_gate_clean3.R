## bench/sn_gate_clean3.R
## Timing gate for WU-CLEAN-3 (sn_inner_serial extraction)
library(robscale)
library(bench)

RATIO_LIMIT   <- 1.05
BASELINE_FILE <- "bench/sn_perf_clean3_baseline.rds"
MODE <- if (length(commandArgs(trailingOnly = TRUE)) > 0 &&
            toupper(commandArgs(trailingOnly = TRUE)[1]) == "CHECK") "CHECK" else "CAPTURE"

stable_sizes <- c(1000L, 5000L)
cat("=== sn_gate_clean3 ===\nMode:", MODE, "\n")

if (MODE == "CAPTURE") {
  pre <- lapply(stable_sizes, function(n) {
    set.seed(42)
    x <- rnorm(n)
    bm <- bench::mark(sn(x), min_iterations = 5L, check = FALSE)
    med_ns <- as.numeric(bm$median[1]) * 1e9
    cat(sprintf("  n=%5d: %.0f ns\n", n, med_ns))
    list(n = n, median_ns = med_ns)
  })
  saveRDS(pre, BASELINE_FILE)
  cat("Baseline saved to", BASELINE_FILE, "\n")
} else {
  if (!file.exists(BASELINE_FILE)) {
    cat("ERROR: no baseline file found at", BASELINE_FILE, "\n")
    quit(status = 1)
  }
  pre <- readRDS(BASELINE_FILE)
  stable_fails <- 0L
  for (b in pre) {
    set.seed(42)
    x <- rnorm(b$n)
    bm <- bench::mark(sn(x), min_iterations = 5L, check = FALSE)
    post_ns <- as.numeric(bm$median[1]) * 1e9
    ratio <- post_ns / b$median_ns
    pass  <- ratio <= RATIO_LIMIT
    if (!pass) stable_fails <- stable_fails + 1L
    cat(sprintf("n=%5d: ratio=%.4f (pre=%.0f ns, post=%.0f ns) %s\n",
                b$n, ratio, b$median_ns, post_ns,
                if (pass) "PASS" else "FAIL"))
  }
  if (stable_fails > 0L) {
    cat(sprintf("GATE FAILED: %d stable size(s) exceeded ratio limit %.2f\n",
                stable_fails, RATIO_LIMIT))
    quit(status = 1)
  } else {
    cat("Gate passed.\n")
  }
}
