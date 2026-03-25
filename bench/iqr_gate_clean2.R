## bench/iqr_gate_clean2.R
## Performance gate for WU-CLEAN-2 (R3 buf2 removal + R4 helper extraction).
## Compares iqr_scaled() on stable sizes n=1000,5000 against a saved baseline.
## Noisy small-n sizes are NOT gating — only stable zone (n>512) blocks.
##
## Usage:
##   Rscript bench/iqr_gate_clean2.R           # CAPTURE baseline
##   Rscript bench/iqr_gate_clean2.R CHECK     # CHECK against baseline
##
## Gate threshold: ratio <= 1.05 (5%) for ALL stable sizes.
## Exit status: 0 = PASS, 1 = FAIL.
library(bench)
devtools::load_all(quiet = TRUE)

RATIO_LIMIT   <- 1.05
BASELINE_FILE <- "bench/iqr_clean2_pre.rds"
MODE <- if (length(commandArgs(trailingOnly = TRUE)) > 0 &&
            toupper(commandArgs(trailingOnly = TRUE)[1]) == "CHECK") "CHECK" else "CAPTURE"

stable_sizes <- c(1000L, 5000L)
N_ITER       <- 5L   # min_iterations per bench::mark() call

cat("=== iqr_gate_clean2 ===\n")
cat("Mode:", MODE, "\n\n")

if (MODE == "CAPTURE") {
  pre <- lapply(stable_sizes, function(n) {
    set.seed(42L)
    x <- rnorm(n)
    bm <- bench::mark(iqr_scaled(x), min_iterations = N_ITER, check = FALSE)
    med_ns <- as.numeric(bm$median[1L]) * 1e9
    cat(sprintf("  n=%5d: %.0f ns\n", n, med_ns))
    list(n = n, median_ns = med_ns)
  })
  saveRDS(pre, BASELINE_FILE)
  cat("\nBaseline saved to", BASELINE_FILE, "\n")
} else {
  if (!file.exists(BASELINE_FILE)) {
    cat("ERROR: baseline file not found:", BASELINE_FILE, "\n")
    cat("Run without CHECK argument first to capture a baseline.\n")
    quit(status = 1L)
  }
  pre <- readRDS(BASELINE_FILE)
  stable_fails <- 0L

  cat(sprintf("%-6s  %12s  %12s  %6s  %s\n",
              "n", "baseline_ns", "current_ns", "ratio", "verdict"))

  for (b in pre) {
    set.seed(42L)
    x <- rnorm(b$n)
    bm <- bench::mark(iqr_scaled(x), min_iterations = N_ITER, check = FALSE)
    curr_ns <- as.numeric(bm$median[1L]) * 1e9
    ratio   <- curr_ns / b$median_ns
    pass    <- ratio <= RATIO_LIMIT
    verdict <- if (pass) "PASS" else sprintf("FAIL (limit=%.2f)", RATIO_LIMIT)
    if (!pass) stable_fails <- stable_fails + 1L
    cat(sprintf("n=%5d  %12.0f  %12.0f  %6.4f  %s\n",
                b$n, b$median_ns, curr_ns, ratio, verdict))
  }

  cat("\n")
  if (stable_fails > 0L) {
    cat(sprintf("GATE FAILED: %d stable-zone size(s) exceeded %.0f%% threshold.\n",
                stable_fails, RATIO_LIMIT * 100))
    quit(status = 1L)
  } else {
    cat("Gate passed.\n")
  }
}
