## bench/ensemble_gate_ws6.R
## Ensemble timing gate for WU-RS6 (OPT-2: vshaped_mad for MAD init).
##
## Usage:
##   Rscript bench/ensemble_gate_ws6.R          # capture pre-change baseline
##   Rscript bench/ensemble_gate_ws6.R CHECK    # compare post-change vs baseline
##
## Gate criterion: post/pre ratio <= 1.05 at ALL sizes (blocking).
## A ratio < 1.0 indicates improvement (expected for OPT-2).
## Noisy-zone: all ensemble sizes are "noisy" (< 1s total); failures
## are informational unless consistently reproduced.

library(robscale)
library(bench)

RATIO_LIMIT <- 1.05
BASELINE_FILE <- "bench/ensemble_ws6_pre.rds"
MODE <- if (length(commandArgs(trailingOnly = TRUE)) > 0 &&
            toupper(commandArgs(trailingOnly = TRUE)[1]) == "CHECK") "CHECK" else "CAPTURE"

ens_sizes <- c(100, 500, 1000)
N_BOOT    <- 500L

cat("=== ensemble gate WU-RS6 ===\n")
cat("Mode:", MODE, "\n")

if (MODE == "CAPTURE") {
  cat("Capturing pre-change baseline (current code)...\n")
  pre <- lapply(ens_sizes, function(n) {
    set.seed(42)
    x <- rnorm(n)
    bm <- bench::mark(
      robscale:::cpp_scale_ensemble(x, N_BOOT),
      min_iterations = 5L, check = FALSE
    )
    med_ns <- as.numeric(bm$median[1]) * 1e9
    cat(sprintf("  n=%4d: %.0f ns\n", n, med_ns))
    list(n = n, median_ns = med_ns)
  })
  saveRDS(pre, BASELINE_FILE)
  cat("Baseline saved to", BASELINE_FILE, "\n")
  cat("Now implement OPT-2 and run: Rscript bench/ensemble_gate_ws6.R CHECK\n")
} else {
  ## CHECK mode: compare post-change vs saved baseline
  if (!file.exists(BASELINE_FILE)) {
    cat("ERROR: baseline not found. Run without CHECK first.\n")
    quit(status = 1)
  }
  pre <- readRDS(BASELINE_FILE)
  cat(sprintf("Gate: ratio <= %.2f\n\n", RATIO_LIMIT))
  cat(sprintf("%-6s  %-12s  %-12s  %-8s  %s\n",
              "n", "post (ns)", "pre (ns)", "ratio", "status"))
  cat(strrep("-", 58), "\n")

  all_pass <- TRUE
  for (b in pre) {
    set.seed(42)
    x <- rnorm(b$n)
    bm <- bench::mark(
      robscale:::cpp_scale_ensemble(x, N_BOOT),
      min_iterations = 5L, check = FALSE
    )
    med_post_ns <- as.numeric(bm$median[1]) * 1e9
    ratio <- med_post_ns / b$median_ns
    pass  <- ratio <= RATIO_LIMIT
    if (!pass) all_pass <- FALSE
    cat(sprintf("%-6d  %12.0f  %12.0f  %8.4f  %s\n",
                b$n, med_post_ns, b$median_ns, ratio,
                if (pass) "PASS" else "FAIL"))
  }

  cat(strrep("-", 58), "\n")
  cat("Overall:", if (all_pass) "PASS" else "FAIL", "\n")

  if (!all_pass) {
    cat("\nGATE FAILED: post/pre ratio > 1.05 — OPT-2 introduced a regression.\n")
    quit(status = 1)
  } else {
    cat("\nGate passed. All post/pre ratios <=", RATIO_LIMIT, "\n")
  }
}
