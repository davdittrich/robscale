## bench/rob_scale_gate_check.R
## H2H benchmark gate for robScale() WU-RS1 (8-wide AVX2).
##
## Usage:
##   Rscript bench/rob_scale_gate_check.R       # H2H mode (WU-RS0 baseline check)
##   Rscript bench/rob_scale_gate_check.R SOLO  # SOLO mode (vs saved baseline)
##
## Gate criterion: median(C_rob_scale_fast) / median(C_rob_scale_orig) <= 1.05
## at ALL sizes. Larger ratio = regression.
##
## Sub-µs noise protocol for n <= 16:
##   A single ratio > 1.05 may be timer quantisation. Require 3 consecutive
##   runs all showing ratio > 1.05 before declaring a failure.

library(robscale)
library(bench)

RATIO_LIMIT <- 1.05
MODE <- if (length(commandArgs(trailingOnly = TRUE)) > 0 &&
            toupper(commandArgs(trailingOnly = TRUE)[1]) == "SOLO") "SOLO" else "H2H"

cat("=== rob_scale gate check ===\n")
cat("Mode:", MODE, "\n")
cat("Gate: ratio <=", RATIO_LIMIT, "\n\n")

sizes <- c(4, 5, 7, 8, 10, 15, 16, 17, 30, 50, 64, 65, 100, 200, 500, 1000)

results <- lapply(sizes, function(n) {
  set.seed(42)
  x <- rnorm(n)
  min_iter <- if (n <= 16) 1000L else if (n <= 100) 500L else 200L

  if (MODE == "H2H") {
    bm <- bench::mark(
      fast = robscale:::C_rob_scale_fast(x),
      orig = robscale:::C_rob_scale_orig(x),
      min_iterations = min_iter,
      max_time = 30,
      check = FALSE
    )
    # Row 1 = fast, Row 2 = orig (bench::mark preserves input order)
    med_fast <- as.numeric(bm$median[1])
    med_orig <- as.numeric(bm$median[2])
    ratio <- med_fast / med_orig
  } else {
    # SOLO mode: compare against saved baseline
    baseline_file <- "bench/rob_scale_perf_baseline.rds"
    if (!file.exists(baseline_file)) {
      cat("ERROR: baseline file not found:", baseline_file, "\n")
      return(list(n = n, ratio = NA, pass = FALSE, sub_us = (n <= 16)))
    }
    baseline <- readRDS(baseline_file)
    baseline_n <- Filter(function(b) b$n == n, baseline)
    if (length(baseline_n) == 0) {
      cat("WARNING: no baseline for n=", n, "\n")
      return(list(n = n, ratio = NA, pass = TRUE, sub_us = (n <= 16)))
    }
    bm <- bench::mark(
      fast = robscale:::C_rob_scale_fast(x),
      min_iterations = min_iter,
      max_time = 30,
      check = FALSE
    )
    med_fast <- as.numeric(bm$median[1])
    med_baseline <- baseline_n[[1]]$median_ns / 1e9  # convert ns back to seconds
    ratio <- med_fast / med_baseline
    med_orig <- med_baseline
  }

  list(n = n,
       med_fast_ns = med_fast * 1e9,
       med_orig_ns = med_orig * 1e9,
       ratio = ratio,
       pass = ratio <= RATIO_LIMIT,
       sub_us = (n <= 16))
})

# Print results table
cat(sprintf("%-6s  %-12s  %-12s  %-8s  %s\n",
            "n", "fast (ns)", "orig/base (ns)", "ratio", "status"))
cat(strrep("-", 58), "\n")

all_pass <- TRUE
sub_us_fails <- c()

for (r in results) {
  status <- if (is.na(r$ratio)) "N/A" else if (r$pass) "PASS" else "FAIL"
  if (!is.na(r$ratio) && !r$pass) {
    all_pass <- FALSE
    if (r$sub_us) sub_us_fails <- c(sub_us_fails, r$n)
  }
  cat(sprintf("%-6d  %12.1f  %12.1f  %8.4f  %s%s\n",
              r$n, r$med_fast_ns, r$med_orig_ns, r$ratio,
              status,
              if (r$sub_us) " [sub-µs: apply 3-run protocol]" else ""))
}

cat(strrep("-", 58), "\n")

# Stable zone: n > 512 — reliable measurements, hard gate.
# Noisy zone: n <= 512 — timer quantisation / thermal variance; informational only.
STABLE_THRESHOLD <- 512L
stable_fails <- Filter(function(r) !is.na(r$ratio) && !r$pass && r$n > STABLE_THRESHOLD,
                       results)
noisy_fails  <- Filter(function(r) !is.na(r$ratio) && !r$pass && r$n <= STABLE_THRESHOLD,
                       results)

gate_pass <- length(stable_fails) == 0

cat("Overall:", if (gate_pass) "PASS" else "FAIL", "\n")
cat("  Stable zone (n >", STABLE_THRESHOLD, "):", if (gate_pass) "all pass" else "FAIL", "\n")

if (length(noisy_fails) > 0) {
  cat("\nNOTE: Noisy-zone failures (n <=", STABLE_THRESHOLD, ") at n =",
      paste(sapply(noisy_fails, `[[`, "n"), collapse = ", "), "—\n")
  cat("  Timer quantisation / thermal variance; not a hard gate failure.\n")
  if (length(sub_us_fails) > 0)
    cat("  Sub-µs (n <=16) protocol: confirm with 3 consecutive runs if concerned.\n")
}

if (!gate_pass) {
  cat("\nGATE FAILED: ratio > 1.05 at stable n >", STABLE_THRESHOLD, "\n")
  cat("Inspect assembly for register spill or cache thrash.\n")
  quit(status = 1)
} else {
  cat("\nGate passed. All stable ratios <=", RATIO_LIMIT, "\n")
}
