## bench/adm_gate_check.R
## H2H + SOLO benchmark gate for ADM optimization WUs (WU-ADM1–4).
##
## Usage:
##   Rscript bench/adm_gate_check.R         # H2H mode
##   Rscript bench/adm_gate_check.R SOLO    # SOLO mode (vs saved baseline)
##
## Gate criterion: median(C_adm_fast) / median(C_adm_orig) <= 1.05 at ALL sizes.
## Larger ratio = regression.
##
## Sub-µs noise protocol for n <= 64:
##   A single ratio > 1.05 may be timer quantisation. Require 3 consecutive
##   runs all showing ratio > 1.05 before declaring a real failure.

library(robscale)
library(bench)

RATIO_LIMIT <- 1.05
MODE <- if (length(commandArgs(trailingOnly = TRUE)) > 0 &&
            toupper(commandArgs(trailingOnly = TRUE)[1]) == "SOLO") "SOLO" else "H2H"

cat("=== adm gate check ===\n")
cat("Mode:", MODE, "\n")
cat("Gate: ratio <=", RATIO_LIMIT, "\n\n")

sizes <- c(8, 16, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 32768, 65536)

results <- lapply(sizes, function(n) {
  set.seed(42)
  x <- rnorm(n)
  min_iter <- if (n <= 64) 5000L else if (n <= 1024) 1000L else 200L

  if (MODE == "H2H") {
    bm <- bench::mark(
      fast = robscale:::C_adm_fast(x),
      orig = robscale:::C_adm_orig(x),
      min_iterations = min_iter,
      check = FALSE
    )
    # Row 1 = fast, Row 2 = orig (bench::mark preserves input order)
    med_fast <- as.numeric(bm$median[1])
    med_orig <- as.numeric(bm$median[2])
    ratio <- med_fast / med_orig
  } else {
    # SOLO mode: compare against saved baseline (C_adm_orig at WU-ADM0)
    baseline_file <- "benchmarks/adm_perf_baseline.rds"
    if (!file.exists(baseline_file)) {
      cat("ERROR: baseline file not found:", baseline_file, "\n")
      return(list(n = n, ratio = NA, pass = FALSE))
    }
    baseline <- readRDS(baseline_file)
    baseline_n <- baseline[baseline$size == n, ]
    if (nrow(baseline_n) == 0) {
      cat("WARNING: no baseline for n=", n, "\n")
      return(list(n = n, ratio = NA, pass = TRUE))
    }
    bm <- bench::mark(
      fast = robscale:::C_adm_fast(x),
      min_iterations = min_iter,
      check = FALSE
    )
    med_fast   <- as.numeric(bm$median[1])
    med_orig   <- mean(baseline_n$median_ns) / 1e9   # avg over seeds, convert ns→s
    ratio      <- med_fast / med_orig
  }

  list(n          = n,
       med_fast_ns = med_fast * 1e9,
       med_orig_ns = med_orig * 1e9,
       ratio       = ratio,
       pass        = ratio <= RATIO_LIMIT,
       sub_us      = (n <= 64))
})

# Print results table
cat(sprintf("%-8s  %-12s  %-14s  %-8s  %s\n",
            "n", "fast (ns)", "orig/base (ns)", "ratio", "status"))
cat(strrep("-", 62), "\n")

all_pass    <- TRUE
sub_us_fails <- c()

for (r in results) {
  status <- if (is.na(r$ratio)) "N/A" else if (r$pass) "PASS" else "FAIL"
  if (!is.na(r$ratio) && !r$pass) {
    all_pass <- FALSE
    if (r$sub_us) sub_us_fails <- c(sub_us_fails, r$n)
  }
  cat(sprintf("%-8d  %12.1f  %14.1f  %8.4f  %s%s\n",
              r$n, r$med_fast_ns, r$med_orig_ns, r$ratio,
              status,
              if (r$sub_us) " [sub-µs: apply 3-run protocol]" else ""))
}

cat(strrep("-", 62), "\n")
cat("Overall:", if (all_pass) "PASS" else "FAIL", "\n")

if (length(sub_us_fails) > 0) {
  cat("\nNOTE: Failures at n =", paste(sub_us_fails, collapse = ", "),
      "are at sub-µs scale.\n")
  cat("Sub-µs protocol: run this script 3 times and check if ALL 3 show\n")
  cat("ratio > 1.05 at these sizes before declaring a real failure.\n")
}

if (!all_pass) {
  cat("\nGATE FAILED: ratio > 1.05 detected.\n")
  quit(status = 1)
} else {
  cat("\nGate passed. All ratios <=", RATIO_LIMIT, "\n")
}
