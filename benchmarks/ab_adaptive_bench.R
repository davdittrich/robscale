#!/usr/bin/env Rscript
# A/B benchmark: pre-adaptive (FR-only) vs adaptive (pdqselect dispatch)
# for robScale, sn, qn across n = 1K..500K.
#
# Requires both tarballs in the working directory:
#   robscale_old_0.2.0.tar.gz  — FR-only (pre-adaptive)
#   robscale_0.2.0.tar.gz      — adaptive dispatch
#
# Usage (from run directory with both tarballs):
#   Rscript --vanilla ab_adaptive_bench.R
#
# Gate: speedup >= 0.95x at all n (no regression).
#
# Finding (Ryzen 9 5900HX, 512KB L2):
#   All three estimators show ~1.00x speedup above threshold.
#   robScale is tanh-dominated; sn is inner_medians-dominated; qn is
#   refinement-loop-dominated. The selection step is a small tail in all
#   three cases, so the pdqselect dispatch is correct but not the bottleneck.

library(bench)

# Install both versions into temp libraries
old_lib <- tempfile("robscale_old_")
new_lib <- tempfile("robscale_new_")
dir.create(old_lib); dir.create(new_lib)

cat("Installing OLD (FR-only) version...\n")
r_bin <- file.path(R.home("bin"), "Rscript")
install_pkg <- function(tar, lib) {
  res <- system2(R.home("bin/R"),
    args = c("CMD", "INSTALL", tar, paste0("--library=", lib),
             "--no-multiarch", "--no-test-load"),
    stdout = TRUE, stderr = TRUE)
  if (!is.null(attr(res, "status")) && attr(res, "status") != 0)
    stop("Install failed: ", paste(tail(res, 5), collapse="\n"))
  invisible(res)
}
install_pkg("robscale_old_0.2.0.tar.gz", old_lib)
cat("Installing NEW (adaptive) version...\n")
install_pkg("robscale_0.2.0.tar.gz", new_lib)

# ---------------------------------------------------------------------------
# Benchmark runner — runs in subprocess to avoid namespace collisions
# ---------------------------------------------------------------------------
bench_version <- function(lib, sizes, n_iter) {
  code <- sprintf('
suppressPackageStartupMessages({
  .libPaths(c(%s, .libPaths()))
  library(robscale)
})
library(bench)
sizes  <- c(%s)
n_iter <- %dL
rows   <- vector("list", length(sizes) * 3L)
k <- 1L
for (n in sizes) {
  set.seed(42); x <- rnorm(n)
  for (est in c("robScale", "sn", "qn")) {
    fn <- get(est, envir = asNamespace("robscale"))
    bm <- bench::mark(fn(x), min_iterations = n_iter, memory = FALSE, check = FALSE)
    rows[[k]] <- data.frame(estimator = est, n = n,
      median_us = as.numeric(bm$median) * 1e6)
    k <- k + 1L
  }
}
df <- do.call(rbind, rows)
write.csv(df, stdout(), row.names = FALSE)
',
    shQuote(lib),
    paste(sizes, collapse = ", "),
    n_iter
  )
  out <- system2(file.path(R.home("bin"), "Rscript"),
    args = c("--vanilla", "-e", shQuote(code)),
    stdout = TRUE, stderr = FALSE)
  read.csv(textConnection(paste(out, collapse = "\n")))
}

sizes  <- c(1000L, 5000L, 10000L, 20000L, 50000L, 100000L, 500000L)
n_iter <- 25L

cat("\nBenchmarking OLD (FR-only)...\n")
old_df <- bench_version(old_lib, sizes, n_iter)
cat("Benchmarking NEW (adaptive)...\n")
new_df <- bench_version(new_lib, sizes, n_iter)

# ---------------------------------------------------------------------------
# Merge and compute speedup
# ---------------------------------------------------------------------------
merged <- merge(old_df, new_df, by = c("estimator", "n"), suffixes = c("_old", "_new"))
merged$speedup <- merged$median_us_old / merged$median_us_new

# Get thresholds from new version
thr_code <- sprintf('
.libPaths(c(%s, .libPaths()))
library(robscale)
cfg <- robscale:::get_qnsn_config()
cat(cfg$pdq_robscale_threshold, cfg$pdq_lowmedian_threshold, cfg$pdq_qn_final_threshold)
', shQuote(new_lib))
thr_raw <- system2(file.path(R.home("bin"), "Rscript"),
  args = c("--vanilla", "-e", shQuote(thr_code)), stdout = TRUE, stderr = FALSE)
thrs <- as.integer(strsplit(trimws(thr_raw[1]), "\\s+")[[1]])
thr_rs <- thrs[1]; thr_sn <- thrs[2]; thr_qn <- thrs[3]

# ---------------------------------------------------------------------------
# Print results table
# ---------------------------------------------------------------------------
cat("\n=============================================================\n")
cat(" A/B: Adaptive pdqselect dispatch vs FR-only (Ryzen)\n")
cat("=============================================================\n")
for (est in c("robScale", "sn", "qn")) {
  sub <- merged[merged$estimator == est, ]
  sub <- sub[order(sub$n), ]
  cat(sprintf("\n%s:\n", est))
  cat(sprintf("  %8s  %10s  %10s  %7s\n", "n", "old_us", "new_us", "speedup"))
  for (i in seq_len(nrow(sub))) {
    r <- sub[i, ]
    flag <- ""
    thr <- switch(est, robScale = thr_rs, sn = thr_sn, qn = thr_qn)
    if (r$n >= thr) flag <- " [pdqselect active]"
    cat(sprintf("  %8d  %10.1f  %10.1f  %6.2fx%s\n",
                r$n, r$median_us_old, r$median_us_new, r$speedup, flag))
  }
}

# ---------------------------------------------------------------------------
# Gate check
# ---------------------------------------------------------------------------
cat("\n--- Gate check ---\n")
all_pass <- TRUE
for (est in c("robScale", "sn", "qn")) {
  sub <- merged[merged$estimator == est, ]
  no_regr <- all(sub$speedup >= 0.95)
  cat(sprintf("  %-10s  no_regression (>=0.95x): %s\n",
              est, ifelse(no_regr, "PASS", "FAIL")))
  if (!no_regr) all_pass <- FALSE
}
cat(sprintf("\nOVERALL: %s\n", ifelse(all_pass, "PASS", "FAIL")))

# Save
platform <- ifelse(Sys.info()["sysname"] == "Darwin", "macos", "linux")
out_file <- sprintf("ab_adaptive_%s_%s.csv", format(Sys.Date(), "%Y%m%d"), platform)
write.csv(merged, out_file, row.names = FALSE)
cat(sprintf("Saved to %s\n", out_file))

if (!all_pass) quit(status = 1)
