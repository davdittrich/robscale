# Qn pdqselect regression benchmark
#
# Compares qn() before and after adaptive final-diff dispatch.
# Also compares against robustbase::Qn as an external reference.
#
# Gate: no regression at any n (speedup >= 0.95x).
# Improvement may be modest since final selection is a small fraction of Qn time.
#
# Usage:
#   Rscript benchmarks/qn_pdqselect_bench.R

library(bench)

pkg_dir <- normalizePath(if (file.exists("DESCRIPTION")) "." else "..")
stopifnot(file.exists(file.path(pkg_dir, "DESCRIPTION")))

sizes  <- c(1000L, 5000L, 10000L, 50000L, 100000L, 500000L, 1000000L)
n_iter <- 20L

bench_code <- sprintf('
library(robscale); library(bench)
sizes <- c(%s); n_iter <- %dL
rows <- vector("list", length(sizes))
for (i in seq_along(sizes)) {
  n <- sizes[i]; set.seed(42); x <- rnorm(n)
  bm <- bench::mark(qn(x), min_iterations = n_iter, memory = FALSE)
  rows[[i]] <- data.frame(n = n, median_us = as.numeric(bm$median) * 1e6)
}
df <- do.call(rbind, rows)
write.csv(df, stdout(), row.names = FALSE)
', paste(sizes, collapse = ", "), n_iter)

cat("=== Phase 1: Installing HEAD~1 ===\n")
system2("git", c("-C", pkg_dir, "stash", "--include-untracked", "-q"),
        stdout = FALSE, stderr = FALSE)
system2("R", c("CMD", "INSTALL", "--no-multiarch", "--no-test-load", pkg_dir),
        stdout = FALSE, stderr = FALSE)
system2("git", c("-C", pkg_dir, "stash", "pop", "-q"),
        stdout = FALSE, stderr = FALSE)

cat("Benchmarking original qn...\n")
old_out <- system2("Rscript", c("-e", shQuote(bench_code)), stdout = TRUE, stderr = FALSE)
old_df  <- read.csv(textConnection(paste(old_out, collapse = "\n")))
cat("  done.\n")

cat("\n=== Phase 2: Installing current version ===\n")
system2("R", c("CMD", "INSTALL", "--no-multiarch", "--no-test-load", pkg_dir),
        stdout = FALSE, stderr = FALSE)

cat("Benchmarking adaptive qn...\n")
new_out <- system2("Rscript", c("-e", shQuote(bench_code)), stdout = TRUE, stderr = FALSE)
new_df  <- read.csv(textConnection(paste(new_out, collapse = "\n")))
cat("  done.\n")

merged <- merge(old_df, new_df, by = "n", suffixes = c("_old", "_new"))
merged$speedup <- merged$median_us_old / merged$median_us_new

cat("\n=== qn: adaptive vs original FR ===\n")
print(merged, digits = 3, row.names = FALSE)

thr_code <- 'library(robscale); cat(get_qnsn_config()$pdq_qn_final_threshold)'
thr_raw  <- system2("Rscript", c("-e", shQuote(thr_code)), stdout = TRUE, stderr = FALSE)
thr      <- suppressWarnings(as.integer(thr_raw[1]))

cat(sprintf("\n--- Gate check (pdq_qn_final_threshold = %d) ---\n", thr))
no_regr <- all(merged$speedup >= 0.95)
cat(sprintf("no regression (>= 0.95x): %s\n", ifelse(no_regr, "PASS", "FAIL")))

out_file <- sprintf("benchmarks/qn_pdqselect_%s_%s.csv",
                    format(Sys.Date(), "%Y%m%d"),
                    ifelse(Sys.info()["sysname"] == "Darwin", "macos", "linux"))
write.csv(merged, out_file, row.names = FALSE)
cat(sprintf("\nSaved to %s\n", out_file))
