# pdqselect A/B regression benchmark
#
# Strategy: install original 0.2.0 (HEAD of branch), time IQR/MAD,
# then install the pdqselect working-tree version, time again.
# Compares within a single R session to control for machine load.

library(bench)

pkg_dir <- normalizePath("..")
if (!file.exists(file.path(pkg_dir, "DESCRIPTION"))) {
  pkg_dir <- normalizePath(".")
}
stopifnot(file.exists(file.path(pkg_dir, "DESCRIPTION")))

sizes <- c(1000L, 5000L, 10000L, 50000L, 100000L, 500000L, 1000000L)
n_iter <- 50L

# --- Phase 1: benchmark the ORIGINAL 0.2.0 (committed HEAD) -----------------
cat("=== Phase 1: Installing original 0.2.0 from HEAD ===\n")
system2("git", c("-C", pkg_dir, "stash", "--include-untracked", "-q"),
        stdout = FALSE, stderr = FALSE)
system2("R", c("CMD", "INSTALL", "--no-multiarch", "--no-test-load", pkg_dir),
        stdout = FALSE, stderr = FALSE)
system2("git", c("-C", pkg_dir, "stash", "pop", "-q"),
        stdout = FALSE, stderr = FALSE)

# Load the old version in a clean subprocess to collect timings
old_bench_code <- sprintf('
library(robscale)
library(bench)
sizes <- c(%s)
n_iter <- %dL
results <- vector("list", length(sizes) * 2L)
k <- 1L
for (n in sizes) {
  set.seed(42)
  x <- rnorm(n)
  bm <- bench::mark(iqr_scaled(x, constant = 1), min_iterations = n_iter)
  results[[k]] <- data.frame(estimator = "iqr", n = n,
    median_us = as.numeric(bm$median) * 1e6)
  k <- k + 1L
  bm <- bench::mark(mad_scaled(x, constant = 1), min_iterations = n_iter)
  results[[k]] <- data.frame(estimator = "mad", n = n,
    median_us = as.numeric(bm$median) * 1e6)
  k <- k + 1L
}
df <- do.call(rbind, results)
write.csv(df, stdout(), row.names = FALSE)
', paste(sizes, collapse = ", "), n_iter)

cat("Benchmarking original...\n")
old_out <- system2("Rscript", c("-e", shQuote(old_bench_code)),
                   stdout = TRUE, stderr = FALSE)
old_df <- read.csv(textConnection(paste(old_out, collapse = "\n")))
cat("  done.\n")

# --- Phase 2: install pdqselect working-tree version -------------------------
cat("\n=== Phase 2: Installing pdqselect version ===\n")
system2("R", c("CMD", "INSTALL", "--no-multiarch", "--no-test-load", pkg_dir),
        stdout = FALSE, stderr = FALSE)

new_bench_code <- old_bench_code  # identical measurement code

cat("Benchmarking pdqselect...\n")
new_out <- system2("Rscript", c("-e", shQuote(new_bench_code)),
                   stdout = TRUE, stderr = FALSE)
new_df <- read.csv(textConnection(paste(new_out, collapse = "\n")))
cat("  done.\n")

# --- Phase 3: compare --------------------------------------------------------
merged <- merge(old_df, new_df, by = c("estimator", "n"),
                suffixes = c("_old", "_new"))
merged$speedup <- merged$median_us_old / merged$median_us_new

cat("\n=== pdqselect vs original 0.2.0 FR ===\n")
print(merged, digits = 3, row.names = FALSE)

cat("\n--- Gate check ---\n")
large <- merged[merged$n >= 50000, ]
iqr_pass <- all(large$speedup[large$estimator == "iqr"] >= 1.5)
mad_pass <- all(large$speedup[large$estimator == "mad"] >= 1.5)
no_regr  <- all(merged$speedup >= 0.95)

cat(sprintf("IQR >= 1.5x at n>=50K:  %s\n", ifelse(iqr_pass, "PASS", "FAIL")))
cat(sprintf("MAD >= 1.5x at n>=50K:  %s\n", ifelse(mad_pass, "PASS", "FAIL")))
cat(sprintf("No regression (<0.95):  %s\n", ifelse(no_regr, "PASS", "FAIL")))

out_file <- sprintf("benchmarks/pdqselect_regression_%s_%s.csv",
                     format(Sys.Date(), "%Y%m%d"),
                     ifelse(Sys.info()["sysname"] == "Darwin", "macos", "linux"))
write.csv(merged, out_file, row.names = FALSE)
cat(sprintf("\nSaved to %s\n", out_file))
