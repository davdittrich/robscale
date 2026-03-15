# robScale pdqselect regression benchmark
#
# Installs HEAD~1 (pre-adaptive), benchmarks robScale, then installs current
# working-tree version, benchmarks again, and compares.
#
# Gates:
#   - speedup >= 1.5x at n >= pdq_robscale_threshold  (Ryzen only; Apple has huge L2)
#   - speedup >= 0.95x at all n  (no regression anywhere)
#
# Usage:
#   Rscript benchmarks/robscale_pdqselect_bench.R
#   # On remote: scp + ssh dd@192.168.1.43 'Rscript ~/benchmarks/robscale_pdqselect_bench.R'

library(bench)

pkg_dir <- normalizePath(if (file.exists("DESCRIPTION")) "." else "..")
stopifnot(file.exists(file.path(pkg_dir, "DESCRIPTION")))

sizes  <- c(1000L, 5000L, 10000L, 50000L, 100000L, 500000L, 1000000L)
n_iter <- 30L

bench_code <- sprintf('
library(robscale); library(bench)
sizes <- c(%s); n_iter <- %dL
rows <- vector("list", length(sizes))
for (i in seq_along(sizes)) {
  n <- sizes[i]; set.seed(42); x <- rnorm(n)
  bm <- bench::mark(robScale(x), min_iterations = n_iter, memory = FALSE)
  rows[[i]] <- data.frame(n = n, median_us = as.numeric(bm$median) * 1e6)
}
df <- do.call(rbind, rows)
write.csv(df, stdout(), row.names = FALSE)
', paste(sizes, collapse = ", "), n_iter)

# --- Phase 1: benchmark HEAD~1 (pre-adaptive version) ---
cat("=== Phase 1: Installing HEAD~1 ===\n")
system2("git", c("-C", pkg_dir, "stash", "--include-untracked", "-q"),
        stdout = FALSE, stderr = FALSE)
system2("R", c("CMD", "INSTALL", "--no-multiarch", "--no-test-load", pkg_dir),
        stdout = FALSE, stderr = FALSE)
system2("git", c("-C", pkg_dir, "stash", "pop", "-q"),
        stdout = FALSE, stderr = FALSE)

cat("Benchmarking original...\n")
old_out <- system2("Rscript", c("-e", shQuote(bench_code)), stdout = TRUE, stderr = FALSE)
old_df  <- read.csv(textConnection(paste(old_out, collapse = "\n")))
cat("  done.\n")

# --- Phase 2: benchmark current working tree ---
cat("\n=== Phase 2: Installing current version ===\n")
system2("R", c("CMD", "INSTALL", "--no-multiarch", "--no-test-load", pkg_dir),
        stdout = FALSE, stderr = FALSE)

cat("Benchmarking adaptive version...\n")
new_out <- system2("Rscript", c("-e", shQuote(bench_code)), stdout = TRUE, stderr = FALSE)
new_df  <- read.csv(textConnection(paste(new_out, collapse = "\n")))
cat("  done.\n")

# --- Compare ---
merged <- merge(old_df, new_df, by = "n", suffixes = c("_old", "_new"))
merged$speedup <- merged$median_us_old / merged$median_us_new

cat("\n=== robScale: adaptive vs original FR ===\n")
print(merged, digits = 3, row.names = FALSE)

# Gate: get threshold from installed package
thr_code <- 'library(robscale); cat(get_qnsn_config()$pdq_robscale_threshold)'
thr_raw  <- system2("Rscript", c("-e", shQuote(thr_code)), stdout = TRUE, stderr = FALSE)
thr      <- suppressWarnings(as.integer(thr_raw[1]))

cat(sprintf("\n--- Gate check (pdq_robscale_threshold = %d) ---\n", thr))
large   <- merged[!is.na(thr) & merged$n >= thr, ]
speedup_pass <- !is.na(thr) && nrow(large) == 0 || all(large$speedup >= 1.5)
no_regr      <- all(merged$speedup >= 0.95)

cat(sprintf("speedup >= 1.5x at n >= thr:  %s\n", ifelse(speedup_pass, "PASS", "FAIL")))
cat(sprintf("no regression (>= 0.95x):     %s\n", ifelse(no_regr, "PASS", "FAIL")))

out_file <- sprintf("benchmarks/robscale_pdqselect_%s_%s.csv",
                    format(Sys.Date(), "%Y%m%d"),
                    ifelse(Sys.info()["sysname"] == "Darwin", "macos", "linux"))
write.csv(merged, out_file, row.names = FALSE)
cat(sprintf("\nSaved to %s\n", out_file))
