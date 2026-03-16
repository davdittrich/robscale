# Threshold calibration: Sn lowmedian selection crossover (FR vs pdqselect)
#
# Purpose: measure the n at which pdqselect overtakes Floyd-Rivest for the
# lowmedian selection on inner_medians inside sn_estimator.cpp.
#
# Working-set context: Sn lowmedian uses 1 active inner_medians array +
# 1 residual sorted_x array (divisor ~2, threshold ~32K on Ryzen).
#
# Usage:
#   Rscript benchmarks/thresh_sn_lowmedian.R
#   # On remote: scp to 192.168.1.43 then ssh dd@192.168.1.43 'Rscript ~/benchmarks/thresh_sn_lowmedian.R'

suppressPackageStartupMessages({
  library(bench)
  library(robscale)
})

pkg_dir <- if (file.exists("DESCRIPTION")) "." else ".."
stopifnot(file.exists(file.path(pkg_dir, "DESCRIPTION")))

message("Installing current package...")
system2("R", c("CMD", "INSTALL", "--no-multiarch", "--no-test-load", pkg_dir),
        stdout = FALSE, stderr = FALSE)

sizes <- c(500L, 1000L, 2000L, 5000L, 10000L, 20000L, 50000L, 100000L, 500000L)
n_iter <- 30L

# inner_medians values depend on sorted_x spacing; use rnorm (representative)
gen_data <- function(n) { set.seed(42); rnorm(n) }

results <- vector("list", length(sizes))
for (i in seq_along(sizes)) {
  n <- sizes[i]
  x <- gen_data(n)
  bm <- bench::mark(
    sn(x),
    min_iterations = n_iter,
    memory = FALSE
  )
  results[[i]] <- data.frame(
    n = n,
    median_us = as.numeric(bm$median) * 1e6
  )
  cat(sprintf("  sn  n=%7d  %.1f us\n", n, results[[i]]$median_us))
}

df <- do.call(rbind, results)

out_file <- sprintf("benchmarks/thresh_sn_lowmedian_%s_%s.csv",
                    format(Sys.Date(), "%Y%m%d"),
                    ifelse(Sys.info()["sysname"] == "Darwin", "macos", "linux"))
write.csv(df, out_file, row.names = FALSE)
cat(sprintf("\nSaved: %s\n", out_file))

cfg <- get_qnsn_config()
cat(sprintf("\nCurrent pdq_lowmedian_threshold: %d\n", cfg$pdq_lowmedian_threshold))
cat(sprintf("L2 per core: %d bytes\n", cfg$l2_per_core))
cat("\nNote: to isolate the lowmedian path, build two versions that force\n")
cat("FR-only vs pdqselect-only in adaptive_lowmedian_select.\n")
