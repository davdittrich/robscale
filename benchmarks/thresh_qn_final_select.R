# Threshold calibration: Qn final diff-window selection crossover (FR vs pdqselect)
#
# Purpose: measure the n at which pdqselect overtakes Floyd-Rivest for the
# final floyd_rivest_select / pdqselect call at qn_estimator.cpp line 438.
#
# Working-set context: Qn final selection has 1 active diffs array plus
# sorted_x/work/bounds warm in cache (divisor ~4, threshold ~16K on Ryzen).
#
# Usage:
#   Rscript benchmarks/thresh_qn_final_select.R
#   # On remote: scp to 192.168.1.43 then ssh dd@192.168.1.43 'Rscript ~/benchmarks/thresh_qn_final_select.R'

suppressPackageStartupMessages({
  library(bench)
  library(robscale)
})

pkg_dir <- if (file.exists("DESCRIPTION")) "." else ".."
stopifnot(file.exists(file.path(pkg_dir, "DESCRIPTION")))

message("Installing current package...")
system2("R", c("CMD", "INSTALL", "--no-multiarch", "--no-test-load", pkg_dir),
        stdout = FALSE, stderr = FALSE)

# Qn final selection only activates on the refinement path (n > qn_exact_threshold ~64).
# Use a sweep that covers the refinement path.
sizes <- c(500L, 1000L, 2000L, 5000L, 10000L, 20000L, 50000L, 100000L, 500000L)
n_iter <- 20L

gen_data <- function(n) { set.seed(42); rnorm(n) }

results <- vector("list", length(sizes))
for (i in seq_along(sizes)) {
  n <- sizes[i]
  x <- gen_data(n)
  bm <- bench::mark(
    qn(x),
    min_iterations = n_iter,
    memory = FALSE
  )
  results[[i]] <- data.frame(
    n = n,
    median_us = as.numeric(bm$median) * 1e6
  )
  cat(sprintf("  qn  n=%7d  %.1f us\n", n, results[[i]]$median_us))
}

df <- do.call(rbind, results)

out_file <- sprintf("benchmarks/thresh_qn_final_%s_%s.csv",
                    format(Sys.Date(), "%Y%m%d"),
                    ifelse(Sys.info()["sysname"] == "Darwin", "macos", "linux"))
write.csv(df, out_file, row.names = FALSE)
cat(sprintf("\nSaved: %s\n", out_file))

cfg <- get_qnsn_config()
cat(sprintf("\nCurrent pdq_qn_final_threshold: %d\n", cfg$pdq_qn_final_threshold))
cat(sprintf("L2 per core: %d bytes\n", cfg$l2_per_core))
cat("\nNote: the final selection is a small fraction of total Qn time at large n.\n")
cat("To isolate it, build two versions forcing FR-only vs pdqselect-only at line 438.\n")
