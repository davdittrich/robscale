#' Comprehensive Legacy Benchmark: robscale vs All Reference Implementations
#' 
#' This script compares robscale estimators (sn, qn, robScale, robLoc, adm) 
#' against legacy counterparts in robustbase, fastqnsn, and revss.
#' 
#' Output: 
#' - benchmarks/legacy_benchmark_raw.csv
#' - benchmarks/legacy_benchmark_summary.csv
#' - benchmarks/legacy_speedup_ribbon.png
#' - benchmarks/legacy_raincloud.png

# ── 0. Dependencies ──────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(robscale)
  library(robustbase)
  library(fastqnsn)
  library(revss)
  library(microbenchmark)
  library(boot)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(sessioninfo)
})

# ── 1. Environment Record ───────────────────────────────────────────────────
cat("=== ENVIRONMENT RECORD ===\n")
si <- sessioninfo::session_info()
print(si)

governor_path <- "/sys/devices/system/cpu/cpu0/cpufreq/scaling_governor"
if (file.exists(governor_path)) {
  cat("CPU governor:", readLines(governor_path), "\n")
} else {
  cat("CPU governor: unknown (non-Linux?)\n")
}

cat("robscale config:\n")
print(robscale:::get_qnsn_config())
cat("\n")

# ── 2. Numerical Equivalence Checks ──────────────────────────────────────────
cat("=== NUMERICAL EQUIVALENCE CHECKS ===\n")
seeds <- 1:20
for (s in seeds) {
  set.seed(s)
  x <- rnorm(100)
  
  # fastqnsn (must be identical)
  stopifnot(all.equal(robscale::sn(x), fastqnsn::sn(x), tolerance = 1e-15))
  stopifnot(all.equal(robscale::qn(x), fastqnsn::qn(x), tolerance = 1e-15))
  
  # robustbase (algorithmic identity, allowing for constant/factor diffs)
  # robustbase uses rounded constants (1.1926, 2.2219) vs precise values in robscale
  # and different correction factors for small n. We check core results without corr.
  stopifnot(all.equal(robscale::sn(x, finite.corr = FALSE), 
                      robustbase::Sn(x, finite.corr = FALSE), tolerance = 1e-5))
  stopifnot(all.equal(robscale::qn(x, finite.corr = FALSE), 
                      robustbase::Qn(x, finite.corr = FALSE), tolerance = 1e-5))
  
  # revss wrappers (should be identical, allowing for floating point noise and iteration diffs)
  stopifnot(all.equal(robscale::robScale(x), revss::robScale(x), tolerance = 1e-7))
  stopifnot(all.equal(robscale::robLoc(x), revss::robLoc(x), tolerance = 1e-7))
  stopifnot(all.equal(robscale::adm(x), revss::adm(x), tolerance = 1e-7))
}
cat("Numerical equivalence confirmed (controlled checks).\n\n")

# ── 3. Benchmark Grid Definition ─────────────────────────────────────────────
comparisons <- list(
  list(name = "Sn", new = robscale::sn,
       legacy = list(
         list(name = "robustbase::Sn", fn = robustbase::Sn),
         list(name = "fastqnsn::sn",   fn = fastqnsn::sn)
       ),
       sizes = c(8, 10, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 
                 12288, 16384, 32768, 65536, 131072, 262144, 524288, 1048576, 10000000)),
  
  list(name = "Qn", new = robscale::qn,
       legacy = list(
         list(name = "robustbase::Qn", fn = robustbase::Qn),
         list(name = "fastqnsn::qn",   fn = fastqnsn::qn)
       ),
       sizes = c(8, 10, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 
                 12288, 16384, 32768, 65536, 131072, 262144, 524288, 1048576, 10000000)),
  
  list(name = "robScale", new = robscale::robScale,
       legacy = list(list(name = "revss::robScale", fn = revss::robScale)),
       sizes = c(8, 10, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 12288, 16384)),
  
  list(name = "robLoc", new = robscale::robLoc,
       legacy = list(list(name = "revss::robLoc", fn = revss::robLoc)),
       sizes = c(8, 10, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 12288, 16384)),
  
  list(name = "adm", new = robscale::adm,
       legacy = list(list(name = "revss::adm", fn = revss::adm)),
       sizes = c(8, 10, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 12288, 16384))
)

get_iters <- function(n) {
  if (n <= 128) 20000L
  else if (n <= 2048) 5000L
  else if (n <= 16384) 500L
  else if (n <= 1048576) 50L
  else 10L
}

# ── 4. Execute Benchmarks ─────────────────────────────────────────────────────
cat("=== EXECUTING BENCHMARK GRID ===\n")
all_raw <- list()

for (comp in comparisons) {
  cat(sprintf("[%s] Starting benchmark sweep...\n", comp$name))
  for (n in comp$sizes) {
    iters <- get_iters(n)
    cat(sprintf("  n = %d (%d iters)... ", n, iters))
    
    set.seed(42)
    x <- rnorm(n)
    
    # Bundle implementations to test
    exprs <- list(robscale = substitute(new_fn(x), list(new_fn = comp$new)))
    for (leg in comp$legacy) {
      exprs[[leg$name]] <- substitute(leg_fn(x), list(leg_fn = leg$fn))
    }
    
    # Run microbenchmark
    bm <- microbenchmark(list = exprs, times = iters)
    
    # Tidy results
    df <- as.data.frame(bm) %>%
      mutate(n = n, estimator = comp$name) %>%
      rename(implementation = expr, time_ns = time) %>%
      mutate(time_us = time_ns / 1000)
    
    all_raw[[paste(comp$name, n, sep="_")]] <- df
    cat("done\n")
  }
}

raw_df <- bind_rows(all_raw)
write.csv(raw_df, "benchmarks/legacy_benchmark_raw.csv", row.names = FALSE)
cat("\nRaw timings saved to benchmarks/legacy_benchmark_raw.csv\n\n")

# ── 5. Bootstrap Analysis ─────────────────────────────────────────────────────
cat("=== PERFORMING BOOTSTRAP ANALYSIS ===\n")

# Ratio of medians function
ratio_median <- function(d, i) {
  rob_median <- median(d[i, "robscale"])
  leg_median <- median(d[i, "legacy"])
  leg_median / rob_median
}

summary_list <- list()

# Group by estimator, legacy target, and n
combos <- raw_df %>% 
  filter(implementation != "robscale") %>%
  distinct(estimator, implementation, n)

for (i in 1:nrow(combos)) {
  row <- combos[i,]
  est <- row$estimator
  leg_name <- row$implementation
  n_val <- row$n
  
  cat(sprintf("  Analyzing %s vs %s at n=%d... ", est, leg_name, n_val))
  
  rob_times <- raw_df %>% filter(estimator == est, n == n_val, implementation == "robscale") %>% pull(time_us)
  leg_times <- raw_df %>% filter(estimator == est, n == n_val, implementation == leg_name) %>% pull(time_us)
  
  # Prepare for boot: a data frame with two columns
  # Since iterations might differ slightly (due to microbenchmark adaptive logic sometimes, 
  # though we fixed it), we take the first min(length) rows.
  len <- min(length(rob_times), length(leg_times))
  boot_data <- data.frame(robscale = rob_times[1:len], legacy = leg_times[1:len])
  
  # Bootstrap
  b <- boot(boot_data, ratio_median, R = 2000)
  ci <- tryCatch(boot.ci(b, type = "perc")$percent[4:5], error = function(e) c(NA, NA))
  
  summary_list[[i]] <- tibble(
    n = n_val,
    estimator = est,
    legacy_name = as.character(leg_name),
    robscale_med_us = median(rob_times),
    legacy_med_us = median(leg_times),
    speedup = median(leg_times) / median(rob_times),
    ci_lo = ci[1],
    ci_hi = ci[2],
    significant = !is.na(ci[1]) && (ci[1] > 1.05 || ci[2] < 0.95)
  )
  cat("done\n")
}

summary_df <- bind_rows(summary_list)
write.csv(summary_df, "benchmarks/legacy_benchmark_summary.csv", row.names = FALSE)
cat("\nSummary table saved to benchmarks/legacy_benchmark_summary.csv\n\n")

# ── 6. Visualizations ────────────────────────────────────────────────────────
cat("=== GENERATING PLOTS ===\n")

# Plot 1: Speedup Ribbon
p1 <- ggplot(summary_df, aes(x = n, y = speedup, color = legacy_name, fill = legacy_name)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
  geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.2, color = NA) +
  geom_line(size = 1) +
  geom_point() +
  facet_wrap(~estimator, scales = "free_y") +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_continuous(labels = label_number(suffix = "x")) +
  labs(title = "robscale Speedup over Legacy Implementations",
       subtitle = "Speedup = Median(Legacy) / Median(robscale). Ribbon is 95% bootstrap CI.",
       x = "Sample Size (n)", y = "Speedup Factor",
       color = "Legacy Source", fill = "Legacy Source") +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave("benchmarks/legacy_speedup_ribbon.png", p1, width = 10, height = 7, dpi = 300)

# Plot 2: Raincloud at selected n
# Selected sizes: 128 (exact/serial transition), 16384 (parallel region)
selected_sizes <- c(128, 16384)
rain_data <- raw_df %>% filter(n %in% selected_sizes)

if (requireNamespace("ggdist", quietly = TRUE)) {
  library(ggdist)
  p2 <- ggplot(rain_data, aes(y = implementation, x = time_us, fill = implementation)) +
    stat_halfeye(alpha = 0.7) +
    stat_dots(side = "bottom", alpha = 0.5) +
    facet_grid(estimator ~ n, scales = "free") +
    labs(title = "Timing Distributions at Selected n",
         subtitle = "x-axis is time (µs), scales vary by row and column.",
         x = "Time (µs)", y = NULL) +
    theme_minimal() +
    theme(legend.position = "none")
  
  ggsave("benchmarks/legacy_raincloud.png", p2, width = 12, height = 10, dpi = 300)
}

cat("Plots saved to benchmarks/legacy_speedup_ribbon.png and benchmarks/legacy_raincloud.png\n")
cat("\n=== BENCHMARK COMPLETE ===\n")
