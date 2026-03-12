source("benchmarks/run_benchmarks.R")
suppressPackageStartupMessages(library(dplyr))

cat("Building fast_mode...\n")
fast_lib <- tempfile("lib_fast_")
res_robscale <- benchmark_robscale(list(ROBSCALE_FAST=1), fast_lib)

cat("Running legacy benchmarks...\n")
res_legacy <- benchmark_legacy()

print_speedup <- function(fast_df, base_df, name, n_val) {
  fast <- fast_df %>% filter(as.character(expr) == name, n == n_val)
  base <- base_df %>% filter(as.character(expr) == name, n == n_val)
  if (nrow(fast) == 0 || nrow(base) == 0) return()
  
  speedup <- as.numeric(base$median) / as.numeric(fast$median)
  cat(sprintf("\n%s (n=%d) Speedup: %.3f\n", name, n_val, speedup))
  if (speedup < 1.0) cat("❌ REGRESSION\n") else cat("✅ OK\n")
  cat("  Base:", as.character(base$median), "\n")
  cat("  Fast:", as.character(fast$median), "\n")
}

print_speedup(res_robscale$m_estimators, res_legacy$revss, "adm", 4)
print_speedup(res_robscale$m_estimators, res_legacy$revss, "robLoc", 5)
print_speedup(res_robscale$m_estimators, res_legacy$revss, "robLoc", 6)
print_speedup(res_robscale$m_estimators, res_legacy$revss, "robScale", 5)
