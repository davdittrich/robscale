library(targets)
library(tarchetypes)

# Set target options:
tar_option_set(
  packages = c("bench", "sessioninfo", "ggplot2", "dplyr", "tidyr", "withr", "robustbase", "revss", "boot", "purrr"),
  format = "rds"
)

# Source the functions
source("R/run_benchmarks.R")
source("R/make_plots.R")
source("R/analyze_results.R")

# Pipeline
list(
  # Benchmark legacy (robustbase, revss)
  tar_target(legacy_bench, benchmark_legacy()),
  
  # Benchmark robscale unoptimized
  tar_target(rob_unoptimized, benchmark_robscale(install_env = list(ROBSCALE_FAST = "0"))),
  
  # Benchmark robscale optimized
  tar_target(rob_optimized, benchmark_robscale(install_env = list(ROBSCALE_FAST = "1"))),
  
  # Rigorous statistical analysis
  tar_target(analyzed_results, analyze_benchmarks(rob_optimized, rob_unoptimized, legacy_bench)),
  
  tar_target(benchmark_report, {
    list(
      fast = rob_optimized,
      slow = rob_unoptimized,
      legacy = legacy_bench,
      analyzed = analyzed_results
    )
  }),
  
  # Create the speedup plot
  tar_target(speedup_figure, plot_benchmarks(analyzed_results)),
  
  # Create the fast vs slow comparison plot
  tar_target(fast_slow_figure, plot_fast_slow(analyzed_results)),
  
  # Render README.qmd
  tar_quarto(readme, "README.qmd")
)
