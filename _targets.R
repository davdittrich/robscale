library(targets)
library(tarchetypes)

# Set target options:
tar_option_set(
  packages = c("bench", "sessioninfo", "ggplot2", "dplyr", "tidyr", "withr", "robustbase", "revss", "boot", "purrr", "callr"),
  format = "rds"
)

# Source the functions
source("benchmarks/run_benchmarks.R")
source("benchmarks/make_plots.R")
source("benchmarks/analyze_results.R")

# Pipeline
list(
  # Track source files
  tar_target(src_files, list.files("src", recursive = TRUE, full.names = TRUE), format = "file"),
  tar_target(r_files, list.files("R", recursive = TRUE, full.names = TRUE), format = "file"),
  
  # Benchmark legacy (robustbase, revss)
  tar_target(legacy_bench, benchmark_legacy()),
  
  # Benchmark robscale unoptimized
  tar_target(rob_unoptimized, {
    force(src_files)
    force(r_files)
    benchmark_robscale(install_env = list(ROBSCALE_FAST = "0"))
  }),
  
  # Benchmark robscale optimized
  tar_target(rob_optimized, {
    force(src_files)
    force(r_files)
    benchmark_robscale(install_env = list(ROBSCALE_FAST = "1"))
  }),
  
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
  tar_target(speedup_figure_obj, plot_benchmarks(analyzed_results)),
  tar_target(speedup_figure, {
    path <- "benchmarks/speedup_fig.png"
    ggsave(path, speedup_figure_obj, width = 12, height = 6)
    path
  }, format = "file"),
  
  # Create the fast vs slow comparison plot
  tar_target(fast_slow_figure_obj, plot_fast_slow(analyzed_results)),
  tar_target(fast_slow_figure, {
    path <- "benchmarks/fast_slow_fig.png"
    ggsave(path, fast_slow_figure_obj, width = 12, height = 6)
    path
  }, format = "file"),

  
  # Render README.qmd
  tar_quarto(readme, "README.qmd")
)
