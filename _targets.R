library(targets)
library(tarchetypes)

# Set target options:
tar_option_set(
  packages = c("bench", "sessioninfo", "ggplot2", "dplyr", "tidyr", "withr", "robustbase", "revss", "boot", "purrr", "callr", "Hmisc", "GiniDistance", "collapse"),
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

  # Benchmark robscale (SIMD and thresholds are now always auto-detected)
  tar_target(rob_optimized, {
    force(src_files)
    force(r_files)
    benchmark_robscale()
  }),

  # Rigorous statistical analysis
  tar_target(analyzed_results, analyze_benchmarks(rob_optimized, legacy_bench)),

  tar_target(benchmark_report, {
    list(
      fast = rob_optimized,
      legacy = legacy_bench,
      analyzed = analyzed_results
    )
  }),

  # Create the speedup plot
  tar_target(speedup_figure_obj, plot_benchmarks(analyzed_results)),
  tar_target(speedup_figure, {
    path <- "benchmarks/speedup_fig.png"
    ggsave(path, speedup_figure_obj, width = 18, height = 6)
    path
  }, format = "file"),


  # Render README.qmd
  tar_quarto(readme_qmd, "README.qmd"),

  # Post-process README.md for GitHub Mermaid compatibility
  tar_target(readme, {
    # Ensure readme_qmd is built
    force(readme_qmd)
    path <- "README.md"
    # Quarto's gfm writer often puts a space in ``` mermaid
    # which GitHub doesn't recognize. We strip it here.
    system2("sed", c("-i", "'s/``` mermaid/```mermaid/g'", path))
    path
  }, format = "file")
)
