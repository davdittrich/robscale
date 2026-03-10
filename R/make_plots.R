library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(scales)

#' Create the speedup plot
#' @param rob_fast Optimized robscale results
#' @param rob_slow Unoptimized robscale results
#' @param legacy Legacy package results
plot_benchmarks <- function(rob_fast, rob_slow, legacy) {
  # Helper to process mark results
  process_mark <- function(bm, label) {
    if (is.null(bm)) return(NULL)
    bm %>%
      mutate(
        expr = as.character(expression),
        median_ms = as.numeric(median) * 1000,
        label = label
      ) %>%
      select(n, expr, median_ms, label)
  }

  # Process all components
  df_fast_m <- process_mark(rob_fast$m_estimators, "robscale (FAST=1)")
  df_fast_scale <- process_mark(rob_fast$scale_estimators, "robscale (FAST=1)")
  
  df_legacy_revss <- process_mark(legacy$revss, "revss")
  df_legacy_robustbase <- process_mark(legacy$robustbase, "robustbase")
  
  # Panel A: M-estimators vs revss
  df_a <- bind_rows(df_fast_m, df_legacy_revss) %>%
    pivot_wider(names_from = label, values_from = median_ms) %>%
    filter(!is.na(revss), !is.na(`robscale (FAST=1)`)) %>%
    mutate(speedup = revss / `robscale (FAST=1)`)
  
  p1 <- ggplot(df_a, aes(x = n, y = speedup, color = expr)) +
    geom_line(size = 0.8) + geom_point(size = 1.5) +
    scale_x_log10(
      breaks = 10^(0:7),
      labels = trans_format("log10", math_format(10^.x))
    ) +
    scale_y_continuous(labels = label_number(suffix = "x")) +
    labs(
      title = "Panel A: Small-sample M-estimators",
      subtitle = "robscale vs revss",
      x = "Sample Size (n)",
      y = "Speedup Factor",
      color = "Estimator"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
    
  # Panel B: Scale estimators vs robustbase
  df_b <- bind_rows(df_fast_scale, df_legacy_robustbase) %>%
    pivot_wider(names_from = label, values_from = median_ms) %>%
    filter(!is.na(robustbase), !is.na(`robscale (FAST=1)`)) %>%
    mutate(speedup = robustbase / `robscale (FAST=1)`)
    
  p2 <- ggplot(df_b, aes(x = n, y = speedup, color = expr)) +
    geom_line(size = 0.8) + geom_point(size = 1.5) +
    scale_x_log10(
      breaks = 10^(0:7),
      labels = trans_format("log10", math_format(10^.x))
    ) +
    scale_y_continuous(labels = label_number(suffix = "x")) +
    labs(
      title = "Panel B: Scale Estimators",
      subtitle = "robscale vs robustbase",
      x = "Sample Size (n)",
      y = "Speedup Factor",
      color = "Estimator"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
    
  p1 + p2 + plot_layout(guides = "collect") & theme(legend.position = "bottom")
}
