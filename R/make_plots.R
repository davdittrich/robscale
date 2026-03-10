library(ggplot2)
library(dplyr)
library(tidyr)

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

#' Create the speedup plot
#' @param rob_fast Optimized robscale results
#' @param rob_slow Unoptimized robscale results
#' @param legacy Legacy package results
plot_benchmarks <- function(rob_fast, rob_slow, legacy) {
  # Helper to process mark results
  process_mark <- function(bm, label) {
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
    mutate(speedup = revss / `robscale (FAST=1)`)
  
  p1 <- ggplot(df_a, aes(x = n, y = speedup, color = expr)) +
    geom_line() + geom_point() +
    scale_x_log10() +
    labs(title = "Panel A: M-estimators vs revss", y = "Speedup (x)") +
    theme_minimal()
    
  # Panel B: Scale estimators vs robustbase
  df_b <- bind_rows(df_fast_scale, df_legacy_robustbase) %>%
    pivot_wider(names_from = label, values_from = median_ms) %>%
    mutate(speedup = robustbase / `robscale (FAST=1)`)
    
  p2 <- ggplot(df_b, aes(x = n, y = speedup, color = expr)) +
    geom_line() + geom_point() +
    scale_x_log10() +
    labs(title = "Panel B: Scale vs robustbase", y = "Speedup (x)") +
    theme_minimal()
    
  p1 + p2 + plot_layout(guides = "collect")
}
