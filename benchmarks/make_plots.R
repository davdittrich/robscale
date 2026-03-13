library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(scales)

#' Create the speedup plot
#' @param analyzed The analyzed results from analyze_benchmarks()
plot_benchmarks <- function(analyzed) {
  # Panel A: M-estimators vs revss (analyzed$leg_small)
  df_a <- analyzed$leg_small
  
  p1 <- ggplot(df_a, aes(x = n, y = median_speedup, color = expr, fill = expr)) +
    geom_ribbon(aes(ymin = ci_low, ymax = ci_high), alpha = 0.2, color = NA) +
    geom_line(linewidth = 0.8) + geom_point(size = 1.5) +
    scale_x_log10(
      breaks = 10^(0:7),
      labels = trans_format("log10", math_format(10^.x))
    ) +
    scale_y_continuous(labels = label_number(suffix = "x")) +
    expand_limits(y = 0) +
    labs(
      title = "Panel A: Small-sample M-estimators",
      subtitle = "robscale vs revss (with 95% BCa confidence bands)",
      x = "Sample Size (n)",
      y = "Median Speedup Factor",
      color = "Estimator",
      fill = "Estimator"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
    
  # Panel B: Scale estimators vs robustbase (analyzed$leg_large)
  df_b <- analyzed$leg_large
    
  p2 <- ggplot(df_b, aes(x = n, y = median_speedup, color = expr, fill = expr)) +
    geom_ribbon(aes(ymin = ci_low, ymax = ci_high), alpha = 0.2, color = NA) +
    geom_line(linewidth = 0.8) + geom_point(size = 1.5) +
    scale_x_log10(
      breaks = 10^(0:7),
      labels = trans_format("log10", math_format(10^.x))
    ) +
    scale_y_continuous(labels = label_number(suffix = "x")) +
    expand_limits(y = 0) +
    labs(
      title = "Panel B: Scale Estimators",
      subtitle = "robscale vs robustbase (with 95% BCa confidence bands)",
      x = "Sample Size (n)",
      y = "Median Speedup Factor",
      color = "Estimator",
      fill = "Estimator"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
    
  # Panel C: New scale estimators vs existing R implementations
  df_c <- analyzed$new_estimators

  p3 <- ggplot(df_c, aes(x = n, y = median_speedup, color = expr, fill = expr)) +
    geom_ribbon(aes(ymin = ci_low, ymax = ci_high), alpha = 0.2, color = NA) +
    geom_line(linewidth = 0.8) + geom_point(size = 1.5) +
    scale_x_log10(
      breaks = 10^(0:7),
      labels = trans_format("log10", math_format(10^.x))
    ) +
    scale_y_continuous(labels = label_number(suffix = "x")) +
    expand_limits(y = 0) +
    labs(
      title = "Panel C: New Scale Estimators",
      subtitle = "robscale vs existing R implementations (with 95% BCa confidence bands)",
      x = "Sample Size (n)",
      y = "Median Speedup Factor",
      color = "Comparison",
      fill = "Comparison"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")

  p1 + p2 + p3 + plot_layout(guides = "collect") & theme(legend.position = "bottom")
}
