library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(scales)

#' Create the speedup plot
#' @param analyzed The analyzed results from analyze_benchmarks()
plot_benchmarks <- function(analyzed) {
  # Parity reference line shared by all three panels
  parity_line <- geom_hline(yintercept = 1, linewidth = 0.3,
                             linetype = "solid", color = "grey50")

  # Panel A: M-estimators vs revss (analyzed$leg_small)
  df_a <- analyzed$leg_small

  p1 <- ggplot(df_a, aes(x = n, y = median_speedup, color = expr, fill = expr)) +
    parity_line +
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
    parity_line +
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
    parity_line +
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

#' Create the absolute-timing figure
#'
#' Combines all robscale estimators into a single panel with log-log axes.
#' The y-axis shows median wall-clock run time; the x-axis shows sample size.
#'
#' @param rob The raw robscale benchmark results list (rob_optimized).
plot_absolute_timings <- function(rob) {
  # bench_time median: as.numeric() returns seconds (bench stores ns but
  # as.numeric.bench_time divides by 1e9 before returning).
  to_secs <- function(df) {
    df %>%
      mutate(
        n          = as.numeric(unlist(n)),
        expr       = as.character(expression),
        median_sec = as.numeric(median)
      ) %>%
      select(n, expr, median_sec)
  }

  df <- bind_rows(
    to_secs(rob$m_estimators),
    to_secs(rob$scale_estimators),
    to_secs(rob$new_estimators)
  )

  # Human-readable time labels for a log-scale axis (input: seconds)
  fmt_time <- function(x) {
    dplyr::case_when(
      x <  1e-6 ~ paste0(signif(x * 1e9, 3), " ns"),
      x <  1e-3 ~ paste0(signif(x * 1e6, 3), " \u03bcs"),
      x <  1    ~ paste0(signif(x * 1e3, 3), " ms"),
      TRUE      ~ paste0(signif(x,         3), " s")
    )
  }

  ggplot(df, aes(x = n, y = median_sec, color = expr)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.5) +
    scale_x_log10(
      breaks = 10^(0:7),
      labels = trans_format("log10", math_format(10^.x))
    ) +
    scale_y_log10(labels = fmt_time) +
    labs(
      title    = "Absolute run times \u2014 all robscale estimators",
      subtitle = "Median wall-clock time across benchmark sample sizes",
      x        = "Sample size (n)",
      y        = "Median run time",
      color    = "Estimator"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
}
