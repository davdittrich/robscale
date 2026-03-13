library(dplyr)
library(ggplot2)
library(scales)
library(patchwork)

# Load optimized mode summary
df <- read.csv("optimized_mode/legacy_benchmark_summary.csv")

# Ensure factors for consistent ordering
df$estimator <- factor(df$estimator, levels = c("robLoc", "robScale", "adm", "Qn", "Sn"))

# Left panel: M-estimators vs revss
df_left <- df %>% filter(grepl("revss", legacy_name))
p1 <- ggplot(df_left, aes(x = n, y = speedup, color = estimator, fill = estimator)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
  geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1) +
  geom_point() +
  facet_wrap(~estimator, scales = "free_y", ncol = 1) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_continuous(labels = label_number(suffix = "x"), limits = c(0, NA)) +
  labs(title = "M-estimators vs. revss",
       x = "Sample Size (n)", y = "Speedup Factor",
       color = "Estimator", fill = "Estimator") +
  theme_minimal() +
  theme(legend.position = "none")

# Right panel: Scale estimators vs robustbase
df_right <- df %>% filter(grepl("robustbase", legacy_name))
p2 <- ggplot(df_right, aes(x = n, y = speedup, color = estimator, fill = estimator)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
  geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1) +
  geom_point() +
  facet_wrap(~estimator, scales = "free_y", ncol = 1) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_continuous(labels = label_number(suffix = "x"), limits = c(0, NA)) +
  labs(title = "Scale estimators vs. robustbase",
       x = "Sample Size (n)", y = "Speedup Factor",
       color = "Estimator", fill = "Estimator") +
  theme_minimal() +
  theme(legend.position = "none")

# Combine with patchwork
p_combined <- p1 + p2 + plot_annotation(
  title = "robscale Speedup over Legacy Implementations",
  subtitle = "Speedup = Median(Legacy) / Median(robscale). Ribbon is 95% bootstrap CI."
)

ggsave("optimized_mode/full_speedup_panel.png", p_combined, width = 10, height = 7, dpi = 300)
cat("Plot saved to optimized_mode/full_speedup_panel.png\n")
