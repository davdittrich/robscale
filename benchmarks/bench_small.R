suppressPackageStartupMessages({
  library(robscale)
  library(robustbase)
  library(fastqnsn)
  library(revss)
  library(microbenchmark)
  library(dplyr)
  library(boot)
  library(ggplot2)
})

comparisons <- list(
  list(name = "Sn", new = robscale::sn, legacy = list(list(name = "robustbase::Sn", fn = robustbase::Sn), list(name = "fastqnsn::sn", fn = fastqnsn::sn)), sizes = 3:7),
  list(name = "Qn", new = robscale::qn, legacy = list(list(name = "robustbase::Qn", fn = robustbase::Qn), list(name = "fastqnsn::qn", fn = fastqnsn::qn)), sizes = 3:7),
  list(name = "robScale", new = robscale::robScale, legacy = list(list(name = "revss::robScale", fn = revss::robScale)), sizes = 3:7),
  list(name = "robLoc", new = robscale::robLoc, legacy = list(list(name = "revss::robLoc", fn = revss::robLoc)), sizes = 3:7),
  list(name = "adm", new = robscale::adm, legacy = list(list(name = "revss::adm", fn = revss::adm)), sizes = 3:7)
)

all_raw <- list()
for (comp in comparisons) {
  for (n in comp$sizes) {
    set.seed(42)
    x <- rnorm(n)
    exprs <- list(robscale = substitute(new_fn(x), list(new_fn = comp$new)))
    for (leg in comp$legacy) exprs[[leg$name]] <- substitute(leg_fn(x), list(leg_fn = leg$fn))
    bm <- microbenchmark(list = exprs, times = 20000L)
    df <- as.data.frame(bm) %>% mutate(n = n, estimator = comp$name) %>% rename(implementation = expr, time_ns = time) %>% mutate(time_us = time_ns / 1000)
    all_raw[[paste(comp$name, n, sep="_")]] <- df
  }
}
raw_small <- bind_rows(all_raw)
write.csv(raw_small, "optimized_mode/small_n_raw.csv", row.names = FALSE)

# Load existing data, drop old 3:7 if existent, compute new summary
raw_existing <- read.csv("optimized_mode/legacy_benchmark_raw.csv") %>% filter(!n %in% 3:7)
raw_full <- bind_rows(raw_existing, raw_small)
write.csv(raw_full, "optimized_mode/legacy_benchmark_raw.csv", row.names = FALSE)

# Recompute summary for everything
ratio_median <- function(d, i) { median(d[i, "legacy"]) / median(d[i, "robscale"]) }
summary_list <- list()
combos <- raw_full %>% filter(implementation != "robscale") %>% distinct(estimator, implementation, n)
for (j in 1:nrow(combos)) {
  row <- combos[j,]
  est <- row$estimator; leg_name <- row$implementation; n_val <- row$n
  rob_times <- raw_full %>% filter(estimator == est, n == n_val, implementation == "robscale") %>% pull(time_us)
  leg_times <- raw_full %>% filter(estimator == est, n == n_val, implementation == leg_name) %>% pull(time_us)
  len <- min(length(rob_times), length(leg_times))
  boot_data <- data.frame(robscale = rob_times[1:len], legacy = leg_times[1:len])
  b <- boot(boot_data, ratio_median, R = 500)
  ci <- tryCatch(boot.ci(b, type = "perc")$percent[4:5], error = function(e) c(NA, NA))
  summary_list[[j]] <- tibble(
    n = n_val, estimator = est, legacy_name = as.character(leg_name),
    robscale_med_us = median(rob_times), legacy_med_us = median(leg_times),
    speedup = median(leg_times) / median(rob_times), ci_lo = ci[1], ci_hi = ci[2]
  )
}
summary_full <- bind_rows(summary_list)
write.csv(summary_full, "optimized_mode/legacy_benchmark_summary.csv", row.names = FALSE)
cat("Successfully updated raw and summary CSVs with n=3 to 7.\n")
