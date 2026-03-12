library(robscale)
library(bench)
library(dplyr)
library(tidyr)
library(purrr)

# Grid from run_benchmarks.R
n_small <- c(3, 4, 5, 6, 7, 8, 10, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 12288, 16384)
n_large <- c(3, 4, 5, 6, 7, 8, 10, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 
             12288, 16384, 32768, 65536, 131072, 262144, 524288, 1048576, 10000000)

get_min_iters <- function(n) {
  if (n <= 128) 10000L
  else if (n <= 2048) 2000L
  else if (n <= 16384) 500L
  else if (n <= 1048576) 20L
  else 5L
}

run_full_baseline <- function() {
  cat("Running full baseline for m-estimators...\n")
  m_res <- bench::press(
    n = n_small,
    {
      set.seed(42 + n)
      x <- rnorm(n)
      bench::mark(
        robLoc = robscale::robLoc(x),
        robScale = robscale::robScale(x),
        adm = robscale::adm(x),
        check = FALSE,
        min_iterations = get_min_iters(n),
        min_time = 1.0
      )
    }
  ) %>% mutate(category = "m-estimators")
  
  cat("Running full baseline for scale estimators...\n")
  scale_res <- bench::press(
    n = n_large,
    {
      set.seed(42 + n)
      x <- rnorm(n)
      bench::mark(
        qn = robscale::qn(x),
        sn = robscale::sn(x),
        check = FALSE,
        min_iterations = get_min_iters(n),
        min_time = 1.0
      )
    }
  ) %>% mutate(category = "scale-estimators")
  
  full_res <- bind_rows(m_res, scale_res)
  saveRDS(full_res, "benchmarks/golden_baseline_0_1_5.rds")
  cat("Baseline saved to benchmarks/golden_baseline_0_1_5.rds\n")
}

run_full_baseline()
