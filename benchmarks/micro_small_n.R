library(bench)
library(callr)

# Paths
GOLD_SRC <- "/tmp/robscale_clean_015"
LIB_GOLD <- file.path(tempdir(), "lib_gold")
dir.create(LIB_GOLD, showWarnings = FALSE)

# Install Gold Standard 0.1.5 into temp lib
if (!file.exists(file.path(LIB_GOLD, "robscale"))) {
  cat("Building Gold Standard v0.1.5...\n")
  rcmd("INSTALL", c(GOLD_SRC, paste0("--library=", LIB_GOLD)))
}

library(robscale)

# Define sample sizes
ns <- c(3, 4, 5, 6, 7, 8, 10, 16, 32, 64)

results <- lapply(ns, function(n) {
  x <- rnorm(n)
  
  cat(sprintf("\n--- n = %d ---\n", n))
  
  # Benchmark adm
  bm_adm <- mark(
    opt = robscale::adm(x),
    gold = callr::r(function(x) robscale::adm(x), args = list(x = x), libpath = LIB_GOLD),
    iterations = 1000,
    check = FALSE
  )
  print(bm_adm)
  
  # Benchmark Sn
  bm_sn <- mark(
    opt = robscale::Sn(x),
    gold = callr::r(function(x) robscale::Sn(x), args = list(x = x), libpath = LIB_GOLD),
    iterations = 1000,
    check = FALSE
  )
  
  list(n = n, adm = bm_adm, sn = bm_sn)
})

saveRDS(results, "benchmarks/micro_small_n_results.rds")
