# benchmarks/do_verification.R
source("benchmarks/rigorous_verify.R")

# 1. Build libs
lib_fast <- build_robscale(fast_mode = 1)
lib_slow <- build_robscale(fast_mode = 0)

# Save for reference if needed
writeLines(c(lib_fast, lib_slow), "benchmarks/verify_libs.txt")

# 2. Run verification
# We only check the small sizes where regressions were seen
n_test <- c(3, 4, 16, 64)
estimators <- c("robLoc", "robScale", "adm", "qn", "sn")

results <- list()
for (est in estimators) {
  for (n in n_test) {
    res <- run_benchmark(lib_fast, lib_slow, est, n, min_time = 1.0)
    results[[length(results) + 1]] <- res
    # Print individual result immediately
    print(res)
  }
}

final_results <- bind_rows(results)
saveRDS(final_results, "benchmarks/final_verification_results.rds")

# Summary Table
cat("\n--- FINAL VERIFICATION SUMMARY ---\n")
print(final_results %>% select(estimator, n, speedup, ci_low, status))

if (any(final_results$status == "❌ FAIL")) {
  cat("\nWARNING: Some regressions still exceed 1% threshold!\n")
} else {
  cat("\nSUCCESS: All estimators meet performance requirements.\n")
}
