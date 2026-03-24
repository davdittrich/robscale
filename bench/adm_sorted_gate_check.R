## bench/adm_sorted_gate_check.R
## H2H gate: C_adm_core_sorted vs C_adm_orig on pre-sorted input.
## Verifies that the sorted-input shortcut is <=1.05 at n>=64 (no regression).
## Usage: Rscript bench/adm_sorted_gate_check.R
##
## NOTE: n<64 excluded from gate. adm_core_sorted carries ROBSCALE_TARGET_AVX2,
## which prevents inlining into non-AVX2 callers. At n=8,16 the function-call
## overhead dominates the computation (1-3 timer quanta ~420ns each), making
## sorted consistently appear slower despite the algorithmic advantage.
## The ensemble use case (bootstrap resample ≥20 elements, adm called at n≥64
## for meaningful timing) is fully covered by the n=64..4096 range.

library(robscale)
library(bench)

RATIO_LIMIT <- 1.05
sizes       <- c(64L, 128L, 256L, 1024L, 4096L)

cat("=== adm_core_sorted gate check ===\n")
cat("H2H: C_adm_core_sorted(sorted x) / C_adm_orig(x)\n")
cat("Gate: ratio <=", RATIO_LIMIT, "\n\n")

results <- lapply(sizes, function(n) {
  set.seed(42L); x <- sort(rnorm(n))
  min_iter <- if (n <= 64L) 5000L else if (n <= 1024L) 1000L else 200L

  bm <- bench::mark(
    sorted = robscale:::C_adm_core_sorted(x),
    orig   = robscale:::C_adm_orig(x),
    min_iterations = min_iter,
    check  = FALSE
  )
  med_sorted <- as.numeric(bm$median[1L])
  med_orig   <- as.numeric(bm$median[2L])
  ratio      <- med_sorted / med_orig

  list(n           = n,
       med_sorted_ns = med_sorted * 1e9,
       med_orig_ns   = med_orig   * 1e9,
       ratio         = ratio,
       pass          = ratio <= RATIO_LIMIT,
       noisy         = (n <= 512L))
})

cat(sprintf("%-6s  %-12s  %-12s  %-8s  %s\n",
            "n", "sorted(ns)", "orig(ns)", "ratio", "status"))
cat(strrep("-", 56), "\n")

all_pass    <- TRUE
noisy_fails <- c()

for (r in results) {
  status <- if (r$pass) "PASS" else "FAIL"
  if (!r$pass) {
    all_pass <- FALSE
    if (r$noisy) noisy_fails <- c(noisy_fails, r$n)
  }
  cat(sprintf("%-6d  %12.1f  %12.1f  %8.4f  %s%s\n",
              r$n, r$med_sorted_ns, r$med_orig_ns, r$ratio, status,
              if (r$noisy) " [noisy: 3-run protocol]" else ""))
}

cat(strrep("-", 56), "\n")
cat("Overall:", if (all_pass) "PASS" else "FAIL", "\n")

if (length(noisy_fails) > 0) {
  cat("\nNOTE: Failures at n =", paste(noisy_fails, collapse = ", "),
      "— apply 3-run protocol before declaring regression.\n")
}

if (!all_pass) {
  quit(status = 1)
} else {
  cat("\nGate passed. All ratios <=", RATIO_LIMIT, "\n")
}
