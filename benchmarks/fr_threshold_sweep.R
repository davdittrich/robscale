#!/usr/bin/env Rscript
# P3 calibration: sweep Floyd-Rivest vs std::nth_element crossover threshold.
# Self-contained — embeds FR algorithm in Rcpp::sourceCpp() for portability.

library(Rcpp)

cpp_code <- '
#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <chrono>
#include <random>
#include <vector>

// Pure std::nth_element (no FR narrowing)
static void nth_element_only(double* x, size_t n, size_t k) {
  std::nth_element(x, x + k, x + n);
}

// Floyd-Rivest with narrowing (always applies narrowing, no fallback)
static void fr_narrowing(double* left_p, double* k_p, double* right_p) {
  double* left = left_p;
  double* right = right_p - 1;

  while (right > left) {
    if (right - left > 0) {  // always narrow
      size_t nn = static_cast<size_t>(right - left + 1);
      size_t i = static_cast<size_t>(k_p - left + 1);
      double z = std::log(static_cast<double>(nn));
      double s = 0.5 * std::exp(2.0 * z / 3.0);
      double sd = 0.5 *
          std::sqrt(z * s * (static_cast<double>(nn) - s) /
                    static_cast<double>(nn)) *
          (static_cast<double>(i) - static_cast<double>(nn) / 2.0 >= 0 ? 1.0 : -1.0);
      double* new_left = std::max(left,
          k_p - static_cast<ptrdiff_t>(static_cast<double>(i) * s / nn + sd));
      double* new_right = std::min(right,
          k_p + static_cast<ptrdiff_t>(static_cast<double>(nn - i) * s / nn + sd));
      // Recurse on narrowed window
      fr_narrowing(new_left, k_p, new_right + 1);
    }

    double pivot = *k_p;
    double* ii = left;
    double* jj = right;
    std::swap(*left, *k_p);
    if (*right > pivot) std::swap(*left, *right);

    while (ii < jj) {
      std::swap(*ii, *jj);
      ++ii; --jj;
      while (*ii < pivot) ++ii;
      while (*jj > pivot) --jj;
    }

    if (*left == pivot) {
      std::swap(*left, *jj);
    } else {
      ++jj;
      std::swap(*jj, *right);
    }

    if (jj <= k_p) left = jj + 1;
    if (k_p <= jj) right = jj - 1;
  }
}

// The hybrid: uses std::nth_element below threshold, FR above
static void fr_hybrid(double* x, size_t n, size_t k, int threshold) {
  if ((int)n < threshold) {
    std::nth_element(x, x + k, x + n);
    return;
  }
  double* left = x;
  double* right = x + n - 1;
  double* kp = x + k;

  while (right > left) {
    if (right - left > threshold) {
      size_t nn = static_cast<size_t>(right - left + 1);
      size_t i = static_cast<size_t>(kp - left + 1);
      double z = std::log(static_cast<double>(nn));
      double s = 0.5 * std::exp(2.0 * z / 3.0);
      double sd = 0.5 *
          std::sqrt(z * s * (static_cast<double>(nn) - s) /
                    static_cast<double>(nn)) *
          (static_cast<double>(i) - static_cast<double>(nn) / 2.0 >= 0 ? 1.0 : -1.0);
      double* new_left = std::max(left,
          kp - static_cast<ptrdiff_t>(static_cast<double>(i) * s / nn + sd));
      double* new_right = std::min(right,
          kp + static_cast<ptrdiff_t>(static_cast<double>(nn - i) * s / nn + sd));
      // Inline the recursion for the narrowed window
      std::nth_element(new_left, kp, new_right + 1);
    }

    double pivot = *kp;
    double* ii = left;
    double* jj = right;
    std::swap(*left, *kp);
    if (*right > pivot) std::swap(*left, *right);

    while (ii < jj) {
      std::swap(*ii, *jj);
      ++ii; --jj;
      while (*ii < pivot) ++ii;
      while (*jj > pivot) --jj;
    }

    if (*left == pivot) {
      std::swap(*left, *jj);
    } else {
      ++jj;
      std::swap(*jj, *right);
    }

    if (jj <= kp) left = jj + 1;
    if (kp <= jj) right = jj - 1;
  }
}

// [[Rcpp::export]]
Rcpp::DataFrame sweep_fr_threshold(int n, Rcpp::IntegerVector thresholds,
                                    int iters, int seed) {
  std::vector<double> master(n);
  std::mt19937 rng(seed);
  std::uniform_real_distribution<double> dist(0.0, 1.0);
  for (int i = 0; i < n; ++i) master[i] = dist(rng);

  size_t k = n / 2;  // median position
  std::vector<double> work(n);

  int nt = thresholds.size();
  Rcpp::IntegerVector thresh_out(nt + 1);
  Rcpp::NumericVector time_out(nt + 1);

  // Baseline: pure std::nth_element (threshold = infinity)
  {
    double total = 0.0;
    for (int iter = 0; iter < iters; ++iter) {
      std::memcpy(work.data(), master.data(), n * sizeof(double));
      auto t0 = std::chrono::high_resolution_clock::now();
      nth_element_only(work.data(), n, k);
      auto t1 = std::chrono::high_resolution_clock::now();
      total += std::chrono::duration<double, std::nano>(t1 - t0).count();
    }
    thresh_out[0] = 999999;
    time_out[0] = total / iters;
  }

  // Sweep thresholds
  for (int ti = 0; ti < nt; ++ti) {
    int thresh = thresholds[ti];
    double total = 0.0;
    for (int iter = 0; iter < iters; ++iter) {
      std::memcpy(work.data(), master.data(), n * sizeof(double));
      auto t0 = std::chrono::high_resolution_clock::now();
      fr_hybrid(work.data(), n, k, thresh);
      auto t1 = std::chrono::high_resolution_clock::now();
      total += std::chrono::duration<double, std::nano>(t1 - t0).count();
    }
    thresh_out[ti + 1] = thresh;
    time_out[ti + 1] = total / iters;
  }

  return Rcpp::DataFrame::create(
    Rcpp::Named("threshold") = thresh_out,
    Rcpp::Named("time_ns") = time_out
  );
}
'

cat("Compiling benchmark harness...\n")
sourceCpp(code = cpp_code, verbose = FALSE)

# Configuration
sample_sizes <- c(200, 400, 600, 800, 1000, 1500, 2000, 4000, 8000, 16000)
thresholds <- as.integer(c(100, 200, 300, 400, 500, 600, 700, 800, 1000, 1200, 1500, 2000))
seed <- 42L

get_iters <- function(n) {
  if (n <= 200)  return(50000L)
  if (n <= 1000) return(20000L)
  if (n <= 4000) return(5000L)
  return(2000L)
}

cat(sprintf("\n  Platform: %s (%s)\n", Sys.info()["sysname"], Sys.info()["machine"]))
cat(sprintf("  R: %s\n\n", R.version.string))

# Run sweep for each sample size
results <- list()
for (n in sample_sizes) {
  iters <- get_iters(n)
  cat(sprintf("  n = %-6d (%d iters)...", n, iters))
  df <- sweep_fr_threshold(n, thresholds, iters, seed)
  df$n <- n

  # Find the fastest threshold for this n
  best_idx <- which.min(df$time_ns)
  nth_time <- df$time_ns[df$threshold == 999999]
  best_time <- df$time_ns[best_idx]
  best_thresh <- df$threshold[best_idx]

  cat(sprintf(" nth_element: %.0fns, best: %.0fns (thresh=%s, %.1fx)\n",
              nth_time, best_time,
              if (best_thresh == 999999) "nth_only" else as.character(best_thresh),
              nth_time / best_time))

  results[[length(results) + 1]] <- df
}

all_results <- do.call(rbind, results)

# Summary: for each n, show the optimal threshold
cat("\n  === Optimal threshold per sample size ===\n")
cat(sprintf("  %-8s %-12s %-12s %-10s %-10s\n",
            "n", "nth_el(ns)", "best(ns)", "threshold", "speedup"))
cat("  ", strrep("-", 56), "\n", sep = "")

for (n in sample_sizes) {
  sub <- all_results[all_results$n == n, ]
  nth_time <- sub$time_ns[sub$threshold == 999999]
  best_idx <- which.min(sub$time_ns)
  best_time <- sub$time_ns[best_idx]
  best_thresh <- sub$threshold[best_idx]

  cat(sprintf("  %-8d %-12.0f %-12.0f %-10s %-10s\n",
              n, nth_time, best_time,
              if (best_thresh == 999999) "nth_only" else as.character(best_thresh),
              sprintf("%.2fx", nth_time / best_time)))
}

# Show full heatmap: time relative to nth_element for each (n, threshold)
cat("\n  === Heatmap: time relative to std::nth_element (lower is better) ===\n")
cat(sprintf("  %-8s", "n"))
for (t in thresholds) cat(sprintf(" t=%-5d", t))
cat("\n")
cat("  ", strrep("-", 8 + length(thresholds) * 8), "\n", sep = "")

for (n in sample_sizes) {
  sub <- all_results[all_results$n == n, ]
  nth_time <- sub$time_ns[sub$threshold == 999999]
  cat(sprintf("  %-8d", n))
  for (t in thresholds) {
    t_time <- sub$time_ns[sub$threshold == t]
    if (length(t_time) == 1) {
      ratio <- t_time / nth_time
      cat(sprintf(" %6.2f ", ratio))
    } else {
      cat("   --   ")
    }
  }
  cat("\n")
}
