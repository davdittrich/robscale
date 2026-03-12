#include <Rcpp.h>
#include <algorithm>
#include <chrono>
#include <cstring>
#include <numeric>
#include <random>
#include <vector>
#include "sort_alternatives.h"

// --------------------------------------------------------------------------
// Method registry
// --------------------------------------------------------------------------
using SortFn = void (*)(double*, size_t);

static void fn_knuth(double* x, size_t n)     { robscale::sort_alt::knuth_network(x, n); }
static void fn_std(double* x, size_t n)        { robscale::sort_alt::std_sort(x, n); }
static void fn_insertion(double* x, size_t n)  { robscale::sort_alt::insertion_sort(x, n); }
static void fn_simd(double* x, size_t n)       { robscale::sort_alt::simd_sort(x, n); }

struct Method {
  const char* name;
  SortFn fn;
};

static const Method methods[] = {
  {"knuth_network",  fn_knuth},
  {"std_sort",       fn_std},
  {"insertion_sort", fn_insertion},
  {"simd_sort",      fn_simd},
};
static constexpr int n_methods = sizeof(methods) / sizeof(methods[0]);

// --------------------------------------------------------------------------
// Test data generators
// --------------------------------------------------------------------------
static void gen_random(double* x, size_t n, std::mt19937& rng) {
  std::normal_distribution<double> dist(0.0, 1.0);
  for (size_t i = 0; i < n; ++i) x[i] = dist(rng);
}
static void gen_sorted(double* x, size_t n, std::mt19937& rng) {
  gen_random(x, n, rng);
  std::sort(x, x + n);
}
static void gen_reverse(double* x, size_t n, std::mt19937& rng) {
  gen_random(x, n, rng);
  std::sort(x, x + n, std::greater<double>());
}
static void gen_equal(double* x, size_t n, std::mt19937&) {
  std::fill(x, x + n, 42.0);
}
static void gen_two_vals(double* x, size_t n, std::mt19937& rng) {
  std::bernoulli_distribution coin(0.5);
  for (size_t i = 0; i < n; ++i) x[i] = coin(rng) ? 1.0 : -1.0;
}
static void gen_outlier(double* x, size_t n, std::mt19937& rng) {
  for (size_t i = 0; i < n; ++i) x[i] = static_cast<double>(i);
  x[n / 2] = 1e6;
  // Shuffle
  for (size_t i = n - 1; i > 0; --i) {
    std::uniform_int_distribution<size_t> d(0, i);
    std::swap(x[i], x[d(rng)]);
  }
}

using GenFn = void (*)(double*, size_t, std::mt19937&);
struct Pattern {
  const char* name;
  GenFn fn;
};
static const Pattern patterns[] = {
  {"random",      gen_random},
  {"sorted",      gen_sorted},
  {"reverse",     gen_reverse},
  {"equal",       gen_equal},
  {"two_vals",    gen_two_vals},
  {"outlier",     gen_outlier},
};
static constexpr int n_patterns = sizeof(patterns) / sizeof(patterns[0]);

// --------------------------------------------------------------------------
// Correctness test
// --------------------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::DataFrame sort_shootout_correctness() {
  static const int sizes[] = {2, 3, 4, 5, 6, 7, 8};
  static constexpr int n_sizes = sizeof(sizes) / sizeof(sizes[0]);
  static constexpr int n_seeds = 50;

  int total_rows = n_methods * n_sizes * n_patterns * n_seeds;
  Rcpp::StringVector col_method(total_rows);
  Rcpp::IntegerVector col_n(total_rows);
  Rcpp::StringVector col_pattern(total_rows);
  Rcpp::IntegerVector col_seed(total_rows);
  Rcpp::LogicalVector col_match(total_rows);

  int row = 0;
  std::vector<double> input(8), work(8), ref(8);

  for (int si = 0; si < n_sizes; ++si) {
    int n = sizes[si];
    for (int pi = 0; pi < n_patterns; ++pi) {
      for (int seed = 1; seed <= n_seeds; ++seed) {
        // Generate input
        std::mt19937 rng(static_cast<unsigned>(seed + n * 1000 + pi * 100));
        patterns[pi].fn(input.data(), static_cast<size_t>(n), rng);

        // Reference sort
        std::memcpy(ref.data(), input.data(), n * sizeof(double));
        std::sort(ref.data(), ref.data() + n);

        for (int mi = 0; mi < n_methods; ++mi) {
          std::memcpy(work.data(), input.data(), n * sizeof(double));
          methods[mi].fn(work.data(), static_cast<size_t>(n));

          bool match = std::equal(work.data(), work.data() + n, ref.data());
          col_method[row] = methods[mi].name;
          col_n[row] = n;
          col_pattern[row] = patterns[pi].name;
          col_seed[row] = seed;
          col_match[row] = match;
          ++row;
        }
      }
    }
  }

  return Rcpp::DataFrame::create(
    Rcpp::Named("method") = col_method,
    Rcpp::Named("n") = col_n,
    Rcpp::Named("pattern") = col_pattern,
    Rcpp::Named("seed") = col_seed,
    Rcpp::Named("matches_reference") = col_match,
    Rcpp::Named("stringsAsFactors") = false
  );
}

// --------------------------------------------------------------------------
// Benchmark — batch timing for sub-nanosecond resolution
// --------------------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::DataFrame sort_shootout_benchmark(int rounds = 200, int batch = 1000) {
  static const int sizes[] = {2, 3, 4, 5, 6, 7, 8};
  static constexpr int n_sizes = sizeof(sizes) / sizeof(sizes[0]);

  int total_rows = n_methods * n_sizes;
  Rcpp::StringVector col_method(total_rows);
  Rcpp::IntegerVector col_n(total_rows);
  Rcpp::NumericVector col_median_ns(total_rows);
  Rcpp::NumericVector col_mean_ns(total_rows);
  Rcpp::NumericVector col_min_ns(total_rows);

  // Pre-generate random data: batch arrays per size, each array = n doubles
  std::mt19937 rng(12345);
  std::vector<std::vector<double>> data_sets(n_sizes);
  for (int si = 0; si < n_sizes; ++si) {
    int n = sizes[si];
    data_sets[si].resize(static_cast<size_t>(n) * batch);
    std::normal_distribution<double> dist(0.0, 1.0);
    for (int b = 0; b < batch; ++b) {
      for (int j = 0; j < n; ++j) {
        data_sets[si][static_cast<size_t>(b) * n + j] = dist(rng);
      }
    }
  }

  // Working copy: one contiguous block per batch
  std::vector<double> work;
  std::vector<double> timings(rounds);

  int row = 0;
  for (int si = 0; si < n_sizes; ++si) {
    int n = sizes[si];
    size_t block = static_cast<size_t>(n) * batch;
    work.resize(block);

    for (int mi = 0; mi < n_methods; ++mi) {
      // Warm up
      std::memcpy(work.data(), data_sets[si].data(), block * sizeof(double));
      for (int b = 0; b < batch; ++b) {
        methods[mi].fn(work.data() + static_cast<size_t>(b) * n,
                       static_cast<size_t>(n));
      }

      // Timed rounds: each round sorts the full batch from fresh copies
      for (int r = 0; r < rounds; ++r) {
        std::memcpy(work.data(), data_sets[si].data(), block * sizeof(double));
        auto t0 = std::chrono::high_resolution_clock::now();
        for (int b = 0; b < batch; ++b) {
          methods[mi].fn(work.data() + static_cast<size_t>(b) * n,
                         static_cast<size_t>(n));
        }
        auto t1 = std::chrono::high_resolution_clock::now();
        double elapsed_ns = std::chrono::duration<double, std::nano>(t1 - t0).count();
        timings[r] = elapsed_ns / batch;
      }

      std::sort(timings.begin(), timings.end());
      double median_ns = timings[rounds / 2];
      double min_ns = timings[0];
      double sum = std::accumulate(timings.begin(), timings.end(), 0.0);
      double mean_ns = sum / rounds;

      col_method[row] = methods[mi].name;
      col_n[row] = n;
      col_median_ns[row] = median_ns;
      col_mean_ns[row] = mean_ns;
      col_min_ns[row] = min_ns;
      ++row;
    }
  }

  return Rcpp::DataFrame::create(
    Rcpp::Named("method") = col_method,
    Rcpp::Named("n") = col_n,
    Rcpp::Named("median_ns") = col_median_ns,
    Rcpp::Named("mean_ns") = col_mean_ns,
    Rcpp::Named("min_ns") = col_min_ns,
    Rcpp::Named("stringsAsFactors") = false
  );
}
