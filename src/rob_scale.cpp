#include "robscale_config.h"
#include "robust_core.h"
#include "pdq_select.h"
#include <Rcpp.h>
#include <memory>
#if defined(ROBSCALE_HAS_SYSTEM_TBB)
#  include <oneapi/tbb/parallel_reduce.h>
#  include <oneapi/tbb/blocked_range.h>
#elif defined(USE_DIRECT_TBB)
#  include <tbb/parallel_reduce.h>
#  include <tbb/blocked_range.h>
#endif

// ---------------------------------------------------------------------------
// Phase 2: Fused per-iteration kernel
//
// The 3-pass loop (scale → tanh → sum-of-squares) writes and reads tmp[]
// 4 times per element per iteration.  For large n this dominates: at n=10^6
// with 32 iterations that is 1 GB of tmp traffic.
//
// The fused kernel collapses all three passes into one, reading data[] once
// and never touching tmp[].  Uses ROBSCALE_TANH4_AVX2 (libmvec preferred,
// SLEEF fallback) + 4-wide FMA accumulation.
//
// D-1: Aitken Δ² acceleration
//
// rob_scale_compute and rob_scale_parallel_compute use Steffensen's method:
// every 2 standard steps, extrapolate to s_acc = s2-(s2-s1)²/(s2-2s1+s0).
// Convergence still exits via the standard |v-1|≤tol check.  This reduces
// iteration count by ~30-50% for small n and ~20% for large n.
// ---------------------------------------------------------------------------

#if defined(ROBSCALE_HAS_SLEEF) && !defined(ROBSCALE_HAS_ACCELERATE) && \
    defined(ROBSCALE_HAS_AVX2_DISPATCH)
/**
 * Fused sum: compute sum_i tanh((data[i]-off)*hinv)^2 in one pass.
 *
 * Uses ROBSCALE_TANH4_AVX2 (glibc libmvec _ZGVdN4v_tanh when available,
 * SLEEF Sleef_tanhd4_u10avx2 otherwise) + FMA accumulation into 4 independent
 * SIMD lanes, avoiding the serial scalar-extraction dependency chain that would
 * otherwise serialise the accumulation loop.  A horizontal sum at the end
 * collapses the 4 lanes.  FP rounding differs from the 3-pass scalar sum
 * by at most n*eps*sum_rho ≈ 1e-12 for the golden test sizes (n ≤ 20),
 * safely within the 1e-12 tolerance of test-cross-check.R.
 *
 * Target attribute: AVX2+FMA instructions are emitted here only; the caller
 * (rob_scale_compute) stays at the baseline ISA, keeping non-AVX2 binaries
 * safe via the runtime-guarded call site.
 */
ROBSCALE_TARGET_AVX2
static double rob_scale_fused_sum_avx2(const double* ROBSCALE_RESTRICT data,
                                        int n, double data_offset,
                                        double half_inv_sc) {
  const __m256d off4  = _mm256_set1_pd(data_offset);
  const __m256d hinv4 = _mm256_set1_pd(half_inv_sc);
  __m256d acc = _mm256_setzero_pd();
  int i = 0;
  for (; i + 4 <= n; i += 4) {
    __m256d d = _mm256_loadu_pd(data + i);
    __m256d t = ROBSCALE_TANH4_AVX2(
        _mm256_mul_pd(_mm256_sub_pd(d, off4), hinv4));
    acc = _mm256_fmadd_pd(t, t, acc);  // acc[lane] += t[lane]^2, 4-wide
  }
  // Horizontal sum: [a0,a1,a2,a3] → scalar
  __m128d lo   = _mm256_castpd256_pd128(acc);       // [a0, a1]
  __m128d hi   = _mm256_extractf128_pd(acc, 1);     // [a2, a3]
  __m128d sum2 = _mm_add_pd(lo, hi);                // [a0+a2, a1+a3]
  __m128d sum1 = _mm_hadd_pd(sum2, sum2);           // [(a0+a2)+(a1+a3), ...]
  double sum_rho = _mm_cvtsd_f64(sum1);
  // Scalar tail for n % 4 remaining elements.
  for (; i < n; ++i) {
    double t = std::tanh((data[i] - data_offset) * half_inv_sc);
    sum_rho += t * t;
  }
  return sum_rho;
}
#endif // ROBSCALE_HAS_SLEEF && ROBSCALE_HAS_AVX2_DISPATCH

// ---------------------------------------------------------------------------
// Phase 4: TBB parallel iteration kernel
//
// Each iteration of rob_scale_compute is a reduction over n elements:
//   sum_rho = sum_i tanh((data[i]-off)*hinv)^2
// The reduction is embarrassingly parallel; the serial dep is only the
// s update between iterations.
//
// This function is called only from rob_scale_core when:
//   (a) USE_DIRECT_TBB is defined (TBB available)
//   (b) ROBSCALE_HAS_SLEEF + ROBSCALE_HAS_AVX2_DISPATCH (fused kernel available)
//   (c) n >= rob_scale_parallel_threshold (amortises TBB overhead)
//   (d) AVX2 confirmed at runtime
//
// rob_scale_compute (single-threaded) is kept unchanged for the ensemble
// bootstrap path which already parallelises over bootstrap replications.
// ---------------------------------------------------------------------------
#if (defined(ROBSCALE_HAS_SYSTEM_TBB) || defined(USE_DIRECT_TBB)) && \
    defined(ROBSCALE_HAS_SLEEF) && !defined(ROBSCALE_HAS_ACCELERATE) && \
    defined(ROBSCALE_HAS_AVX2_DISPATCH)
static double rob_scale_parallel_compute(const double* ROBSCALE_RESTRICT data,
                                          size_t n, double data_offset,
                                          double s, int maxit, double tol) {
  const double inv_n = 1.0 / static_cast<double>(n);
  const size_t grain  = robscale::qnsn::RuntimeConfig::get().grain_size;

  // Compute sum_rho via TBB parallel_reduce for a given scale.
  auto rho_sum = [&](double sc) -> double {
    const double hisc = 0.5 * robscale::INV_RHO_SCALE_CONST / sc;
    return tbb::parallel_reduce(
      tbb::blocked_range<size_t>(0, n, grain),
      0.0,
      [data, data_offset, hisc]
      (const tbb::blocked_range<size_t>& r, double acc) -> double {
        return acc + rob_scale_fused_sum_avx2(
            data + r.begin(), static_cast<int>(r.size()),
            data_offset, hisc);
      },
      std::plus<double>{}
    );
  };

  // Aitken Δ² accelerated iteration (see rob_scale_compute for rationale).
  int k = 0;
  while (k < maxit) {
    const double s0 = s;
    const double v0 = std::sqrt(2.0 * rho_sum(s0) * inv_n);
    ++k;
    const double s1 = s0 * v0;
    if (std::abs(v0 - 1.0) <= tol) { s = s1; break; }
    if (k >= maxit) { s = s1; break; }

    const double v1 = std::sqrt(2.0 * rho_sum(s1) * inv_n);
    ++k;
    const double s2 = s1 * v1;
    if (std::abs(v1 - 1.0) <= tol) { s = s2; break; }

    const double d1 = s1 - s0;
    const double d2 = s2 - s1;
    const double denom = d2 - d1;
    if (d1 * d2 > 0.0 && std::abs(d2) < std::abs(d1) &&
        std::abs(denom) > 1e-30 * s0 && k < maxit) {
      const double candidate = s2 - d2 * d2 / denom;
      if (candidate > 0.0 && std::abs(candidate - s2) < std::abs(s2 - s0)) {
        s = candidate;
        continue;
      }
    }
    s = s2;
  }
  return s;
}
#elif defined(ROBSCALE_HAS_OMP_PARALLEL) && \
      defined(ROBSCALE_HAS_SLEEF) && !defined(ROBSCALE_HAS_ACCELERATE) && \
      defined(ROBSCALE_HAS_AVX2_DISPATCH)
// OpenMP fallback: same grain-chunked structure, parallel_for + reduction.
static double rob_scale_parallel_compute(const double* ROBSCALE_RESTRICT data,
                                          size_t n, double data_offset,
                                          double s, int maxit, double tol) {
  const double inv_n   = 1.0 / static_cast<double>(n);
  const size_t grain   = robscale::qnsn::RuntimeConfig::get().grain_size;
  const size_t nchunks = (n + grain - 1) / grain;

  // Aitken Δ² accelerated iteration (OpenMP reduction, inlined for pragma compat).
  int k = 0;
  while (k < maxit) {
    const double s0 = s;
    double sr0 = 0.0;
    { const double hisc = 0.5 * robscale::INV_RHO_SCALE_CONST / s0;
#pragma omp parallel for reduction(+:sr0) schedule(static)
      for (size_t c = 0; c < nchunks; ++c) {
        size_t begin = c * grain;
        size_t end   = std::min(begin + grain, n);
        sr0 += rob_scale_fused_sum_avx2(data + begin, (int)(end - begin),
                                         data_offset, hisc);
      }
    }
    const double v0 = std::sqrt(2.0 * sr0 * inv_n);
    ++k;
    const double s1 = s0 * v0;
    if (std::abs(v0 - 1.0) <= tol) { s = s1; break; }
    if (k >= maxit) { s = s1; break; }

    double sr1 = 0.0;
    { const double hisc = 0.5 * robscale::INV_RHO_SCALE_CONST / s1;
#pragma omp parallel for reduction(+:sr1) schedule(static)
      for (size_t c = 0; c < nchunks; ++c) {
        size_t begin = c * grain;
        size_t end   = std::min(begin + grain, n);
        sr1 += rob_scale_fused_sum_avx2(data + begin, (int)(end - begin),
                                         data_offset, hisc);
      }
    }
    const double v1 = std::sqrt(2.0 * sr1 * inv_n);
    ++k;
    const double s2 = s1 * v1;
    if (std::abs(v1 - 1.0) <= tol) { s = s2; break; }

    const double d1 = s1 - s0;
    const double d2 = s2 - s1;
    const double denom = d2 - d1;
    if (d1 * d2 > 0.0 && std::abs(d2) < std::abs(d1) &&
        std::abs(denom) > 1e-30 * s0 && k < maxit) {
      const double candidate = s2 - d2 * d2 / denom;
      if (candidate > 0.0 && std::abs(candidate - s2) < std::abs(s2 - s0)) {
        s = candidate;
        continue;
      }
    }
    s = s2;
  }
  return s;
}
#endif // parallel backends

/**
 * Portably optimized robScale kernel.
 * Non-static: also called by estimators_internal.h for the ensemble.
 */
double rob_scale_compute(const double* ROBSCALE_RESTRICT data,
                         size_t n, double data_offset, double s,
                         int maxit, double tol,
                         double* ROBSCALE_RESTRICT tmp) {
  const double inv_n = 1.0 / (double)n;

  // Select the fused AVX2 path once, outside the hot loop.
  // Condition: SLEEF+AVX2 compiled in, n large enough for vectorisation,
  // and AVX2 confirmed at runtime.
#if defined(ROBSCALE_HAS_SLEEF) && !defined(ROBSCALE_HAS_ACCELERATE) && \
    defined(ROBSCALE_HAS_AVX2_DISPATCH)
  const bool use_fused = (n >= 8) &&
    (robscale::qnsn::RuntimeConfig::get().hw.simd_level >=
     robscale::qnsn::SIMDLevel::AVX2);
#else
  const bool use_fused = false;
#endif

  // Compute sum_i tanh((data[i]-offset)*half_inv_sc)^2 for a given scale.
  // Dispatches to fused AVX2 kernel or 3-pass fallback.
  auto rho_sum = [&](double sc) -> double {
    const double hisc = 0.5 * robscale::INV_RHO_SCALE_CONST / sc;
    if (use_fused) {
#if defined(ROBSCALE_HAS_SLEEF) && !defined(ROBSCALE_HAS_ACCELERATE) && \
    defined(ROBSCALE_HAS_AVX2_DISPATCH)
      return rob_scale_fused_sum_avx2(data, (int)n, data_offset, hisc);
#endif
    }
    // 3-pass fallback
    for (size_t i = 0; i < n; ++i) tmp[i] = (data[i] - data_offset) * hisc;
    robscale::bulk_tanh(tmp, (int)n);
    double sr = 0.0;
    for (size_t i = 0; i < n; ++i) sr += tmp[i] * tmp[i];
    return sr;
  };

  // Aitken Δ² (Steffensen) accelerated iteration.
  // Collects triplets (s0, s1, s2) via 2 standard steps then extrapolates:
  //   s_acc = s2 - (s2-s1)^2 / (s2 - 2*s1 + s0)
  // s_acc seeds the next round; convergence always exits via |v-1| <= tol.
  //
  // Acceptance guards:
  //   (a) Monotone + contracting: d1*d2 > 0 && |d2| < |d1| — Steffensen
  //       is unreliable for oscillating or diverging sequences (small n).
  //   (b) k < maxit: never exit the loop with an unverified Aitken value.
  int k = 0;
  while (k < maxit) {
    const double s0 = s;
    const double v0 = std::sqrt(2.0 * rho_sum(s0) * inv_n);
    ++k;
    const double s1 = s0 * v0;
    if (std::abs(v0 - 1.0) <= tol) { s = s1; break; }
    if (k >= maxit) { s = s1; break; }

    const double v1 = std::sqrt(2.0 * rho_sum(s1) * inv_n);
    ++k;
    const double s2 = s1 * v1;
    if (std::abs(v1 - 1.0) <= tol) { s = s2; break; }

    const double d1 = s1 - s0;
    const double d2 = s2 - s1;
    const double denom = d2 - d1;  // = s2 - 2*s1 + s0
    if (d1 * d2 > 0.0 && std::abs(d2) < std::abs(d1) &&
        std::abs(denom) > 1e-30 * s0 && k < maxit) {
      const double candidate = s2 - d2 * d2 / denom;
      if (candidate > 0.0 && std::abs(candidate - s2) < std::abs(s2 - s0)) {
        s = candidate;
        continue;
      }
    }
    s = s2;
  }
  return s;
}

/**
 * Shared core: median, MAD, fallback, iteration.
 * Called by both the small-n and large-n entry points.
 */
static double rob_scale_core(const double* xp, size_t n,
                             double* w, double* dev,
                             bool has_loc, double loc_val,
                             double implbound, int maxit,
                             double tol, int fallback) {
  // For n <= ROBSCALE_SORT_MEDIAN_THRESHOLD (64), median_net is always the
  // right path.  Calling it directly avoids: RuntimeConfig::get() (TLS),
  // the pdq_robscale_threshold comparison, and two function frames — ~3 ns
  // per call, ~6 ns total (median + MAD), on a ~50–200 ns small-n call.
  const bool is_small = (n <= ROBSCALE_SORT_MEDIAN_THRESHOLD);

  double t, s_init;
  if (has_loc) {
    t = loc_val;
    for (size_t i = 0; i < n; ++i) dev[i] = std::abs(xp[i] - t);
    s_init = robscale::MAD_CONSISTENCY *
        (is_small ? robscale::median_net(dev, n)
                  : robscale::adaptive_robscale_median_select(dev, n));
  } else {
    std::memcpy(w, xp, n * sizeof(double));
    t = is_small ? robscale::median_net(w, n)
                 : robscale::adaptive_robscale_median_select(w, n);
    // After median selection, w[] is permuted but holds the same multiset as
    // xp[].  Reading deviations from warm w[] instead of cold xp[] avoids an
    // extra scan of the original data; MAD is permutation-invariant.
    for (size_t i = 0; i < n; ++i) dev[i] = std::abs(w[i] - t);
    s_init = robscale::MAD_CONSISTENCY *
        (is_small ? robscale::median_net(dev, n)
                  : robscale::adaptive_robscale_median_select(dev, n));
  }

  int minobs = has_loc ? 3 : 4;
  if (ROBSCALE_UNLIKELY(n < (size_t)minobs)) {
    if (s_init <= implbound) {
      if (ROBSCALE_UNLIKELY(fallback == 1)) return R_NaReal;
      if (has_loc) {
        std::memcpy(w, xp, n * sizeof(double));
        double med_orig = robscale::adaptive_robscale_median_select(w, n);
        double mad_orig = robscale::adaptive_mad_select(xp, (int)n, med_orig, dev);
        return (mad_orig <= implbound)
          ? robscale::adm_core(xp, (int)n, med_orig, robscale::ADM_CONSISTENCY)
          : mad_orig;
      } else {
        return robscale::adm_core(xp, (int)n, t, robscale::ADM_CONSISTENCY);
      }
    }
    return s_init;
  }

  if (ROBSCALE_UNLIKELY(s_init <= implbound && fallback == 1)) return R_NaReal;
  if (ROBSCALE_UNLIKELY(s_init == 0.0)) {
    return robscale::adm_core(xp, (int)n, t, robscale::ADM_CONSISTENCY);
  }

  // Phase 4: dispatch to parallel kernel for large n.
  // Only from this call site (the user-facing path); ensemble_one_replicate
  // calls rob_scale_compute directly to avoid nesting inside its own TBB loop.
#if (defined(ROBSCALE_HAS_SYSTEM_TBB) || defined(USE_DIRECT_TBB) || \
     defined(ROBSCALE_HAS_OMP_PARALLEL)) && \
    defined(ROBSCALE_HAS_SLEEF) && !defined(ROBSCALE_HAS_ACCELERATE) && \
    defined(ROBSCALE_HAS_AVX2_DISPATCH)
  {
    const auto& cfg = robscale::qnsn::RuntimeConfig::get();
    if (n >= cfg.rob_scale_parallel_threshold &&
        cfg.hw.simd_level >= robscale::qnsn::SIMDLevel::AVX2)
      return rob_scale_parallel_compute(xp, n, t, s_init, maxit, tol);
  }
#endif

  return rob_scale_compute(xp, n, t, s_init, maxit, tol, w);
}

/**
 * Small-n entry point (n <= 64): minimal stack frame (~1KB).
 * Noinline ensures the compiler gives this its own frame, so the large
 * buf_stack in the main function doesn't penalise small-n calls.
 */
ROBSCALE_NOINLINE
static double rob_scale_impl_small(const double* xp, size_t n,
                                   bool has_loc, double loc_val,
                                   double implbound, int maxit,
                                   double tol, int fallback) {
  double arena[ROBSCALE_MICRO_BUFFER_SIZE]; // 128 doubles = 1KB
  return rob_scale_core(xp, n, arena, arena + n,
                        has_loc, loc_val, implbound, maxit, tol, fallback);
}

// [[Rcpp::export]]
double rob_scale_impl(Rcpp::NumericVector x, bool has_loc, double loc_val,
                      double implbound, int maxit, double tol, int fallback) {
  size_t n = (size_t)x.size();
  if (ROBSCALE_UNLIKELY(n == 0)) return 0.0;

  const double* xp = x.begin();

  // Small-n: dispatch to noinline helper with minimal stack frame
  if (n <= 64)
    return rob_scale_impl_small(xp, n, has_loc, loc_val,
                                implbound, maxit, tol, fallback);

  // Large-n: stack or heap arena
  constexpr size_t SCALE_STACK_SIZE = 2048;
  double buf_stack[SCALE_STACK_SIZE * 2];
  std::unique_ptr<double[]> heap;
  double* arena;

  if (ROBSCALE_LIKELY(n <= SCALE_STACK_SIZE)) {
    arena = buf_stack;
  } else {
    heap.reset(new double[n * 2]);
    arena = heap.get();
  }

  return rob_scale_core(xp, n, arena, arena + n,
                        has_loc, loc_val, implbound, maxit, tol, fallback);
}
