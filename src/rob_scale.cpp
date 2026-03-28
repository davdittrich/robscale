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

#include "validate_finite.h"

// ---------------------------------------------------------------------------
// Fused single-pass NR step kernels
//
// nr_scale_step_scalar — portable fallback, std::tanh per element
// nr_scale_step_avx2   — AVX2+FMA, 4-wide, scalar tail for n%4
//
// Each reads data[] once, computes u_i, tanh(u_i), and accumulates
// sum_tanh2 and sum_u_tanh_sech2 via FMA — no scratch buffer.
// ---------------------------------------------------------------------------

// Scalar single-pass NR step.
static void nr_scale_step_scalar(const double* ROBSCALE_RESTRICT data,
                                  int n, double data_offset, double hisc,
                                  double* out_sum_tanh2,
                                  double* out_sum_u_tanh_sech2) {
  double sum_tanh2 = 0.0, sum_u_tanh_sech2 = 0.0;
  for (int i = 0; i < n; ++i) {
    const double ui  = (data[i] - data_offset) * hisc;
    const double tv  = std::tanh(ui);
    const double tv2 = tv * tv;
    sum_tanh2        += tv2;
    sum_u_tanh_sech2 += ui * tv * (1.0 - tv2);
  }
  *out_sum_tanh2        = sum_tanh2;
  *out_sum_u_tanh_sech2 = sum_u_tanh_sech2;
}

#if defined(ROBSCALE_HAS_AVX2_TANH) && !defined(ROBSCALE_HAS_ACCELERATE)
// AVX2 fused single-pass NR step.
// Loads data[] once, computes u = (d - off) * h inline, calls
// ROBSCALE_TANH4_AVX2 (libmvec preferred, SLEEF fallback) in 4-wide batches,
// and accumulates via _mm256_fmadd_pd. Scalar tail handles n%4 elements.
ROBSCALE_TARGET_AVX2
static void nr_scale_step_avx2(const double* ROBSCALE_RESTRICT data,
                                int n, double data_offset, double hisc,
                                double* out_sum_tanh2,
                                double* out_sum_u_tanh_sech2) {
  const __m256d off4 = _mm256_set1_pd(data_offset);
  const __m256d h4   = _mm256_set1_pd(hisc);
  const __m256d one4 = _mm256_set1_pd(1.0);
  __m256d acc_t2  = _mm256_setzero_pd();
  __m256d acc_uts = _mm256_setzero_pd();

  int i = 0;
  for (; i + 4 <= n; i += 4) {
    __m256d d    = _mm256_loadu_pd(data + i);
    __m256d u    = _mm256_mul_pd(_mm256_sub_pd(d, off4), h4);
    __m256d tv   = ROBSCALE_TANH4_AVX2(u);
    __m256d tv2  = _mm256_mul_pd(tv, tv);
    acc_t2       = _mm256_add_pd(acc_t2, tv2);
    __m256d s2   = _mm256_sub_pd(one4, tv2);   // sech²(u) = 1 - tanh²(u)
    __m256d u_tv = _mm256_mul_pd(u, tv);
    acc_uts      = _mm256_fmadd_pd(u_tv, s2, acc_uts);
  }

  // Horizontal reduction: acc_t2
  __m128d lo = _mm256_castpd256_pd128(acc_t2);
  __m128d hi = _mm256_extractf128_pd(acc_t2, 1);
  __m128d s128 = _mm_add_pd(lo, hi);
  s128 = _mm_hadd_pd(s128, s128);
  double sum_tanh2 = _mm_cvtsd_f64(s128);

  // Horizontal reduction: acc_uts
  lo = _mm256_castpd256_pd128(acc_uts);
  hi = _mm256_extractf128_pd(acc_uts, 1);
  s128 = _mm_add_pd(lo, hi);
  s128 = _mm_hadd_pd(s128, s128);
  double sum_u_tanh_sech2 = _mm_cvtsd_f64(s128);

  // Scalar tail for n%4 remaining elements
  for (; i < n; ++i) {
    const double ui  = (data[i] - data_offset) * hisc;
    const double tv  = std::tanh(ui);
    const double tv2 = tv * tv;
    sum_tanh2        += tv2;
    sum_u_tanh_sech2 += ui * tv * (1.0 - tv2);
  }

  *out_sum_tanh2        = sum_tanh2;
  *out_sum_u_tanh_sech2 = sum_u_tanh_sech2;
}
#endif  // ROBSCALE_HAS_AVX2_TANH

#if (defined(ROBSCALE_HAS_SYSTEM_TBB) || defined(USE_DIRECT_TBB)) && \
    defined(ROBSCALE_HAS_AVX2_TANH) && !defined(ROBSCALE_HAS_ACCELERATE)
static double rob_scale_parallel_compute(const double* ROBSCALE_RESTRICT data,
                                          size_t n, double data_offset,
                                          double s, int maxit, double tol) {
  const double inv_n    = 1.0 / static_cast<double>(n);
  const double hisc_pre = 0.5 * robscale::INV_RHO_SCALE_CONST;
  const size_t grain    = robscale::qnsn::RuntimeConfig::get().grain_size;

  struct Accum { double s1 = 0.0; double s2 = 0.0; };

  for (int iter = 0; iter < maxit; ) {
    const double hisc = hisc_pre / s;

    Accum a = tbb::parallel_reduce(
      tbb::blocked_range<size_t>(0, n, grain),
      Accum{},
      [data, data_offset, hisc](const tbb::blocked_range<size_t>& r, Accum acc) -> Accum {
        double st2 = 0.0, suts = 0.0;
        nr_scale_step_avx2(data + r.begin(), static_cast<int>(r.size()),
                            data_offset, hisc, &st2, &suts);
        acc.s1 += st2;
        acc.s2 += suts;
        return acc;
      },
      [](Accum x, Accum y) -> Accum { return {x.s1 + y.s1, x.s2 + y.s2}; }
    );

    const double numer   = a.s1 * inv_n - 0.5;
    const double denom   = 2.0 * inv_n * a.s2;
    if (std::abs(denom) <= 1e-14 * s) { s *= std::sqrt(2.0 * a.s1 * inv_n); ++iter; continue; }
    const double delta_s = s * numer / denom;
    if (s + delta_s <= 0.0)           { s /= 2.0;                           ++iter; continue; }
    s += delta_s;
    ++iter;
    if (std::abs(delta_s) / s <= tol) break;
  }
  return s;
}
#elif defined(ROBSCALE_HAS_OMP_PARALLEL) && \
      defined(ROBSCALE_HAS_AVX2_TANH) && !defined(ROBSCALE_HAS_ACCELERATE)
// OpenMP fallback: NR two-sum reduction; each grain calls nr_scale_step_avx2.
static double rob_scale_parallel_compute(const double* ROBSCALE_RESTRICT data,
                                          size_t n, double data_offset,
                                          double s, int maxit, double tol) {
  const double inv_n    = 1.0 / static_cast<double>(n);
  const double hisc_pre = 0.5 * robscale::INV_RHO_SCALE_CONST;
  const size_t grain    = robscale::qnsn::RuntimeConfig::get().grain_size;
  const size_t nchunks  = (n + grain - 1) / grain;

  for (int iter = 0; iter < maxit; ) {
    const double hisc = hisc_pre / s;
    double sum_tanh2 = 0.0, sum_u_tanh_sech2 = 0.0;

#pragma omp parallel for reduction(+:sum_tanh2,sum_u_tanh_sech2) schedule(static)
    for (size_t c = 0; c < nchunks; ++c) {
      const size_t begin = c * grain;
      const size_t end   = std::min(begin + grain, n);
      double st2 = 0.0, suts = 0.0;
      nr_scale_step_avx2(data + begin, static_cast<int>(end - begin),
                          data_offset, hisc, &st2, &suts);
      sum_tanh2        += st2;
      sum_u_tanh_sech2 += suts;
    }

    const double numer   = sum_tanh2 * inv_n - 0.5;
    const double denom   = 2.0 * inv_n * sum_u_tanh_sech2;
    if (std::abs(denom) <= 1e-14 * s) { s *= std::sqrt(2.0 * sum_tanh2 * inv_n); ++iter; continue; }
    const double delta_s = s * numer / denom;
    if (s + delta_s <= 0.0)           { s /= 2.0;                                ++iter; continue; }
    s += delta_s;
    ++iter;
    if (std::abs(delta_s) / s <= tol) break;
  }
  return s;
}
#endif // parallel backends

/**
 * Newton-Raphson scale iteration for M-scale. Scratch-free.
 *
 * Each iteration dispatches nr_scale_step_avx2 (AVX2, libmvec/SLEEF) or
 * nr_scale_step_scalar — a single-pass fused kernel that reads data[] once,
 * computes u_i = (data[i] - data_offset) * hisc inline, and returns
 * (sum_tanh2, sum_u_tanh_sech2) without any scratch buffer.
 *
 * data:    w[] = |x_i - t| when called from rob_scale_core (data_offset=0, OPT-F).
 *          sorted_x when called via rob_scale_compute from ensemble path.
 * use_avx2: pre-cached flag (OPT-L2); avoids RuntimeConfig::get() per iter.
 */
static double nr_scale_compute(const double* ROBSCALE_RESTRICT data,
                                size_t n, double data_offset,
                                double s, int maxit, double tol,
                                bool use_avx2) {
  const double inv_n    = 1.0 / static_cast<double>(n);
  const double hisc_pre = 0.5 * robscale::INV_RHO_SCALE_CONST;
  int iter = 0;

  for (; iter < maxit; ) {
    const double hisc = hisc_pre / s;
    double sum_tanh2 = 0.0, sum_u_tanh_sech2 = 0.0;

    // Fused single-pass: compute u_i and tanh(u_i) inline, accumulate both sums.
#if defined(ROBSCALE_HAS_AVX2_TANH) && !defined(ROBSCALE_HAS_ACCELERATE)
    if (use_avx2 && static_cast<int>(n) >= 4)
      nr_scale_step_avx2(data, static_cast<int>(n), data_offset, hisc,
                         &sum_tanh2, &sum_u_tanh_sech2);
    else
#endif
      nr_scale_step_scalar(data, static_cast<int>(n), data_offset, hisc,
                           &sum_tanh2, &sum_u_tanh_sech2);

    const double numer = sum_tanh2 * inv_n - 0.5;
    const double denom = 2.0 * inv_n * sum_u_tanh_sech2;

    // GP guard: denominator degenerate → multiplicative fallback
    if (std::abs(denom) <= 1e-14 * s) {
      s *= std::sqrt(2.0 * sum_tanh2 * inv_n);
      ++iter; continue;
    }

    const double delta_s = s * numer / denom;

    // Neg-s guard: proposed update non-positive → halve
    if (s + delta_s <= 0.0) { s /= 2.0; ++iter; continue; }

    s += delta_s;
    ++iter;
    if (std::abs(delta_s) / s <= tol) break;
  }
  return s;
}

/**
 * robScale kernel — delegates to nr_scale_compute (NR iteration).
 * Non-static: called from estimators_internal.h (ensemble path) and
 * rob_scale_sorted. No scratch allocation — fused kernel is scratch-free.
 */
ROBSCALE_HIDDEN
double rob_scale_compute(const double* ROBSCALE_RESTRICT data,
                         size_t n, double data_offset, double s,
                         int maxit, double tol,
                         bool use_avx2) {
  return nr_scale_compute(data, n, data_offset, s, maxit, tol, use_avx2);
}

/**
 * Shared core: median, MAD, fallback, iteration.
 * Called by both the small-n and large-n entry points.
 */
static double rob_scale_core(const double* ROBSCALE_RESTRICT xp, size_t n,
                             double* ROBSCALE_RESTRICT w,
                             bool has_loc, double loc_val,
                             double implbound, int maxit,
                             double tol, int fallback) {
  // For n <= ROBSCALE_SORT_NETWORK_THRESHOLD (16), median_net is always the
  // right path.  Calling it directly avoids: RuntimeConfig::get() (TLS),
  // the pdq_robscale_threshold comparison, and two function frames — ~3 ns
  // per call, ~6 ns total (median + MAD), on a ~50–200 ns small-n call.
  const bool is_small = (n <= ROBSCALE_SORT_NETWORK_THRESHOLD);

  // OPT-E: cache RuntimeConfig once — used for both the parallel threshold
  // check (below) and the use_avx2 flag passed to rob_scale_compute.
  // Eliminates a second TLS lookup on every rob_scale_core call.
#if defined(ROBSCALE_HAS_AVX2_TANH) && !defined(ROBSCALE_HAS_ACCELERATE)
  const auto& cfg = robscale::qnsn::RuntimeConfig::get();
  const bool avx2 = (cfg.hw.simd_level >= robscale::qnsn::SIMDLevel::AVX2);
#else
  const bool avx2 = false;
#endif

  double t, s_init;
  if (has_loc) {
    t = loc_val;
    // OPT-3: SIMD-annotated abs-diff kernel (no scalar loop).
    // OPT-F: deviations stored in w[] (offset=0 dispatch below).
    robscale::bulk_abs_diff(w, xp, (int)n, t);
    s_init = robscale::MAD_CONSISTENCY *
        (is_small ? robscale::median_net(w, n)
                  : robscale::adaptive_robscale_median_select(w, n));
  } else {
    std::memcpy(w, xp, n * sizeof(double));
    t = is_small ? robscale::median_net(w, n)
                 : robscale::adaptive_robscale_median_select(w, n);
    // OPT-3: compute deviations from xp (source), overwrite w (destination).
    // w was permuted by median_select; xp preserves original order.
    // OPT-F: pre-subtract so rob_scale_compute gets offset=0, saving one
    // VSUB per element per NR iteration.  MAD is permutation-invariant.
    robscale::bulk_abs_diff(w, xp, (int)n, t);
    s_init = robscale::MAD_CONSISTENCY *
        (is_small ? robscale::median_net(w, n)
                  : robscale::adaptive_robscale_median_select(w, n));
  }

  int minobs = has_loc ? 3 : 4;
  if (ROBSCALE_UNLIKELY(n < (size_t)minobs)) {
    if (s_init <= implbound) {
      if (ROBSCALE_UNLIKELY(fallback == 1)) return R_NaReal;
      if (has_loc) {
        std::memcpy(w, xp, n * sizeof(double));
        double med_orig = robscale::adaptive_robscale_median_select(w, n);
        double dev_local[4]; // n < minobs (3 when has_loc), always fits
        double mad_orig = robscale::adaptive_mad_select(xp, (int)n, med_orig, dev_local);
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
    defined(ROBSCALE_HAS_AVX2_TANH) && !defined(ROBSCALE_HAS_ACCELERATE)
  // OPT-6: UNLIKELY — parallel threshold only crossed for large n (rare path).
  // Uses xp with data_offset=t (signed deviations); serial path below uses w
  // with data_offset=0 (OPT-F abs-deviations). Both equivalent: the NR sums
  // (tanh² and u·tanh·sech²) are even functions of u.
  if (ROBSCALE_UNLIKELY(n >= cfg.rob_scale_parallel_threshold && avx2))
    return rob_scale_parallel_compute(xp, n, t, s_init, maxit, tol);
#endif

  // OPT-F: pass w[] (abs-deviations) with data_offset=0 — eliminates one VSUB
  // per element per iteration. dev[] no longer needed as NR scratch (fused kernel).
  return nr_scale_compute(w, n, 0.0, s_init, maxit, tol, avx2);
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
  // OPT-5: 32-byte alignment for AVX2 load/store paths.
  // Arena holds w[] only — dev eliminated (allocated locally on rare fallback path)
  alignas(32) double arena[ROBSCALE_MICRO_BUFFER_SIZE]; // 128 doubles = 1KB
  return rob_scale_core(xp, n, arena,
                        has_loc, loc_val, implbound, maxit, tol, fallback);
}

// [[Rcpp::export(rng = false)]]
double rob_scale_impl(Rcpp::NumericVector x, bool has_loc, double loc_val,
                      double implbound, int maxit, double tol, int fallback) {
  size_t n = (size_t)x.size();
  if (ROBSCALE_UNLIKELY(n == 0)) return 0.0;

  const double* xp = x.begin();
  validate_finite(xp, (int)n);

  // Small-n: dispatch to noinline helper with minimal stack frame
  if (n <= 64)
    return rob_scale_impl_small(xp, n, has_loc, loc_val,
                                implbound, maxit, tol, fallback);

  // Large-n: stack or heap arena
  constexpr size_t SCALE_STACK_SIZE = 2048;
  // OPT-5: 32-byte alignment for AVX2 load/store paths.
  // Arena holds w[] only — dev eliminated (allocated locally on rare fallback path)
  alignas(32) double buf_stack[SCALE_STACK_SIZE];
  std::unique_ptr<double[]> heap;
  double* arena;

  if (ROBSCALE_LIKELY(n <= SCALE_STACK_SIZE)) {
    arena = buf_stack;
  } else {
    heap.reset(new double[n]);
    arena = heap.get();
  }

  return rob_scale_core(xp, n, arena,
                        has_loc, loc_val, implbound, maxit, tol, fallback);
}

// Thin wrapper for gate benchmarks — calls production rob_scale_impl directly.
// [[Rcpp::export(rng = false)]]
double C_rob_scale_fast(Rcpp::NumericVector x) {
  // Thin wrapper → production rob_scale_impl (has_loc=false, default params).
  // Always reflects the current production path; after WU-RS1 calls 8-wide kernel.
  return rob_scale_impl(x, false, 0.0, 1e-4, 80, 1.4901161e-8, 0);
}
