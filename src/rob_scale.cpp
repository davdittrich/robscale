#include "robscale_config.h"
#include "robust_core.h"
#include "pdq_select.h"
#include <Rcpp.h>
#include <memory>

// ---------------------------------------------------------------------------
// Phase 2: Fused per-iteration kernel
//
// The 3-pass loop (scale → tanh → sum-of-squares) writes and reads tmp[]
// 4 times per element per iteration.  For large n this dominates: at n=10^6
// with 32 iterations that is 1 GB of tmp traffic.
//
// The fused kernel collapses all three passes into one, reading data[] once
// and never touching tmp[].  Uses Sleef_tanhd4_u10avx2 + 4-wide FMA
// accumulation; FP results are within 1e-12 of the 3-pass path (golden
// test tolerance) even though the accumulation order differs slightly.
// ---------------------------------------------------------------------------

#if defined(ROBSCALE_HAS_SLEEF) && !defined(ROBSCALE_HAS_ACCELERATE) && \
    defined(ROBSCALE_HAS_AVX2_DISPATCH)
/**
 * Fused sum: compute sum_i tanh((data[i]-off)*hinv)^2 in one pass.
 *
 * Uses Sleef_tanhd4_u10avx2 + FMA accumulation into 4 independent SIMD
 * lanes, avoiding the serial scalar-extraction dependency chain that would
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
    __m256d t = Sleef_tanhd4_u10avx2(
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

  for (int k = 0; k < maxit; ++k) {
    const double half_inv_sc = 0.5 * robscale::INV_RHO_SCALE_CONST / s;
    double sum_rho;

    if (use_fused) {
#if defined(ROBSCALE_HAS_SLEEF) && !defined(ROBSCALE_HAS_ACCELERATE) && \
    defined(ROBSCALE_HAS_AVX2_DISPATCH)
      sum_rho = rob_scale_fused_sum_avx2(data, (int)n, data_offset, half_inv_sc);
#endif
    } else {
      // 3-pass fallback: unchanged from pre-Phase-2, bit-identical on all
      // non-AVX2 targets.
      for (size_t i = 0; i < n; ++i)
        tmp[i] = (data[i] - data_offset) * half_inv_sc;
      robscale::bulk_tanh(tmp, (int)n);
      sum_rho = 0.0;
      for (size_t i = 0; i < n; ++i)
        sum_rho += tmp[i] * tmp[i];
    }

    double v = std::sqrt(2.0 * sum_rho * inv_n);
    s *= v;
    if (std::abs(v - 1.0) <= tol) break;
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
  double t, s_init;
  if (has_loc) {
    t = loc_val;
    for (size_t i = 0; i < n; ++i) dev[i] = std::abs(xp[i] - t);
    s_init = robscale::MAD_CONSISTENCY * robscale::adaptive_robscale_median_select(dev, n);
  } else {
    std::memcpy(w, xp, n * sizeof(double));
    t = robscale::adaptive_robscale_median_select(w, n);
    s_init = robscale::adaptive_mad_select(xp, (int)n, t, dev);
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
