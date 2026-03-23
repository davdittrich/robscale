#include "robscale_config.h"
#include "robust_core.h"
#include <Rcpp.h>
#include <limits>
#include <memory>
#if defined(ROBSCALE_HAS_SYSTEM_TBB)
#  include <oneapi/tbb/parallel_reduce.h>
#  include <oneapi/tbb/blocked_range.h>
#elif defined(USE_DIRECT_TBB)
#  include <tbb/parallel_reduce.h>
#  include <tbb/blocked_range.h>
#endif

// ---------------------------------------------------------------------------
// OPT-L1: Fused AVX2 NR step
//
// Computes sum_psi  = Σ tanh(u_i)  and
//          sum_dpsi = Σ (1 - tanh²(u_i))   [= Σ sech²(u_i)]
// where u_i = (xp[i] - t) * half_inv_s, in a single pass over xp[].
//
// Two accumulators:
//   acc_psi  tracks Σ tanh(u)   — signed, in (-n, n)
//   acc_p2   tracks Σ tanh²(u)  — non-negative; sum_dpsi = range_n - sum_p2
//
// fmadd for p²: acc_p2 = fmadd(p, p, acc_p2) = p*p + acc_p2.
// One FMA rounding step preserves non-negativity; avoids the two-rounding
// error of mul(p,p) + add.
//
// Precision: 4-wide accumulation differs from scalar by ≤ n*eps.
// At n=1000: max diff ≈ 2.2e-13, well within 2*sqrt(eps) = 2.98e-8 tolerance.
//
// @param range_n  number of elements to process
// @param out_psi   [out] Σ tanh(u_i)
// @param out_dpsi  [out] Σ (1 - tanh²(u_i)) = range_n - Σ tanh²(u_i)
// ---------------------------------------------------------------------------
#if defined(ROBSCALE_HAS_SLEEF) && !defined(ROBSCALE_HAS_ACCELERATE) && \
    defined(ROBSCALE_HAS_AVX2_DISPATCH)
ROBSCALE_TARGET_AVX2
static void rob_loc_nr_step_avx2(const double* ROBSCALE_RESTRICT xp, int range_n,
                                  double t, double half_inv_s,
                                  double* out_psi, double* out_dpsi) {
  const __m256d t4    = _mm256_set1_pd(t);
  const __m256d hinv4 = _mm256_set1_pd(half_inv_s);
  __m256d acc_psi = _mm256_setzero_pd();
  __m256d acc_p2  = _mm256_setzero_pd();

  int i = 0;
  for (; i + 4 <= range_n; i += 4) {
    __m256d d = _mm256_loadu_pd(xp + i);
    __m256d u = _mm256_mul_pd(_mm256_sub_pd(d, t4), hinv4);
    __m256d p = ROBSCALE_TANH4_AVX2(u);
    acc_psi = _mm256_add_pd(acc_psi, p);
    acc_p2  = _mm256_fmadd_pd(p, p, acc_p2);   // acc_p2 += p*p (one rounding step)
  }

  // Horizontal sum of 4-wide SIMD accumulators
  __m128d lo_psi = _mm256_castpd256_pd128(acc_psi);
  __m128d hi_psi = _mm256_extractf128_pd(acc_psi, 1);
  __m128d s2_psi = _mm_add_pd(lo_psi, hi_psi);
  double sum_psi = _mm_cvtsd_f64(_mm_hadd_pd(s2_psi, s2_psi));

  __m128d lo_p2  = _mm256_castpd256_pd128(acc_p2);
  __m128d hi_p2  = _mm256_extractf128_pd(acc_p2, 1);
  __m128d s2_p2  = _mm_add_pd(lo_p2, hi_p2);
  double sum_p2  = _mm_cvtsd_f64(_mm_hadd_pd(s2_p2, s2_p2));

  // Scalar tail for range_n % 4 remaining elements
  for (; i < range_n; ++i) {
    double p  = std::tanh((xp[i] - t) * half_inv_s);
    sum_psi  += p;
    sum_p2   += p * p;
  }

  *out_psi  = sum_psi;
  *out_dpsi = (double)range_n - sum_p2;   // Σ(1-p²) = range_n - Σp²
}
#endif // ROBSCALE_HAS_SLEEF && ROBSCALE_HAS_AVX2_DISPATCH

// ---------------------------------------------------------------------------
// OPT-L3: TBB parallel NR reduction
//
// For n >= rob_scale_parallel_threshold, the per-NR-iteration sum over n
// elements is embarrassingly parallel.  Each TBB chunk calls rob_loc_nr_step_avx2
// on its sub-range; the NRAccum struct carries (psi, dpsi) partial sums.
//
// The NR loop itself is serial across iterations (each update to t depends on
// the previous t).  Only the inner sum per iteration is parallelized.
// ---------------------------------------------------------------------------
#if (defined(ROBSCALE_HAS_SYSTEM_TBB) || defined(USE_DIRECT_TBB)) && \
    defined(ROBSCALE_HAS_SLEEF) && !defined(ROBSCALE_HAS_ACCELERATE) && \
    defined(ROBSCALE_HAS_AVX2_DISPATCH)

struct NRAccum {
  double psi{0.0}, dpsi{0.0};
  void operator+=(const NRAccum& o) { psi += o.psi; dpsi += o.dpsi; }
};

static double rob_loc_parallel_compute(const double* ROBSCALE_RESTRICT xp,
                                        size_t n, double t, double s,
                                        int maxit, double tol) {
  const double half_inv_s = 0.5 / s;
  const size_t grain = robscale::qnsn::RuntimeConfig::get().grain_size;

  for (int k = 0; k < maxit; ++k) {
    // Parallel reduction: each grain calls the fused AVX2 step on its sub-range.
    // NRAccum combines (partial_psi, partial_dpsi) from each chunk.
    NRAccum acc = tbb::parallel_reduce(
      tbb::blocked_range<size_t>(0, n, grain),
      NRAccum{},
      [xp, t, half_inv_s](const tbb::blocked_range<size_t>& r, NRAccum init) -> NRAccum {
        double psi, dpsi;
        rob_loc_nr_step_avx2(xp + r.begin(), static_cast<int>(r.size()),
                              t, half_inv_s, &psi, &dpsi);
        init.psi  += psi;
        init.dpsi += dpsi;
        return init;
      },
      [](NRAccum a, const NRAccum& b) -> NRAccum { a += b; return a; }
    );

    // sum_dpsi guard: Σ sech²(u_i) → 0 when all |u_i| >> 1.
    if (ROBSCALE_UNLIKELY(acc.dpsi < std::numeric_limits<double>::min())) break;

    double v = 2.0 * s * acc.psi / acc.dpsi;
    t += v;
    if (std::abs(v) <= tol) break;
  }
  return t;
}
#endif // TBB + SLEEF + AVX2

/**
 * Portably optimized robLoc NR kernel (serial).
 *
 * OPT-L1: dispatch to fused AVX2 single-pass kernel when available.
 *   Scalar 3-pass fallback: scale → tanh → accumulate.
 * OPT-L2: RuntimeConfig::get() hoisted once before the NR loop.
 *   bulk_tanh_dispatched() skips the repeated dispatch check per iteration.
 *
 * sum_dpsi guard: if Σ sech²(u_i) underflows (all |u_i| >> 1, degenerate
 * scale), the NR step blows up. Break early with the current t.
 */
static ROBSCALE_INLINE double rob_loc_compute(const double* ROBSCALE_RESTRICT xp,
                                              size_t n, double t, double s,
                                              int maxit, double tol,
                                              double* ROBSCALE_RESTRICT tmp) {
  const double half_inv_s = 0.5 / s;

  // OPT-L1+L2: hoist SIMD dispatch check once before the NR loop.
  // use_fused:     dispatch to rob_loc_nr_step_avx2 (fused single-pass kernel).
  // use_avx2_tanh: same flag — if AVX2+n>=8, also use it for the 3-pass fallback
  //               (only reached when use_fused=false, i.e., n<8 or non-AVX2 build).
#if defined(ROBSCALE_HAS_SLEEF) && !defined(ROBSCALE_HAS_ACCELERATE) && \
    defined(ROBSCALE_HAS_AVX2_DISPATCH)
  const bool use_fused = (n >= 8) &&
    (robscale::qnsn::RuntimeConfig::get().hw.simd_level >=
     robscale::qnsn::SIMDLevel::AVX2);
  const bool use_avx2_tanh = use_fused;
#else
  const bool use_fused     = false;
  const bool use_avx2_tanh = false;
#endif

  for (int k = 0; k < maxit; ++k) {
    double sum_psi, sum_dpsi;

#if defined(ROBSCALE_HAS_SLEEF) && !defined(ROBSCALE_HAS_ACCELERATE) && \
    defined(ROBSCALE_HAS_AVX2_DISPATCH)
    if (use_fused) {
      rob_loc_nr_step_avx2(xp, (int)n, t, half_inv_s, &sum_psi, &sum_dpsi);
    } else
#endif
    {
      // 3-pass scalar fallback
      for (size_t i = 0; i < n; ++i)
        tmp[i] = (xp[i] - t) * half_inv_s;
      robscale::bulk_tanh_dispatched(tmp, (int)n, use_avx2_tanh);
      sum_psi = 0.0; sum_dpsi = 0.0;
      for (size_t i = 0; i < n; ++i) {
        double p = tmp[i];
        sum_psi  += p;
        sum_dpsi += 1.0 - p * p;
      }
    }

    // sum_dpsi guard: Σ sech²(u_i) → 0 when all |u_i| >> 1 (degenerate scale).
    if (ROBSCALE_UNLIKELY(sum_dpsi < std::numeric_limits<double>::min())) break;

    double v = 2.0 * s * sum_psi / sum_dpsi;
    t += v;
    if (std::abs(v) <= tol) break;
  }
  return t;
}

/**
 * Shared core: median, MAD, iteration.
 * Called by both the small-n and large-n entry points.
 */
static double rob_loc_core(const double* xp, size_t n,
                           double* buf, double* dev,
                           bool has_scale, double scale_val,
                           int maxit, double tol) {
  std::memcpy(buf, xp, n * sizeof(double));
  double med = robscale::median_select(buf, n);

  int minobs = has_scale ? 3 : 4;
  if (ROBSCALE_UNLIKELY(n < (size_t)minobs)) return med;

  // OPT-L4: pass buf[] (warm in L1/L2 after memcpy+median_select) instead of
  // xp[] (potentially evicted from cache).  MAD is permutation-invariant, so
  // mad_select(buf, ...) == mad_select(xp, ...) by construction.
  double s = has_scale ? scale_val : robscale::mad_select(buf, (int)n, med, dev);
  if (ROBSCALE_UNLIKELY(s == 0.0)) return med;

  // OPT-L3: dispatch to parallel NR for large n.
  // Condition: TBB+SLEEF+AVX2 compiled in, n >= rob_scale_parallel_threshold
  // (reusing robScale's threshold — same per-element work), AVX2 at runtime.
#if (defined(ROBSCALE_HAS_SYSTEM_TBB) || defined(USE_DIRECT_TBB)) && \
    defined(ROBSCALE_HAS_SLEEF) && !defined(ROBSCALE_HAS_ACCELERATE) && \
    defined(ROBSCALE_HAS_AVX2_DISPATCH)
  {
    const auto& cfg = robscale::qnsn::RuntimeConfig::get();
    if (n >= cfg.rob_scale_parallel_threshold &&
        cfg.hw.simd_level >= robscale::qnsn::SIMDLevel::AVX2)
      return rob_loc_parallel_compute(xp, n, med, s, maxit, tol);
  }
#endif

  return rob_loc_compute(xp, n, med, s, maxit, tol, buf);
}

/**
 * Small-n entry point (n <= 64): minimal stack frame (~1KB).
 */
ROBSCALE_NOINLINE
static double rob_loc_impl_small(const double* xp, size_t n,
                                 bool has_scale, double scale_val,
                                 int maxit, double tol) {
  double arena[ROBSCALE_MICRO_BUFFER_SIZE]; // 128 doubles = 1KB
  return rob_loc_core(xp, n, arena, arena + n,
                      has_scale, scale_val, maxit, tol);
}

// [[Rcpp::export]]
double rob_loc_impl(Rcpp::NumericVector x, bool has_scale, double scale_val,
                    int maxit, double tol) {
  size_t n = (size_t)x.size();
  if (ROBSCALE_UNLIKELY(n == 0)) return 0.0;

  const double* xp = x.begin();

  // Small-n: dispatch to noinline helper with minimal stack frame
  if (n <= 64)
    return rob_loc_impl_small(xp, n, has_scale, scale_val, maxit, tol);

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

  return rob_loc_core(xp, n, arena, arena + n,
                      has_scale, scale_val, maxit, tol);
}

// ---------------------------------------------------------------------------
// Diagnostic exports — Phase 3+4 TDD gates
// ---------------------------------------------------------------------------

// OPT-RL0: H2H baseline — verbatim snapshot of robLoc dispatch before RL1..RL3.
//
// Frozen behaviors (for H2H comparison against optimised rob_loc_impl):
//   - use_fused threshold: n >= 8  (RL1 changes to n >= 4)
//   - scalar fallback:    3-pass via buf[]  (RL2 replaces with single-pass)
//   - NR data source:     xp (cold original)  (RL3 changes to warm buf)
//
// Does NOT include TBB dispatch (gate sizes n <= 1000 are all below the ≥ 4096
// parallel threshold, so omitting it produces an identical result for all gate inputs).
//
// Remove after WU-RL3 gate passes; run Rcpp::compileAttributes() afterwards.
// [[Rcpp::export]]
double rob_loc_fast_orig(Rcpp::NumericVector x, bool has_scale, double scale_val,
                         int maxit, double tol) {
  size_t n = (size_t)x.size();
  if (ROBSCALE_UNLIKELY(n == 0)) return 0.0;
  const double* xp = x.begin();

  // Arena setup — mirrors rob_loc_impl (micro ≤ 64, stack ≤ 2048, heap)
  constexpr size_t MICRO_SIZE = 64;
  constexpr size_t STACK_SIZE = 2048;
  double micro_buf[ROBSCALE_MICRO_BUFFER_SIZE]; // 128 doubles (n*2 at n ≤ 64)
  double stack_buf[STACK_SIZE * 2];
  std::unique_ptr<double[]> heap_buf;
  double* buf;
  double* dev;

  if (n <= MICRO_SIZE) {
    buf = micro_buf;      dev = micro_buf  + n;
  } else if (ROBSCALE_LIKELY(n <= STACK_SIZE)) {
    buf = stack_buf;      dev = stack_buf  + n;
  } else {
    heap_buf.reset(new double[n * 2]);
    buf = heap_buf.get(); dev = heap_buf.get() + n;
  }

  // Median (OPT-L4: warm buf, not cold xp)
  std::memcpy(buf, xp, n * sizeof(double));
  double med = robscale::median_select(buf, n);

  int minobs = has_scale ? 3 : 4;
  if (ROBSCALE_UNLIKELY(n < (size_t)minobs)) return med;

  double s = has_scale ? scale_val : robscale::mad_select(buf, (int)n, med, dev);
  if (ROBSCALE_UNLIKELY(s == 0.0)) return med;

  // NR iteration — FROZEN pre-RL1..RL3 state
  const double half_inv_s = 0.5 / s;
  double t = med;

#if defined(ROBSCALE_HAS_SLEEF) && !defined(ROBSCALE_HAS_ACCELERATE) && \
    defined(ROBSCALE_HAS_AVX2_DISPATCH)
  const bool use_fused = (n >= 8) &&  // FROZEN: n >= 8  (RL1 changes to n >= 4)
    (robscale::qnsn::RuntimeConfig::get().hw.simd_level >=
     robscale::qnsn::SIMDLevel::AVX2);
#else
  const bool use_fused = false;
#endif

  for (int k = 0; k < maxit; ++k) {
    double sum_psi, sum_dpsi;

#if defined(ROBSCALE_HAS_SLEEF) && !defined(ROBSCALE_HAS_ACCELERATE) && \
    defined(ROBSCALE_HAS_AVX2_DISPATCH)
    if (use_fused) {
      rob_loc_nr_step_avx2(xp, (int)n, t, half_inv_s, &sum_psi, &sum_dpsi); // FROZEN: xp
    } else
#endif
    {
      // 3-pass scalar — FROZEN (RL2 replaces with single-pass; RL3 uses buf not xp)
      for (size_t i = 0; i < n; ++i)
        buf[i] = (xp[i] - t) * half_inv_s;              // FROZEN: reads xp
      robscale::bulk_tanh_dispatched(buf, (int)n, false);
      sum_psi = 0.0; sum_dpsi = 0.0;
      for (size_t i = 0; i < n; ++i) {
        double p = buf[i];
        sum_psi  += p;
        sum_dpsi += 1.0 - p * p;
      }
    }

    if (ROBSCALE_UNLIKELY(sum_dpsi < std::numeric_limits<double>::min())) break;
    double v = 2.0 * s * sum_psi / sum_dpsi;
    t += v;
    if (std::abs(v) <= tol) break;
  }
  return t;
}

// Phase 3 gate (test 0.3): always uses scalar 3-pass NR path.
// robLoc() dispatches to rob_loc_nr_step_avx2 when AVX2 is present.
// Assert: |rob_loc_scalar_impl(x) - robLoc(x)| < 2*sqrt(eps).
// [[Rcpp::export]]
double rob_loc_scalar_impl(Rcpp::NumericVector x) {
  static const double TOL   = std::sqrt(std::numeric_limits<double>::epsilon());
  static const int    MAXIT = 80;

  const size_t n = (size_t)x.size();
  if (n == 0) return 0.0;

  const double* xp = x.begin();

  std::unique_ptr<double[]> arena(new double[n * 2]);
  double* buf = arena.get();
  double* dev = arena.get() + n;

  std::memcpy(buf, xp, n * sizeof(double));
  double med = robscale::median_select(buf, n);
  if (n < 4) return med;

  double s = robscale::mad_select(buf, (int)n, med, dev);
  if (s == 0.0) return med;

  // 3-pass scalar NR — no SIMD dispatch
  const double half_inv_s = 0.5 / s;
  double t = med;

  for (int k = 0; k < MAXIT; ++k) {
    for (size_t i = 0; i < n; ++i)
      buf[i] = (xp[i] - t) * half_inv_s;
    for (size_t i = 0; i < n; ++i)
      buf[i] = std::tanh(buf[i]);
    double sum_psi = 0.0, sum_dpsi = 0.0;
    for (size_t i = 0; i < n; ++i) {
      double p = buf[i];
      sum_psi  += p;
      sum_dpsi += 1.0 - p * p;
    }
    if (ROBSCALE_UNLIKELY(sum_dpsi < std::numeric_limits<double>::min())) break;
    double v = 2.0 * s * sum_psi / sum_dpsi;
    t += v;
    if (std::abs(v) <= TOL) break;
  }
  return t;
}

// Phase 4 gate (test 0.5): returns TRUE if TBB parallel path is compiled and
// AVX2 is confirmed at runtime.
// [[Rcpp::export]]
bool rob_loc_has_parallel() {
#if (defined(ROBSCALE_HAS_SYSTEM_TBB) || defined(USE_DIRECT_TBB)) && \
    defined(ROBSCALE_HAS_SLEEF) && !defined(ROBSCALE_HAS_ACCELERATE) && \
    defined(ROBSCALE_HAS_AVX2_DISPATCH)
  return robscale::qnsn::RuntimeConfig::get().hw.simd_level >=
         robscale::qnsn::SIMDLevel::AVX2;
#else
  return false;
#endif
}

// Phase 4 gate (test 0.5): forces the serial rob_loc_compute path, bypassing
// the parallel dispatch in rob_loc_core.  Used to cross-check parallel result.
// [[Rcpp::export]]
double rob_loc_serial_impl(Rcpp::NumericVector x) {
  static const double TOL   = std::sqrt(std::numeric_limits<double>::epsilon());
  static const int    MAXIT = 80;

  const size_t n = (size_t)x.size();
  if (n == 0) return 0.0;

  const double* xp = x.begin();

  std::unique_ptr<double[]> arena(new double[n * 2]);
  double* buf = arena.get();
  double* dev = arena.get() + n;

  std::memcpy(buf, xp, n * sizeof(double));
  double med = robscale::median_select(buf, n);
  if (n < 4) return med;

  double s = robscale::mad_select(buf, (int)n, med, dev);
  if (s == 0.0) return med;

  // Call serial rob_loc_compute directly (no TBB dispatch)
  return rob_loc_compute(xp, n, med, s, MAXIT, TOL, buf);
}
