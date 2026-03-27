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
#if defined(ROBSCALE_HAS_AVX2_TANH) && !defined(ROBSCALE_HAS_ACCELERATE)
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
#endif // ROBSCALE_HAS_AVX2_TANH

// ---------------------------------------------------------------------------
// OPT-L3: TBB parallel quadratic NR reduction
//
// Same quadratic NR as rob_loc_compute: sum_dpsi = Σ sech²((x_i−t)/(2s)) is
// the observed Hessian, recomputed each iteration from the current t.
// T'(t*) = 0 → quadratic local convergence (not linear IRLS).
//
// For n >= rob_scale_parallel_threshold, the per-NR-iteration inner sum over
// n elements is embarrassingly parallel.  Each TBB chunk calls
// rob_loc_nr_step_avx2 on its sub-range; NRAccum carries (psi, dpsi) partials.
//
// The outer NR loop is serial across iterations (each update to t depends on
// the previous t).  Only the inner sum per iteration is parallelized.
// ---------------------------------------------------------------------------
#if (defined(ROBSCALE_HAS_SYSTEM_TBB) || defined(USE_DIRECT_TBB)) && \
    defined(ROBSCALE_HAS_AVX2_TANH) && !defined(ROBSCALE_HAS_ACCELERATE)

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
 * CONVERGENCE: Quadratic Newton-Raphson, not linear IRLS.
 *   sum_dpsi = Σ sech²((x_i−t)/(2s)) is recomputed each iteration from the
 *   current t (observed Hessian / observed Fisher information).  This gives
 *   T'(t*) = 0 and quadratic local convergence, unlike IRLS which uses a
 *   fixed expected-information denominator and converges linearly with
 *   |T'(t*)| < 1.  Consequence: typically 2–4 iterations on real data;
 *   maxit=80 is a safety limit only.  Aitken/Steffensen acceleration
 *   provides no benefit and can increase evaluation count.
 *
 * OPT-L1: dispatch to fused AVX2 single-pass kernel when available.
 * OPT-L2: RuntimeConfig::get() hoisted once before the NR loop.
 * OPT-RL2: scalar fallback is now single-pass (no tmp[] writes; accumulators
 *   stay in registers). Eliminates the 3-pass scale->tanh->accumulate pattern.
 *
 * sum_dpsi guard: if Σ sech²(u_i) underflows (all |u_i| >> 1, degenerate
 * scale), the NR step blows up. Break early with the current t.
 */
static ROBSCALE_INLINE double rob_loc_compute(const double* ROBSCALE_RESTRICT xp,
                                              size_t n, double t, double s,
                                              int maxit, double tol) {
  const double half_inv_s = 0.5 / s;

  // OPT-L1+L2: hoist SIMD dispatch check once before the NR loop.
  // use_fused: dispatch to rob_loc_nr_step_avx2 (fused single-pass kernel).
  //   Scalar fallback (n<4 or non-AVX2): single-pass std::tanh loop.
#if defined(ROBSCALE_HAS_AVX2_TANH) && !defined(ROBSCALE_HAS_ACCELERATE)
  const bool use_fused = (n >= 4) &&  // OPT-RL1: lowered from n>=8 to n>=4
    (robscale::qnsn::RuntimeConfig::get().hw.simd_level >=
     robscale::qnsn::SIMDLevel::AVX2);
#endif

  for (int k = 0; k < maxit; ++k) {
    double sum_psi, sum_dpsi;

#if defined(ROBSCALE_HAS_AVX2_TANH) && !defined(ROBSCALE_HAS_ACCELERATE)
    if (use_fused) {
      rob_loc_nr_step_avx2(xp, (int)n, t, half_inv_s, &sum_psi, &sum_dpsi);
    } else
#endif
    {
      // OPT-RL2: single-pass fused scalar — no tmp[] writes, accumulators in registers.
      // OPT-RL5: omp simd hint for non-AVX2 builds (no-op when AVX2+SLEEF are active).
      sum_psi = 0.0; sum_dpsi = 0.0;
#if defined(_OPENMP) || defined(ROBSCALE_HAS_OMP_SIMD)
#pragma omp simd reduction(+:sum_psi,sum_dpsi)
#endif
      for (size_t i = 0; i < n; ++i) {
        double p = std::tanh((xp[i] - t) * half_inv_s);
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
static double rob_loc_core(const double* ROBSCALE_RESTRICT xp, size_t n,
                           double* ROBSCALE_RESTRICT buf,
                           double* ROBSCALE_RESTRICT dev,
                           bool has_scale, double scale_val,
                           int maxit, double tol) {
  std::memcpy(buf, xp, n * sizeof(double));

  // OPT-RL5: hoist small-n dispatch check — avoids repeated branch inside
  // median_select.  For n <= ROBSCALE_SORT_NETWORK_THRESHOLD (16), call
  // median_net directly (mirrors rob_scale_core pattern).
  const bool is_small = (n <= ROBSCALE_SORT_NETWORK_THRESHOLD);
  double med = is_small ? robscale::median_net(buf, n)
                        : robscale::median_select(buf, n);

  int minobs = has_scale ? 3 : 4;
  if (ROBSCALE_UNLIKELY(n < (size_t)minobs)) return med;

  // OPT-L4 + OPT-RL3: pass buf[] (warm in L1/L2 after memcpy+median_select) instead
  // of xp[] (potentially evicted from cache).  Both MAD and NR sums are
  // permutation-invariant, so buf and xp produce identical results.
  double s = has_scale ? scale_val : robscale::mad_select(buf, (int)n, med, dev);
  if (ROBSCALE_UNLIKELY(s == 0.0)) return med;

  // OPT-L3: dispatch to parallel NR for large n.
  // Condition: TBB+SLEEF+AVX2 compiled in, n >= rob_scale_parallel_threshold
  // (reusing robScale's threshold — same per-element work), AVX2 at runtime.
#if (defined(ROBSCALE_HAS_SYSTEM_TBB) || defined(USE_DIRECT_TBB)) && \
    defined(ROBSCALE_HAS_AVX2_TANH) && !defined(ROBSCALE_HAS_ACCELERATE)
  {
    const auto& cfg = robscale::qnsn::RuntimeConfig::get();
    if (n >= cfg.rob_scale_parallel_threshold &&
        cfg.hw.simd_level >= robscale::qnsn::SIMDLevel::AVX2)
      return rob_loc_parallel_compute(buf, n, med, s, maxit, tol);  // OPT-RL3: warm buf
  }
#endif

  return rob_loc_compute(buf, n, med, s, maxit, tol);  // OPT-RL3: warm buf
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
    defined(ROBSCALE_HAS_AVX2_TANH) && !defined(ROBSCALE_HAS_ACCELERATE)
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

  // Call serial rob_loc_compute directly (no TBB dispatch).
  // Intentionally passes xp (original data order) — serial reference for cross-checking.
  return rob_loc_compute(xp, n, med, s, MAXIT, TOL);
}

// Aitken correctness cross-check — plain NR, no Aitken acceleration.
// Self-contained: does NOT call rob_loc_compute, so it remains an unmodified
// plain-NR reference after WU-RL-A2 applies Aitken to that function.
// Assert: |rob_loc_noaitken_impl(x) - robLoc(x)| < 2*sqrt(eps)
// after WU-RL-A2 (same fixed point, different convergence path).
// [[Rcpp::export]]
double rob_loc_noaitken_impl(Rcpp::NumericVector x) {
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

  // Single-pass scalar NR — matches the production scalar fallback in rob_loc_compute
  // (OPT-RL2 pattern: accumulators in registers, no tmp[] writes), no SIMD, no Aitken.
  const double half_inv_s = 0.5 / s;
  double t = med;

  for (int k = 0; k < MAXIT; ++k) {
    double sum_psi = 0.0, sum_dpsi = 0.0;
    for (size_t i = 0; i < n; ++i) {
      double p = std::tanh((xp[i] - t) * half_inv_s);
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

// Returns the number of NR evaluations for plain NR.
// Finding (WU-RL-A2): robLoc uses quadratic NR (observed Hessian sum_dpsi = Σ sech²),
// NOT linearly converging IRLS. Aitken was therefore not added to the production path.
// This diagnostic confirms ≤ 4 evaluations on typical data.
// Not exported (test 0.46 guards against re-export).
int rob_loc_noaitken_iters(Rcpp::NumericVector x) {
  static const double TOL   = std::sqrt(std::numeric_limits<double>::epsilon());
  static const int    MAXIT = 80;

  const size_t n = (size_t)x.size();
  if (n == 0) return 0;

  const double* xp = x.begin();

  std::unique_ptr<double[]> arena(new double[n * 2]);
  double* buf = arena.get();
  double* dev = arena.get() + n;

  std::memcpy(buf, xp, n * sizeof(double));
  double med = robscale::median_select(buf, n);
  if (n < 4) return 0;

  double s = robscale::mad_select(buf, (int)n, med, dev);
  if (s == 0.0) return 0;

  const double half_inv_s = 0.5 / s;
  double t = med;
  int neval = 0;

  for (int k = 0; k < MAXIT; ++k) {
    ++neval;
    double sum_psi = 0.0, sum_dpsi = 0.0;
    for (size_t i = 0; i < n; ++i) {
      double p = std::tanh((xp[i] - t) * half_inv_s);
      sum_psi  += p;
      sum_dpsi += 1.0 - p * p;
    }
    if (ROBSCALE_UNLIKELY(sum_dpsi < std::numeric_limits<double>::min())) break;
    double v = 2.0 * s * sum_psi / sum_dpsi;
    t += v;
    if (std::abs(v) <= TOL) break;
  }
  return neval;
}
