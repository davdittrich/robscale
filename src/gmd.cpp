#include "robscale_config.h"
#include "robust_core.h"
#include "qnsn_sort_utils.h"
#include <Rcpp.h>
#include <cstring>
#include <memory>

#include "validate_finite.h"

// GMD sorted kernel: assumes x is already sorted ascending.
// OPT-G6: ROBSCALE_RESTRICT allows the compiler to assume no aliasing.
// OPT-G5: caller precomputes scale = constant * 2 / (n*(n-1)) once.
// WU-GMD-1: delegates to gmd_weighted_sum in robust_core.h which dispatches
// to AVX2 FMA 4-wide kernel or scalar fallback.
// use_avx2: pre-hoisted flag from caller; avoids TLS read per call.
static ROBSCALE_INLINE double gmd_sorted(const double* ROBSCALE_RESTRICT x,
                                          int n, double scale, bool use_avx2) {
  return robscale::gmd_weighted_sum(x, n, scale, use_avx2);
}

// Shared sort+kernel logic.
// OPT-G5: computes scale once before the loop.
// OPT-G6: ROBSCALE_RESTRICT on both pointer parameters.
static double gmd_core(const double* ROBSCALE_RESTRICT xp, int n,
                       double* ROBSCALE_RESTRICT buf, double constant) {
  const double scale = constant * 2.0 / (static_cast<double>(n) * (n - 1));
  std::memcpy(buf, xp, n * sizeof(double));
  if (static_cast<size_t>(n) <= ROBSCALE_SORT_NET_THRESHOLD)
    robscale::small_sort(buf, static_cast<size_t>(n));
  else
    robscale::qnsn::optimized_sort(buf, buf + n);
#if defined(ROBSCALE_HAS_AVX2_TANH) && !defined(ROBSCALE_HAS_ACCELERATE)
  const bool use_avx2 = (robscale::qnsn::RuntimeConfig::get().hw.simd_level >=
                         robscale::qnsn::SIMDLevel::AVX2);
#else
  const bool use_avx2 = false;
#endif
  return gmd_sorted(buf, n, scale, use_avx2);
}

// OPT-G1: small-n path with minimal stack frame (1 KB arena only).
// ROBSCALE_NOINLINE isolates this frame so buf_stack[2048] is never
// allocated when n <= ROBSCALE_MICRO_BUFFER_SIZE.
ROBSCALE_NOINLINE
static double gmd_impl_small(const double* xp, int n, double constant) {
  double arena[ROBSCALE_MICRO_BUFFER_SIZE];  // 128 doubles = 1 KB
  return gmd_core(xp, n, arena, constant);
}

// [[Rcpp::export(rng = false)]]
double gmd_impl(Rcpp::NumericVector x, double constant) {
  int n = x.size();
  if (n < 2) return 0.0;
  validate_finite(x.begin(), n);
  if (n <= ROBSCALE_MICRO_BUFFER_SIZE)
    return gmd_impl_small(x.begin(), n, constant);
  constexpr int STACK_SIZE = ROBSCALE_SN_STACK_THRESHOLD;
  double buf_stack[STACK_SIZE];
  std::unique_ptr<double[]> heap;
  double* buf = (n <= STACK_SIZE)
    ? buf_stack
    : (heap.reset(new double[n]), heap.get());
  return gmd_core(x.begin(), n, buf, constant);
}

