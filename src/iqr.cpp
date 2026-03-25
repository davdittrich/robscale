#include "robscale_config.h"
#include "estimators_internal.h"
#include "pdq_select.h"
#include <Rcpp.h>
#include <cstring>
#include <memory>

// OPT-I1: Extract large-n path as ROBSCALE_NOINLINE so buf_stack[2048] is
// never allocated in the entry frame when n <= IQR_INLINE_LIMIT.
// IQR_INLINE_LIMIT=256 pushes the NOINLINE boundary above the noisy ~6µs
// region (n=129) where call overhead is disproportionate to computation.
// OPT-I2: STACK_SIZE reduced from 4096 to 2048 (IQR needs one working copy).
static constexpr int IQR_INLINE_LIMIT = 256;

// OPT-I3: symmetric Q1 selection for frac1>0 — pdqselect to lo1+1 then
// max-scan [0..lo1] O(0.25n) instead of min-scan [lo1+1..n-1] O(0.75n).
static ROBSCALE_NOINLINE
double iqr_impl_large(const double* ROBSCALE_RESTRICT xp, int n, double constant) {
  constexpr int STACK_SIZE = ROBSCALE_SN_STACK_THRESHOLD;
  double buf_stack[STACK_SIZE];
  std::unique_ptr<double[]> heap;
  double* buf = (n <= STACK_SIZE)
    ? buf_stack
    : (heap.reset(new double[n]), heap.get());

  std::memcpy(buf, xp, n * sizeof(double));

  // R4: delegate to shared helper (eliminates copy-paste with micro path below).
  return robscale::internal::iqr_select_and_interp(buf, n) * constant;
}

// OPT-I6: n<=16 sort-once-then-index helper.
// NOINLINE keeps the code for the n>16 micro-path compact (avoids I-cache
// pressure from inlining the sorting-network code into the pdqselect path).
static ROBSCALE_NOINLINE
double iqr_impl_small(const double* ROBSCALE_RESTRICT xp, int n, double constant) {
  double buf[16];
  std::memcpy(buf, xp, static_cast<size_t>(n) * sizeof(double));
  double h1 = (n - 1.0) * 0.25;
  int lo1 = static_cast<int>(h1);
  double frac1 = h1 - lo1;
  double h3 = (n - 1.0) * 0.75;
  int lo3 = static_cast<int>(h3);
  double frac3 = h3 - lo3;
  robscale::small_sort(buf, static_cast<size_t>(n));
  double q1 = buf[lo1];
  if (frac1 > 0.0) q1 += frac1 * (buf[lo1 + 1] - q1);
  double q3 = buf[lo3];
  if (frac3 > 0.0 && lo3 + 1 < n) q3 += frac3 * (buf[lo3 + 1] - q3);
  return (q3 - q1) * constant;
}

// [[Rcpp::export]]
double iqr_impl(Rcpp::NumericVector x, double constant) {
  int n = x.size();
  if (n < 2) return 0.0;

  const double* xp = x.begin();

  // OPT-I6: LIKELY branch first so the n>16 pdqselect/micro-path code is laid
  // out inline (immediately after function entry). The n<=16 sort path is a
  // NOINLINE out-of-line call — no I-cache pressure on the hot n>16 path.
  if (ROBSCALE_LIKELY(n > 16)) {
    if (n > IQR_INLINE_LIMIT) return iqr_impl_large(xp, n, constant);

    // Micro path (17 <= n <= IQR_INLINE_LIMIT):
    // R4: delegate to shared helper (eliminates copy-paste with iqr_impl_large).
    double buf_micro[IQR_INLINE_LIMIT];
    std::memcpy(buf_micro, xp, n * sizeof(double));
    return robscale::internal::iqr_select_and_interp(buf_micro, n) * constant;
  }

  return iqr_impl_small(xp, n, constant);
}


