#include "robscale_config.h"
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

static ROBSCALE_NOINLINE
double iqr_impl_large(const double* ROBSCALE_RESTRICT xp, int n, double constant) {
  constexpr int STACK_SIZE = 2048;
  double buf_stack[STACK_SIZE];
  std::unique_ptr<double[]> heap;
  double* buf = (n <= STACK_SIZE)
    ? buf_stack
    : (heap.reset(new double[n]), heap.get());

  std::memcpy(buf, xp, n * sizeof(double));

  // Q1/Q3 target indices (Type 7)
  double h1 = (n - 1.0) * 0.25;
  int lo1 = static_cast<int>(h1);
  double frac1 = h1 - lo1;

  double h3 = (n - 1.0) * 0.75;
  int lo3 = static_cast<int>(h3);
  double frac3 = h3 - lo3;

  // Q1: select on full array
  miniselect::pdqselect(buf, buf + lo1, buf + n);
  double q1 = robscale::interp_q7(buf, n, lo1, frac1);

  // Q3: select on buf[lo1+1 .. n-1] only (everything <= Q1 is irrelevant)
  int start = lo1 + 1;
  miniselect::pdqselect(buf + start, buf + lo3, buf + n);
  double q3 = robscale::interp_q7(buf, n, lo3, frac3);

  return (q3 - q1) * constant;
}

// [[Rcpp::export]]
double iqr_impl(Rcpp::NumericVector x, double constant) {
  int n = x.size();
  if (n < 2) return 0.0;

  if (n <= IQR_INLINE_LIMIT) {
    // Micro path: buf_micro is the only stack allocation in this frame.
    // IQR_INLINE_LIMIT=256 (2KB) vs old frame (buf_micro[128]+buf_stack[4096]=33KB).
    double buf_micro[IQR_INLINE_LIMIT];
    const double* xp = x.begin();
    std::memcpy(buf_micro, xp, n * sizeof(double));

    double h1 = (n - 1.0) * 0.25;
    int lo1 = static_cast<int>(h1);
    double frac1 = h1 - lo1;

    double h3 = (n - 1.0) * 0.75;
    int lo3 = static_cast<int>(h3);
    double frac3 = h3 - lo3;

    miniselect::pdqselect(buf_micro, buf_micro + lo1, buf_micro + n);
    double q1 = robscale::interp_q7(buf_micro, n, lo1, frac1);

    int start = lo1 + 1;
    miniselect::pdqselect(buf_micro + start, buf_micro + lo3, buf_micro + n);
    double q3 = robscale::interp_q7(buf_micro, n, lo3, frac3);

    return (q3 - q1) * constant;
  }
  return iqr_impl_large(x.begin(), n, constant);
}

