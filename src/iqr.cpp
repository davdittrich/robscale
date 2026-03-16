#include "robscale_config.h"
#include "pdq_select.h"
#include <Rcpp.h>
#include <cstring>
#include <memory>

// [[Rcpp::export]]
double iqr_impl(Rcpp::NumericVector x, double constant) {
  int n = x.size();
  if (n < 2) return 0.0;

  // Arena: single copy (incremental Q3 reuses Q1's partition)
  double buf_micro[ROBSCALE_MICRO_BUFFER_SIZE];
  constexpr int STACK_SIZE = 4096;
  double buf_stack[STACK_SIZE];
  std::unique_ptr<double[]> heap;
  double* buf;

  if (n <= ROBSCALE_MICRO_BUFFER_SIZE) {
    buf = buf_micro;
  } else if (n <= STACK_SIZE) {
    buf = buf_stack;
  } else {
    heap.reset(new double[n]);
    buf = heap.get();
  }

  std::memcpy(buf, x.begin(), n * sizeof(double));

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
