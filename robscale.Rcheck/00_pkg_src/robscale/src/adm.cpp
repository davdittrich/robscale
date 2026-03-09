#include "robust_core.h"
#include <Rcpp.h>

// Called when center is explicitly provided by user
// [[Rcpp::export]]
double adm_impl(Rcpp::NumericVector x, double center, double constant) {
  return adm_core(x.begin(), x.size(), center, constant);
}

// Called when center is not provided — compute median via O(n) selection
// [[Rcpp::export]]
double adm_impl_auto(Rcpp::NumericVector x, double constant) {
  int n = x.size();
  double buf_stack[STACK_BUF_SIZE];
  std::unique_ptr<double[]> heap;
  double* buf;
  if (n <= STACK_BUF_SIZE) {
    buf = buf_stack;
  } else {
    heap.reset(new double[n]);
    buf = heap.get();
  }
  std::memcpy(buf, x.begin(), n * sizeof(double));
  double med = median_select(buf, n);
  return adm_core(x.begin(), n, med, constant);
}
