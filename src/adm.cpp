#include <Rcpp.h>
#include "robust_core.h"

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
  double* buf = (n <= STACK_BUF_SIZE) ? buf_stack : new double[n];
  std::memcpy(buf, x.begin(), n * sizeof(double));
  double med = median_select(buf, n);
  if (n > STACK_BUF_SIZE) delete[] buf;
  return adm_core(x.begin(), n, med, constant);
}
