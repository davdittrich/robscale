#include "robust_core.h"
#include <Rcpp.h>
#include <memory>
#ifndef STACK_BUF_SIZE
#define STACK_BUF_SIZE 1024
#endif

// [[Rcpp::export]]
double adm_impl(Rcpp::NumericVector x, double center, double constant) {
  return robscale::adm_core(x.begin(), (int)x.size(), center, constant);
}

// [[Rcpp::export]]
double adm_impl_auto(Rcpp::NumericVector x, double constant) {
  int n = x.size();
  if (n == 0) return 0.0;
  
  double med;
  if (n <= 64) {
    double buf_micro[64];
    std::memcpy(buf_micro, x.begin(), n * sizeof(double));
    med = robscale::median_select(buf_micro, n);
  } else {
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
    med = robscale::median_select(buf, n);
  }
  return robscale::adm_core(x.begin(), n, med, constant);
}
