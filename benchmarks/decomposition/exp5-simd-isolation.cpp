// [[Rcpp::plugins(openmp)]]

// Platform-specific SIMD detection
#if defined(__APPLE__) && defined(__MACH__)
  #include <vecLib/vForce.h>
  #define ROBSCALE_HAS_ACCELERATE 1
#endif

#include <Rcpp.h>
#include <cmath>
#include <cstring>
#include <algorithm>
#include <memory>

constexpr double MAD_CONSISTENCY = 1.4826;
constexpr int STACK_BUF = 512;

// ---- Sorting networks (identical to production robscale) ----
#define SWAP_IF(a, b) \
  do { if ((a) > (b)) { double t_ = (a); (a) = (b); (b) = t_; } } while(0)

inline void sort_net_2(double* x) { SWAP_IF(x[0], x[1]); }
inline void sort_net_3(double* x) {
  SWAP_IF(x[0], x[1]); SWAP_IF(x[0], x[2]); SWAP_IF(x[1], x[2]);
}
inline void sort_net_4(double* x) {
  SWAP_IF(x[0], x[1]); SWAP_IF(x[2], x[3]);
  SWAP_IF(x[0], x[2]); SWAP_IF(x[1], x[3]); SWAP_IF(x[1], x[2]);
}
inline void sort_net_5(double* x) {
  SWAP_IF(x[0], x[1]); SWAP_IF(x[2], x[3]);
  SWAP_IF(x[0], x[2]); SWAP_IF(x[1], x[3]); SWAP_IF(x[1], x[2]);
  SWAP_IF(x[0], x[4]); SWAP_IF(x[2], x[4]);
  SWAP_IF(x[1], x[2]); SWAP_IF(x[3], x[4]);
}
inline void sort_net_6(double* x) {
  SWAP_IF(x[0], x[5]); SWAP_IF(x[1], x[3]); SWAP_IF(x[2], x[4]);
  SWAP_IF(x[1], x[2]); SWAP_IF(x[3], x[4]);
  SWAP_IF(x[0], x[3]); SWAP_IF(x[2], x[5]);
  SWAP_IF(x[0], x[1]); SWAP_IF(x[2], x[3]); SWAP_IF(x[4], x[5]);
  SWAP_IF(x[1], x[2]); SWAP_IF(x[3], x[4]);
}
inline void sort_net_7(double* x) {
  SWAP_IF(x[0], x[1]); SWAP_IF(x[2], x[3]); SWAP_IF(x[4], x[5]);
  SWAP_IF(x[0], x[2]); SWAP_IF(x[1], x[3]); SWAP_IF(x[4], x[6]);
  SWAP_IF(x[5], x[6]);
  SWAP_IF(x[0], x[4]); SWAP_IF(x[1], x[5]); SWAP_IF(x[2], x[6]);
  SWAP_IF(x[2], x[4]); SWAP_IF(x[3], x[5]);
  SWAP_IF(x[1], x[2]); SWAP_IF(x[3], x[4]); SWAP_IF(x[5], x[6]);
  SWAP_IF(x[2], x[3]);
}
inline void sort_net_8(double* x) {
  SWAP_IF(x[0], x[1]); SWAP_IF(x[2], x[3]); SWAP_IF(x[4], x[5]);
  SWAP_IF(x[6], x[7]);
  SWAP_IF(x[0], x[2]); SWAP_IF(x[1], x[3]); SWAP_IF(x[4], x[6]);
  SWAP_IF(x[5], x[7]);
  SWAP_IF(x[1], x[2]); SWAP_IF(x[5], x[6]);
  SWAP_IF(x[0], x[4]); SWAP_IF(x[1], x[5]); SWAP_IF(x[2], x[6]);
  SWAP_IF(x[3], x[7]);
  SWAP_IF(x[2], x[4]); SWAP_IF(x[3], x[5]);
  SWAP_IF(x[1], x[2]); SWAP_IF(x[3], x[4]); SWAP_IF(x[5], x[6]);
}

inline void small_sort(double* x, int n) {
  switch (n) {
    case 0: case 1: return;
    case 2: sort_net_2(x); return;
    case 3: sort_net_3(x); return;
    case 4: sort_net_4(x); return;
    case 5: sort_net_5(x); return;
    case 6: sort_net_6(x); return;
    case 7: sort_net_7(x); return;
    case 8: sort_net_8(x); return;
    default: std::sort(x, x + n); return;
  }
}

// ---- Median / MAD ----
inline double median_sorted(const double* x, int n) {
  if (n & 1) return x[n >> 1];
  return (x[(n >> 1) - 1] + x[n >> 1]) * 0.5;
}

inline double median_select(double* buf, int n) {
  if (n <= 8) { small_sort(buf, n); return median_sorted(buf, n); }
  int k = (n - 1) / 2;
  std::nth_element(buf, buf + k, buf + n);
  if (n & 1) return buf[k];
  double hi = buf[k + 1];
  for (int i = k + 2; i < n; ++i) if (buf[i] < hi) hi = buf[i];
  return (buf[k] + hi) * 0.5;
}

inline double mad_select(const double* x, int n, double med, double* dev) {
  for (int i = 0; i < n; ++i) dev[i] = std::abs(x[i] - med);
  return MAD_CONSISTENCY * median_select(dev, n);
}

// ==== TANH VARIANTS (the only dimension that differs) ====

// --- Scalar: plain loop, no vectorization hints ---
inline void bulk_tanh_scalar(double* inout, int n) {
  for (int i = 0; i < n; ++i)
    inout[i] = std::tanh(inout[i]);
}

// --- SIMD: platform-default vectorization ---
#if defined(ROBSCALE_HAS_ACCELERATE)
inline void bulk_tanh_simd(double* inout, int n) {
  vvtanh(inout, inout, &n);
}
#else
  #ifdef _OPENMP
    #pragma omp declare simd notinbranch
  #endif
  inline double vec_tanh(double x) { return std::tanh(x); }
  inline void bulk_tanh_simd(double* inout, int n) {
    #ifdef _OPENMP
      #pragma omp simd
    #endif
    for (int i = 0; i < n; ++i)
      inout[i] = vec_tanh(inout[i]);
  }
#endif

// ==== TEMPLATED robLoc (parameterized on tanh function) ====
// All instantiations share identical machine code for sorting networks,
// median, MAD, NR structure, and stack arena. The ONLY difference is
// the tanh function pointer.
template<void (*TanhFn)(double*, int)>
static double robLoc_impl(const double* xp, int n, int maxit, double tol) {
  double stk[STACK_BUF * 2];
  std::unique_ptr<double[]> heap;
  double* arena;
  if (n * 2 <= STACK_BUF * 2) { arena = stk; }
  else { heap.reset(new double[n * 2]); arena = heap.get(); }
  double* buf = arena;
  double* dev = arena + n;

  std::memcpy(buf, xp, n * sizeof(double));
  double med = median_select(buf, n);
  if (n < 4) return med;

  double s = mad_select(xp, n, med, dev);
  if (s == 0.0) return med;

  double t = med, half_inv_s = 0.5 / s;
  for (int k = 0; k < maxit; ++k) {
    for (int i = 0; i < n; ++i)
      buf[i] = (xp[i] - t) * half_inv_s;
    TanhFn(buf, n);               // *** ONLY DIFFERENCE ***
    double sp = 0.0, sd = 0.0;
    for (int i = 0; i < n; ++i) {
      double p = buf[i];
      sp += p;
      sd += 1.0 - p * p;
    }
    double v = 2.0 * s * sp / sd;
    t += v;
    if (std::abs(v) <= tol) break;
  }
  return t;
}

// [[Rcpp::export]]
double robLoc_scalar(Rcpp::NumericVector x, int maxit = 200, double tol = 1e-12) {
  return robLoc_impl<bulk_tanh_scalar>(x.begin(), x.size(), maxit, tol);
}

// [[Rcpp::export]]
double robLoc_simd(Rcpp::NumericVector x, int maxit = 200, double tol = 1e-12) {
  return robLoc_impl<bulk_tanh_simd>(x.begin(), x.size(), maxit, tol);
}

// ==== ADM (no tanh involved, single implementation) ====
// [[Rcpp::export]]
double adm_sourcecpp(Rcpp::NumericVector x) {
  int n = x.size();
  double stk[STACK_BUF];
  std::unique_ptr<double[]> heap;
  double* buf;
  if (n <= STACK_BUF) { buf = stk; }
  else { heap.reset(new double[n]); buf = heap.get(); }
  std::memcpy(buf, x.begin(), n * sizeof(double));
  double med = median_select(buf, n);
  double sum = 0.0;
  const double* xp = x.begin();
  for (int i = 0; i < n; ++i)
    sum += std::abs(xp[i] - med);
  return 1.2533141373155001 * sum / n;
}

// ==== SLEEF: hand-optimized AVX2 vector tanh (x86_64) ====
#include <sleef.h>
#include <immintrin.h>

inline void bulk_tanh_sleef(double* inout, int n) {
  int i = 0;
  for (; i + 4 <= n; i += 4) {
    __m256d v = _mm256_loadu_pd(inout + i);
    v = Sleef_tanhd4_u10avx2(v);
    _mm256_storeu_pd(inout + i, v);
  }
  for (; i < n; i++)
    inout[i] = std::tanh(inout[i]);
}

// [[Rcpp::export]]
double robLoc_sleef(Rcpp::NumericVector x, int maxit = 200, double tol = 1e-12) {
  return robLoc_impl<bulk_tanh_sleef>(x.begin(), x.size(), maxit, tol);
}

