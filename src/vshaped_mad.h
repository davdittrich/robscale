#ifndef ROBSCALE_VSHAPED_MAD_H
#define ROBSCALE_VSHAPED_MAD_H

#include <algorithm>

namespace robscale {

// Find the 0-indexed k-th smallest of A[0..la-1] union B[0..lb-1].
// Both arrays sorted ascending. O(log(min(la, lb))).
// Precondition: 0 <= k < la + lb.
inline double kth_of_two_sorted(
    const double* A, int la, const double* B, int lb, int k) {
  if (la > lb) return kth_of_two_sorted(B, lb, A, la, k);
  if (la == 0) return B[k];
  if (k == 0)  return std::min(A[0], B[0]);
  int ia = std::min(la - 1, k / 2);
  int ib = k - ia - 1;
  if (A[ia] < B[ib])
    return kth_of_two_sorted(A + ia + 1, la - ia - 1, B, lb, k - ia - 1);
  if (A[ia] > B[ib])
    return kth_of_two_sorted(A, la, B + ib + 1, lb - ib - 1, k - ib - 1);
  return A[ia];
}

// V-shaped deviation median: MAD of sorted s[0..n-1] with known median m.
// tmp: caller-supplied scratch of >= n doubles.
// Complexity: O(log n). Returns 0.0 for n < 2.
inline double vshaped_mad(const double* s, int n, double m, double* tmp) {
  if (n < 2) return 0.0;
  const int k = (n - 1) / 2;
  for (int i = 0; i <= k; ++i)     tmp[i] = m - s[k - i];
  for (int i = k + 1; i < n; ++i)  tmp[i] = s[i] - m;
  const double* L = tmp;
  const double* R = tmp + k + 1;
  const int la = k + 1;
  const int lb = n - k - 1;
  if (n & 1) {
    return kth_of_two_sorted(L, la, R, lb, k);
  } else {
    double lo = kth_of_two_sorted(L, la, R, lb, k);
    double hi = kth_of_two_sorted(L, la, R, lb, k + 1);
    return lo + (hi - lo) * 0.5;  // overflow-safe averaging
  }
}

} // namespace robscale

#endif // ROBSCALE_VSHAPED_MAD_H
