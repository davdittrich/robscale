#include <iostream>
#include <vector>
#include <algorithm>
#include <stdint.h>
#include <stddef.h>

#define ROBSCALE_UNLIKELY(x) (x)
#define NA_REAL 0.0

template <typename T>
inline T whimed_cpp(T* a, int32_t* iw, size_t n, int64_t target) {
  if (n == 0) return T(0);
  if (n == 1) return a[0];
  size_t l = 0, r = n - 1;
  int64_t t = target;
  while (l < r) {
    T pivot = a[l + (r - l) / 2];
    size_t i = l, j = l;
    while (j <= r) {
      if (a[j] < pivot) {
        std::swap(a[i], a[j]); std::swap(iw[i], iw[j]); i++;
      }
      j++;
    }
    int64_t wleft = 0;
    for (size_t idx = l; idx < i; ++idx) wleft += iw[idx];
    if (wleft > t) {
      r = (i > l) ? i - 1 : l;
    } else {
      size_t i_eq = i, j_eq = i;
      while (j_eq <= r) {
        if (a[j_eq] == pivot) {
          std::swap(a[i_eq], a[j_eq]); std::swap(iw[i_eq], iw[j_eq]); i_eq++;
        }
        j_eq++;
      }
      int64_t weq = 0;
      for (size_t idx = i; idx < i_eq; ++idx) weq += iw[idx];
      if (wleft + weq > t) return pivot;
      else { t -= (wleft + weq); l = i_eq; }
    }
  }
  return a[l];
}

int main() {
    std::cout << "Starting test..." << std::endl;
    // Just a basic test to see if it links and runs
    return 0;
}
