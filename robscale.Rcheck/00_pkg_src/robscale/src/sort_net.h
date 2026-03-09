#ifndef ROBSCALE_SORT_NET_H
#define ROBSCALE_SORT_NET_H

#include <algorithm>

namespace robscale {

#define SWAP_IF_GREATER(a, b) \
  do {                        \
    if ((a) > (b)) {          \
      double tmp_ = (a);      \
      (a) = (b);              \
      (b) = tmp_;             \
    }                         \
  } while (0)

// Optimal sorting networks for small n.
// References: Knuth TAOCP Vol 3, optimal networks from sorting-network.org

inline void sort_net_2(double* x) {
  SWAP_IF_GREATER(x[0], x[1]);
}

inline void sort_net_3(double* x) {
  SWAP_IF_GREATER(x[0], x[1]);
  SWAP_IF_GREATER(x[0], x[2]);
  SWAP_IF_GREATER(x[1], x[2]);
}

inline void sort_net_4(double* x) {
  SWAP_IF_GREATER(x[0], x[1]);
  SWAP_IF_GREATER(x[2], x[3]);
  SWAP_IF_GREATER(x[0], x[2]);
  SWAP_IF_GREATER(x[1], x[3]);
  SWAP_IF_GREATER(x[1], x[2]);
}

inline void sort_net_5(double* x) {
  SWAP_IF_GREATER(x[0], x[1]);
  SWAP_IF_GREATER(x[2], x[3]);
  SWAP_IF_GREATER(x[0], x[2]);
  SWAP_IF_GREATER(x[1], x[3]);
  SWAP_IF_GREATER(x[1], x[2]);
  SWAP_IF_GREATER(x[0], x[4]);
  SWAP_IF_GREATER(x[2], x[4]);
  SWAP_IF_GREATER(x[1], x[2]);
  SWAP_IF_GREATER(x[3], x[4]);
}

// 12-comparator optimal network for n=6 (verified: Knuth TAOCP 5.3.4)
inline void sort_net_6(double* x) {
  SWAP_IF_GREATER(x[0], x[5]);
  SWAP_IF_GREATER(x[1], x[3]);
  SWAP_IF_GREATER(x[2], x[4]);
  SWAP_IF_GREATER(x[1], x[2]);
  SWAP_IF_GREATER(x[3], x[4]);
  SWAP_IF_GREATER(x[0], x[3]);
  SWAP_IF_GREATER(x[2], x[5]);
  SWAP_IF_GREATER(x[0], x[1]);
  SWAP_IF_GREATER(x[2], x[3]);
  SWAP_IF_GREATER(x[4], x[5]);
  SWAP_IF_GREATER(x[1], x[2]);
  SWAP_IF_GREATER(x[3], x[4]);
}

// 16-comparator optimal network for n=7
inline void sort_net_7(double* x) {
  SWAP_IF_GREATER(x[0], x[1]);
  SWAP_IF_GREATER(x[2], x[3]);
  SWAP_IF_GREATER(x[4], x[5]);
  SWAP_IF_GREATER(x[0], x[2]);
  SWAP_IF_GREATER(x[1], x[3]);
  SWAP_IF_GREATER(x[4], x[6]);
  SWAP_IF_GREATER(x[5], x[6]);
  SWAP_IF_GREATER(x[0], x[4]);
  SWAP_IF_GREATER(x[1], x[5]);
  SWAP_IF_GREATER(x[2], x[6]);
  SWAP_IF_GREATER(x[2], x[4]);
  SWAP_IF_GREATER(x[3], x[5]);
  SWAP_IF_GREATER(x[1], x[2]);
  SWAP_IF_GREATER(x[3], x[4]);
  SWAP_IF_GREATER(x[5], x[6]);
  SWAP_IF_GREATER(x[2], x[3]);
}

// 19-comparator optimal network for n=8
inline void sort_net_8(double* x) {
  SWAP_IF_GREATER(x[0], x[1]);
  SWAP_IF_GREATER(x[2], x[3]);
  SWAP_IF_GREATER(x[4], x[5]);
  SWAP_IF_GREATER(x[6], x[7]);
  SWAP_IF_GREATER(x[0], x[2]);
  SWAP_IF_GREATER(x[1], x[3]);
  SWAP_IF_GREATER(x[4], x[6]);
  SWAP_IF_GREATER(x[5], x[7]);
  SWAP_IF_GREATER(x[1], x[2]);
  SWAP_IF_GREATER(x[5], x[6]);
  SWAP_IF_GREATER(x[0], x[4]);
  SWAP_IF_GREATER(x[1], x[5]);
  SWAP_IF_GREATER(x[2], x[6]);
  SWAP_IF_GREATER(x[3], x[7]);
  SWAP_IF_GREATER(x[2], x[4]);
  SWAP_IF_GREATER(x[3], x[5]);
  SWAP_IF_GREATER(x[1], x[2]);
  SWAP_IF_GREATER(x[3], x[4]);
  SWAP_IF_GREATER(x[5], x[6]);
}

// Dispatcher: use sorting networks for n<=8, std::sort for n>8
inline void small_sort(double* x, size_t n) {
  switch (n) {
    case 0:
    case 1:
      return;
    case 2:
      sort_net_2(x);
      return;
    case 3:
      sort_net_3(x);
      return;
    case 4:
      sort_net_4(x);
      return;
    case 5:
      sort_net_5(x);
      return;
    case 6:
      sort_net_6(x);
      return;
    case 7:
      sort_net_7(x);
      return;
    case 8:
      sort_net_8(x);
      return;
    default:
      std::sort(x, x + n);
      return;
  }
}

#undef SWAP_IF_GREATER

} // namespace robscale

#endif // ROBSCALE_SORT_NET_H
