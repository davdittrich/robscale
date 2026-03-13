#ifndef ROBSCALE_SORT_NET_H
#define ROBSCALE_SORT_NET_H

#include <algorithm>

namespace robscale {

#define SWAP_IF_GREATER(a, b) \
  do {                        \
    auto min_val = (std::min)((a), (b)); \
    auto max_val = (std::max)((a), (b)); \
    (a) = min_val;            \
    (b) = max_val;            \
  } while (0)

// Optimal sorting networks for small n.
// References: Knuth TAOCP Vol 3, optimal networks from sorting-network.org

template <typename T>
inline void sort_net_2(T* x) {
  SWAP_IF_GREATER(x[0], x[1]);
}

template <typename T>
inline void sort_net_3(T* x) {
  SWAP_IF_GREATER(x[0], x[1]);
  SWAP_IF_GREATER(x[0], x[2]);
  SWAP_IF_GREATER(x[1], x[2]);
}

template <typename T>
inline void sort_net_4(T* x) {
  SWAP_IF_GREATER(x[0], x[1]);
  SWAP_IF_GREATER(x[2], x[3]);
  SWAP_IF_GREATER(x[0], x[2]);
  SWAP_IF_GREATER(x[1], x[3]);
  SWAP_IF_GREATER(x[1], x[2]);
}

template <typename T>
inline void sort_net_5(T* x) {
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
template <typename T>
inline void sort_net_6(T* x) {
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
// Source: Bert Dobbelaere's verified optimal sorting networks
template <typename T>
inline void sort_net_7(T* x) {
  SWAP_IF_GREATER(x[0], x[6]);
  SWAP_IF_GREATER(x[2], x[3]);
  SWAP_IF_GREATER(x[4], x[5]);
  SWAP_IF_GREATER(x[0], x[2]);
  SWAP_IF_GREATER(x[1], x[4]);
  SWAP_IF_GREATER(x[3], x[6]);
  SWAP_IF_GREATER(x[0], x[1]);
  SWAP_IF_GREATER(x[2], x[5]);
  SWAP_IF_GREATER(x[3], x[4]);
  SWAP_IF_GREATER(x[1], x[2]);
  SWAP_IF_GREATER(x[4], x[6]);
  SWAP_IF_GREATER(x[2], x[3]);
  SWAP_IF_GREATER(x[4], x[5]);
  SWAP_IF_GREATER(x[1], x[2]);
  SWAP_IF_GREATER(x[3], x[4]);
  SWAP_IF_GREATER(x[5], x[6]);
}

// 19-comparator optimal network for n=8
template <typename T>
inline void sort_net_8(T* x) {
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

// Sorting networks for n=9..16
// Source: Bert Dobbelaere's verified optimal sorting networks
// https://bertdobbelaere.github.io/sorting_networks.html

template <typename T>
inline void sort_net_9(T* x) { // 25 comparators
  SWAP_IF_GREATER(x[0],x[3]); SWAP_IF_GREATER(x[1],x[7]); SWAP_IF_GREATER(x[2],x[5]); SWAP_IF_GREATER(x[4],x[8]);
  SWAP_IF_GREATER(x[0],x[7]); SWAP_IF_GREATER(x[2],x[4]); SWAP_IF_GREATER(x[3],x[8]); SWAP_IF_GREATER(x[5],x[6]);
  SWAP_IF_GREATER(x[0],x[2]); SWAP_IF_GREATER(x[1],x[3]); SWAP_IF_GREATER(x[4],x[5]); SWAP_IF_GREATER(x[7],x[8]);
  SWAP_IF_GREATER(x[1],x[4]); SWAP_IF_GREATER(x[3],x[6]); SWAP_IF_GREATER(x[5],x[7]);
  SWAP_IF_GREATER(x[0],x[1]); SWAP_IF_GREATER(x[2],x[4]); SWAP_IF_GREATER(x[3],x[5]); SWAP_IF_GREATER(x[6],x[8]);
  SWAP_IF_GREATER(x[2],x[3]); SWAP_IF_GREATER(x[4],x[5]); SWAP_IF_GREATER(x[6],x[7]);
  SWAP_IF_GREATER(x[1],x[2]); SWAP_IF_GREATER(x[3],x[4]); SWAP_IF_GREATER(x[5],x[6]);
}

template <typename T>
inline void sort_net_10(T* x) { // 29 comparators
  SWAP_IF_GREATER(x[0],x[8]); SWAP_IF_GREATER(x[1],x[9]); SWAP_IF_GREATER(x[2],x[7]); SWAP_IF_GREATER(x[3],x[5]); SWAP_IF_GREATER(x[4],x[6]);
  SWAP_IF_GREATER(x[0],x[2]); SWAP_IF_GREATER(x[1],x[4]); SWAP_IF_GREATER(x[5],x[8]); SWAP_IF_GREATER(x[7],x[9]);
  SWAP_IF_GREATER(x[0],x[3]); SWAP_IF_GREATER(x[2],x[4]); SWAP_IF_GREATER(x[5],x[7]); SWAP_IF_GREATER(x[6],x[9]);
  SWAP_IF_GREATER(x[0],x[1]); SWAP_IF_GREATER(x[3],x[6]); SWAP_IF_GREATER(x[8],x[9]);
  SWAP_IF_GREATER(x[1],x[5]); SWAP_IF_GREATER(x[2],x[3]); SWAP_IF_GREATER(x[4],x[8]); SWAP_IF_GREATER(x[6],x[7]);
  SWAP_IF_GREATER(x[1],x[2]); SWAP_IF_GREATER(x[3],x[5]); SWAP_IF_GREATER(x[4],x[6]); SWAP_IF_GREATER(x[7],x[8]);
  SWAP_IF_GREATER(x[2],x[3]); SWAP_IF_GREATER(x[4],x[5]); SWAP_IF_GREATER(x[6],x[7]);
  SWAP_IF_GREATER(x[3],x[4]); SWAP_IF_GREATER(x[5],x[6]);
}

template <typename T>
inline void sort_net_11(T* x) { // 35 comparators
  SWAP_IF_GREATER(x[0],x[9]); SWAP_IF_GREATER(x[1],x[6]); SWAP_IF_GREATER(x[2],x[4]); SWAP_IF_GREATER(x[3],x[7]); SWAP_IF_GREATER(x[5],x[8]);
  SWAP_IF_GREATER(x[0],x[1]); SWAP_IF_GREATER(x[3],x[5]); SWAP_IF_GREATER(x[4],x[10]); SWAP_IF_GREATER(x[6],x[9]); SWAP_IF_GREATER(x[7],x[8]);
  SWAP_IF_GREATER(x[1],x[3]); SWAP_IF_GREATER(x[2],x[5]); SWAP_IF_GREATER(x[4],x[7]); SWAP_IF_GREATER(x[8],x[10]);
  SWAP_IF_GREATER(x[0],x[4]); SWAP_IF_GREATER(x[1],x[2]); SWAP_IF_GREATER(x[3],x[7]); SWAP_IF_GREATER(x[5],x[9]); SWAP_IF_GREATER(x[6],x[8]);
  SWAP_IF_GREATER(x[0],x[1]); SWAP_IF_GREATER(x[2],x[6]); SWAP_IF_GREATER(x[4],x[5]); SWAP_IF_GREATER(x[7],x[8]); SWAP_IF_GREATER(x[9],x[10]);
  SWAP_IF_GREATER(x[2],x[4]); SWAP_IF_GREATER(x[3],x[6]); SWAP_IF_GREATER(x[5],x[7]); SWAP_IF_GREATER(x[8],x[9]);
  SWAP_IF_GREATER(x[1],x[2]); SWAP_IF_GREATER(x[3],x[4]); SWAP_IF_GREATER(x[5],x[6]); SWAP_IF_GREATER(x[7],x[8]);
  SWAP_IF_GREATER(x[2],x[3]); SWAP_IF_GREATER(x[4],x[5]); SWAP_IF_GREATER(x[6],x[7]);
}

template <typename T>
inline void sort_net_12(T* x) { // 39 comparators
  SWAP_IF_GREATER(x[0],x[8]); SWAP_IF_GREATER(x[1],x[7]); SWAP_IF_GREATER(x[2],x[6]); SWAP_IF_GREATER(x[3],x[11]); SWAP_IF_GREATER(x[4],x[10]); SWAP_IF_GREATER(x[5],x[9]);
  SWAP_IF_GREATER(x[0],x[1]); SWAP_IF_GREATER(x[2],x[5]); SWAP_IF_GREATER(x[3],x[4]); SWAP_IF_GREATER(x[6],x[9]); SWAP_IF_GREATER(x[7],x[8]); SWAP_IF_GREATER(x[10],x[11]);
  SWAP_IF_GREATER(x[0],x[2]); SWAP_IF_GREATER(x[1],x[6]); SWAP_IF_GREATER(x[5],x[10]); SWAP_IF_GREATER(x[9],x[11]);
  SWAP_IF_GREATER(x[0],x[3]); SWAP_IF_GREATER(x[1],x[2]); SWAP_IF_GREATER(x[4],x[6]); SWAP_IF_GREATER(x[5],x[7]); SWAP_IF_GREATER(x[8],x[11]); SWAP_IF_GREATER(x[9],x[10]);
  SWAP_IF_GREATER(x[1],x[4]); SWAP_IF_GREATER(x[3],x[5]); SWAP_IF_GREATER(x[6],x[8]); SWAP_IF_GREATER(x[7],x[10]);
  SWAP_IF_GREATER(x[1],x[3]); SWAP_IF_GREATER(x[2],x[5]); SWAP_IF_GREATER(x[6],x[9]); SWAP_IF_GREATER(x[8],x[10]);
  SWAP_IF_GREATER(x[2],x[3]); SWAP_IF_GREATER(x[4],x[5]); SWAP_IF_GREATER(x[6],x[7]); SWAP_IF_GREATER(x[8],x[9]);
  SWAP_IF_GREATER(x[4],x[6]); SWAP_IF_GREATER(x[5],x[7]);
  SWAP_IF_GREATER(x[3],x[4]); SWAP_IF_GREATER(x[5],x[6]); SWAP_IF_GREATER(x[7],x[8]);
}

template <typename T>
inline void sort_net_13(T* x) { // 45 comparators
  SWAP_IF_GREATER(x[0],x[12]); SWAP_IF_GREATER(x[1],x[10]); SWAP_IF_GREATER(x[2],x[9]); SWAP_IF_GREATER(x[3],x[7]); SWAP_IF_GREATER(x[5],x[11]); SWAP_IF_GREATER(x[6],x[8]);
  SWAP_IF_GREATER(x[1],x[6]); SWAP_IF_GREATER(x[2],x[3]); SWAP_IF_GREATER(x[4],x[11]); SWAP_IF_GREATER(x[7],x[9]); SWAP_IF_GREATER(x[8],x[10]);
  SWAP_IF_GREATER(x[0],x[4]); SWAP_IF_GREATER(x[1],x[2]); SWAP_IF_GREATER(x[3],x[6]); SWAP_IF_GREATER(x[7],x[8]); SWAP_IF_GREATER(x[9],x[10]); SWAP_IF_GREATER(x[11],x[12]);
  SWAP_IF_GREATER(x[4],x[6]); SWAP_IF_GREATER(x[5],x[9]); SWAP_IF_GREATER(x[8],x[11]); SWAP_IF_GREATER(x[10],x[12]);
  SWAP_IF_GREATER(x[0],x[5]); SWAP_IF_GREATER(x[3],x[8]); SWAP_IF_GREATER(x[4],x[7]); SWAP_IF_GREATER(x[6],x[11]); SWAP_IF_GREATER(x[9],x[10]);
  SWAP_IF_GREATER(x[0],x[1]); SWAP_IF_GREATER(x[2],x[5]); SWAP_IF_GREATER(x[6],x[9]); SWAP_IF_GREATER(x[7],x[8]); SWAP_IF_GREATER(x[10],x[11]);
  SWAP_IF_GREATER(x[1],x[3]); SWAP_IF_GREATER(x[2],x[4]); SWAP_IF_GREATER(x[5],x[6]); SWAP_IF_GREATER(x[9],x[10]);
  SWAP_IF_GREATER(x[1],x[2]); SWAP_IF_GREATER(x[3],x[4]); SWAP_IF_GREATER(x[5],x[7]); SWAP_IF_GREATER(x[6],x[8]);
  SWAP_IF_GREATER(x[2],x[3]); SWAP_IF_GREATER(x[4],x[5]); SWAP_IF_GREATER(x[6],x[7]); SWAP_IF_GREATER(x[8],x[9]);
  SWAP_IF_GREATER(x[3],x[4]); SWAP_IF_GREATER(x[5],x[6]);
}

template <typename T>
inline void sort_net_14(T* x) { // 51 comparators
  SWAP_IF_GREATER(x[0],x[1]); SWAP_IF_GREATER(x[2],x[3]); SWAP_IF_GREATER(x[4],x[5]); SWAP_IF_GREATER(x[6],x[7]); SWAP_IF_GREATER(x[8],x[9]); SWAP_IF_GREATER(x[10],x[11]); SWAP_IF_GREATER(x[12],x[13]);
  SWAP_IF_GREATER(x[0],x[2]); SWAP_IF_GREATER(x[1],x[3]); SWAP_IF_GREATER(x[4],x[8]); SWAP_IF_GREATER(x[5],x[9]); SWAP_IF_GREATER(x[10],x[12]); SWAP_IF_GREATER(x[11],x[13]);
  SWAP_IF_GREATER(x[0],x[4]); SWAP_IF_GREATER(x[1],x[2]); SWAP_IF_GREATER(x[3],x[7]); SWAP_IF_GREATER(x[5],x[8]); SWAP_IF_GREATER(x[6],x[10]); SWAP_IF_GREATER(x[9],x[13]); SWAP_IF_GREATER(x[11],x[12]);
  SWAP_IF_GREATER(x[0],x[6]); SWAP_IF_GREATER(x[1],x[5]); SWAP_IF_GREATER(x[3],x[9]); SWAP_IF_GREATER(x[4],x[10]); SWAP_IF_GREATER(x[7],x[13]); SWAP_IF_GREATER(x[8],x[12]);
  SWAP_IF_GREATER(x[2],x[10]); SWAP_IF_GREATER(x[3],x[11]); SWAP_IF_GREATER(x[4],x[6]); SWAP_IF_GREATER(x[7],x[9]);
  SWAP_IF_GREATER(x[1],x[3]); SWAP_IF_GREATER(x[2],x[8]); SWAP_IF_GREATER(x[5],x[11]); SWAP_IF_GREATER(x[6],x[7]); SWAP_IF_GREATER(x[10],x[12]);
  SWAP_IF_GREATER(x[1],x[4]); SWAP_IF_GREATER(x[2],x[6]); SWAP_IF_GREATER(x[3],x[5]); SWAP_IF_GREATER(x[7],x[11]); SWAP_IF_GREATER(x[8],x[10]); SWAP_IF_GREATER(x[9],x[12]);
  SWAP_IF_GREATER(x[2],x[4]); SWAP_IF_GREATER(x[3],x[6]); SWAP_IF_GREATER(x[5],x[8]); SWAP_IF_GREATER(x[7],x[10]); SWAP_IF_GREATER(x[9],x[11]);
  SWAP_IF_GREATER(x[3],x[4]); SWAP_IF_GREATER(x[5],x[6]); SWAP_IF_GREATER(x[7],x[8]); SWAP_IF_GREATER(x[9],x[10]);
  SWAP_IF_GREATER(x[6],x[7]);
}

template <typename T>
inline void sort_net_15(T* x) { // 56 comparators
  SWAP_IF_GREATER(x[1],x[2]); SWAP_IF_GREATER(x[3],x[10]); SWAP_IF_GREATER(x[4],x[14]); SWAP_IF_GREATER(x[5],x[8]); SWAP_IF_GREATER(x[6],x[13]); SWAP_IF_GREATER(x[7],x[12]); SWAP_IF_GREATER(x[9],x[11]);
  SWAP_IF_GREATER(x[0],x[14]); SWAP_IF_GREATER(x[1],x[5]); SWAP_IF_GREATER(x[2],x[8]); SWAP_IF_GREATER(x[3],x[7]); SWAP_IF_GREATER(x[6],x[9]); SWAP_IF_GREATER(x[10],x[12]); SWAP_IF_GREATER(x[11],x[13]);
  SWAP_IF_GREATER(x[0],x[7]); SWAP_IF_GREATER(x[1],x[6]); SWAP_IF_GREATER(x[2],x[9]); SWAP_IF_GREATER(x[4],x[10]); SWAP_IF_GREATER(x[5],x[11]); SWAP_IF_GREATER(x[8],x[13]); SWAP_IF_GREATER(x[12],x[14]);
  SWAP_IF_GREATER(x[0],x[6]); SWAP_IF_GREATER(x[2],x[4]); SWAP_IF_GREATER(x[3],x[5]); SWAP_IF_GREATER(x[7],x[11]); SWAP_IF_GREATER(x[8],x[10]); SWAP_IF_GREATER(x[9],x[12]); SWAP_IF_GREATER(x[13],x[14]);
  SWAP_IF_GREATER(x[0],x[3]); SWAP_IF_GREATER(x[1],x[2]); SWAP_IF_GREATER(x[4],x[7]); SWAP_IF_GREATER(x[5],x[9]); SWAP_IF_GREATER(x[6],x[8]); SWAP_IF_GREATER(x[10],x[11]); SWAP_IF_GREATER(x[12],x[13]);
  SWAP_IF_GREATER(x[0],x[1]); SWAP_IF_GREATER(x[2],x[3]); SWAP_IF_GREATER(x[4],x[6]); SWAP_IF_GREATER(x[7],x[9]); SWAP_IF_GREATER(x[10],x[12]); SWAP_IF_GREATER(x[11],x[13]);
  SWAP_IF_GREATER(x[1],x[2]); SWAP_IF_GREATER(x[3],x[5]); SWAP_IF_GREATER(x[8],x[10]); SWAP_IF_GREATER(x[11],x[12]);
  SWAP_IF_GREATER(x[3],x[4]); SWAP_IF_GREATER(x[5],x[6]); SWAP_IF_GREATER(x[7],x[8]); SWAP_IF_GREATER(x[9],x[10]);
  SWAP_IF_GREATER(x[2],x[3]); SWAP_IF_GREATER(x[4],x[5]); SWAP_IF_GREATER(x[6],x[7]); SWAP_IF_GREATER(x[8],x[9]); SWAP_IF_GREATER(x[10],x[11]);
  SWAP_IF_GREATER(x[5],x[6]); SWAP_IF_GREATER(x[7],x[8]);
}

template <typename T>
inline void sort_net_16(T* x) { // 60 comparators
  SWAP_IF_GREATER(x[0],x[13]); SWAP_IF_GREATER(x[1],x[12]); SWAP_IF_GREATER(x[2],x[15]); SWAP_IF_GREATER(x[3],x[14]); SWAP_IF_GREATER(x[4],x[8]); SWAP_IF_GREATER(x[5],x[6]); SWAP_IF_GREATER(x[7],x[11]); SWAP_IF_GREATER(x[9],x[10]);
  SWAP_IF_GREATER(x[0],x[5]); SWAP_IF_GREATER(x[1],x[7]); SWAP_IF_GREATER(x[2],x[9]); SWAP_IF_GREATER(x[3],x[4]); SWAP_IF_GREATER(x[6],x[13]); SWAP_IF_GREATER(x[8],x[14]); SWAP_IF_GREATER(x[10],x[15]); SWAP_IF_GREATER(x[11],x[12]);
  SWAP_IF_GREATER(x[0],x[1]); SWAP_IF_GREATER(x[2],x[3]); SWAP_IF_GREATER(x[4],x[5]); SWAP_IF_GREATER(x[6],x[8]); SWAP_IF_GREATER(x[7],x[9]); SWAP_IF_GREATER(x[10],x[11]); SWAP_IF_GREATER(x[12],x[13]); SWAP_IF_GREATER(x[14],x[15]);
  SWAP_IF_GREATER(x[0],x[2]); SWAP_IF_GREATER(x[1],x[3]); SWAP_IF_GREATER(x[4],x[10]); SWAP_IF_GREATER(x[5],x[11]); SWAP_IF_GREATER(x[6],x[7]); SWAP_IF_GREATER(x[8],x[9]); SWAP_IF_GREATER(x[12],x[14]); SWAP_IF_GREATER(x[13],x[15]);
  SWAP_IF_GREATER(x[1],x[2]); SWAP_IF_GREATER(x[3],x[12]); SWAP_IF_GREATER(x[4],x[6]); SWAP_IF_GREATER(x[5],x[7]); SWAP_IF_GREATER(x[8],x[10]); SWAP_IF_GREATER(x[9],x[11]); SWAP_IF_GREATER(x[13],x[14]);
  SWAP_IF_GREATER(x[1],x[4]); SWAP_IF_GREATER(x[2],x[6]); SWAP_IF_GREATER(x[5],x[8]); SWAP_IF_GREATER(x[7],x[10]); SWAP_IF_GREATER(x[9],x[13]); SWAP_IF_GREATER(x[11],x[14]);
  SWAP_IF_GREATER(x[2],x[4]); SWAP_IF_GREATER(x[3],x[6]); SWAP_IF_GREATER(x[9],x[12]); SWAP_IF_GREATER(x[11],x[13]);
  SWAP_IF_GREATER(x[3],x[5]); SWAP_IF_GREATER(x[6],x[8]); SWAP_IF_GREATER(x[7],x[9]); SWAP_IF_GREATER(x[10],x[12]);
  SWAP_IF_GREATER(x[3],x[4]); SWAP_IF_GREATER(x[5],x[6]); SWAP_IF_GREATER(x[7],x[8]); SWAP_IF_GREATER(x[9],x[10]); SWAP_IF_GREATER(x[11],x[12]);
  SWAP_IF_GREATER(x[6],x[7]); SWAP_IF_GREATER(x[8],x[9]);
}

// Dispatcher: use sorting networks for n<=16, std::sort for n>16
template <typename T>
inline void small_sort(T* x, size_t n) {
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
    case 9:
      sort_net_9(x);
      return;
    case 10:
      sort_net_10(x);
      return;
    case 11:
      sort_net_11(x);
      return;
    case 12:
      sort_net_12(x);
      return;
    case 13:
      sort_net_13(x);
      return;
    case 14:
      sort_net_14(x);
      return;
    case 15:
      sort_net_15(x);
      return;
    case 16:
      sort_net_16(x);
      return;
    default:
      std::sort(x, x + n);
      return;
  }
}

#undef SWAP_IF_GREATER

} // namespace robscale

#endif // ROBSCALE_SORT_NET_H
