#ifndef ROBSCALE_SELECTION_H
#define ROBSCALE_SELECTION_H

#include <algorithm>
#include <cmath>
#include <iterator>

namespace robscale {

template <typename Iter>
void floyd_rivest_select(Iter left_in, Iter k, Iter right_in) {
  using T = typename std::iterator_traits<Iter>::value_type;
  size_t n = static_cast<size_t>(std::distance(left_in, right_in));
  if (n < 600) {
    std::nth_element(left_in, k, right_in);
    return;
  }

  Iter left = left_in;
  Iter right = right_in - 1;

  while (right > left) {
    if (right - left > 600) {
      size_t nn = static_cast<size_t>(right - left + 1);
      size_t i = static_cast<size_t>(std::distance(left, k) + 1);
      double z = std::log(static_cast<double>(nn));
      double s = 0.5 * std::exp(2.0 * z / 3.0);
      double sd = 0.5 *
          std::sqrt(z * s * (static_cast<double>(nn) - s) /
                    static_cast<double>(nn)) *
          (static_cast<double>(i) - static_cast<double>(nn) / 2.0 >= 0 ? 1.0
                                                                        : -1.0);
      Iter new_left =
          (std::max)(left,
                     k - static_cast<ptrdiff_t>(static_cast<double>(i) * s /
                                                      static_cast<double>(nn) +
                                                  sd));
      Iter new_right =
          (std::min)(right,
                     k + static_cast<ptrdiff_t>(static_cast<double>(nn - i) * s /
                                                      static_cast<double>(nn) +
                                                  sd));
      floyd_rivest_select(new_left, k, new_right + 1);
    }

    T pivot = *k;
    Iter i = left;
    Iter j = right;
    std::swap(*left, *k);
    if (*right > pivot) std::swap(*left, *right);

    while (i < j) {
      std::swap(*i, *j);
      ++i;
      --j;
      while (*i < pivot) ++i;
      while (*j > pivot) --j;
    }

    if (*left == pivot) {
      std::swap(*left, *j);
    } else {
      ++j;
      std::swap(*j, *right);
    }

    if (j <= k) left = j + 1;
    if (k <= j) right = j - 1;
  }
}

} // namespace robscale

#endif // ROBSCALE_SELECTION_H
