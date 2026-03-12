#ifndef ROBSCALE_QNSN_SORT_UTILS_H
#define ROBSCALE_QNSN_SORT_UTILS_H

#include "robscale_config.h"
#include "qnsn_runtime_config.h"
#include <algorithm>
#include <boost/sort/spreadsort/spreadsort.hpp>
// Include TBB parallel_sort
#include <tbb/parallel_sort.h>
#include <type_traits>

namespace robscale::qnsn {

template <typename Iterator> void optimized_sort(Iterator begin, Iterator end) {
  using T = typename std::iterator_traits<Iterator>::value_type;
  size_t n = static_cast<size_t>(std::distance(begin, end));

  if (n < 2)
    return;

  if constexpr (std::is_floating_point_v<T>) {
    if (n <= ROBSCALE_SORT_BOOST_THRESHOLD) {
      std::sort(begin, end);
    } else if (n < RuntimeConfig::get().sort_tbb_threshold) {
      boost::sort::spreadsort::float_sort(begin, end);
    } else {
      tbb::parallel_sort(begin, end);
    }
  } else if constexpr (std::is_integral_v<T>) {
    if (n <= ROBSCALE_SORT_BOOST_THRESHOLD) {
      std::sort(begin, end);
    } else if (n < RuntimeConfig::get().sort_tbb_threshold) {
      boost::sort::spreadsort::integer_sort(begin, end);
    } else {
      tbb::parallel_sort(begin, end);
    }
  } else {
    if (n < RuntimeConfig::get().sort_tbb_threshold) {
      std::sort(begin, end);
    } else {
      tbb::parallel_sort(begin, end);
    }
  }
}

} // namespace robscale::qnsn

#endif // ROBSCALE_QNSN_SORT_UTILS_H
