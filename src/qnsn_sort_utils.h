#ifndef ROBSCALE_QNSN_SORT_UTILS_H
#define ROBSCALE_QNSN_SORT_UTILS_H

#include "robscale_config.h"
#include "qnsn_runtime_config.h"
#include <algorithm>
#include <boost/sort/spreadsort/spreadsort.hpp>
// RcppParallel.h must be included before any system TBB headers so that
// RcppParallel's bundled TBB 2019 include guards fire first, preventing
// symbol clashes when the system oneTBB headers are subsequently included.
#include <RcppParallel.h>
#if defined(ROBSCALE_HAS_SYSTEM_TBB)
#include <oneapi/tbb/parallel_sort.h>
#elif defined(USE_DIRECT_TBB)
#include <tbb/parallel_sort.h>
#endif
#include <type_traits>

namespace robscale::qnsn {

template <typename Iterator> void optimized_sort(Iterator begin, Iterator end) {
  using T = typename std::iterator_traits<Iterator>::value_type;
  auto& config = RuntimeConfig::get();
  size_t n = static_cast<size_t>(std::distance(begin, end));

  if (n < 2)
    return;

  if constexpr (std::is_floating_point_v<T>) {
    if (n <= config.sort_boost_threshold) {
      std::sort(begin, end);
    } else if (n < config.sort_tbb_threshold) {
      boost::sort::spreadsort::float_sort(begin, end);
    } else {
#if defined(ROBSCALE_HAS_SYSTEM_TBB) || defined(USE_DIRECT_TBB)
      tbb::parallel_sort(begin, end);
#else
      boost::sort::spreadsort::float_sort(begin, end);
#endif
    }
  } else if constexpr (std::is_integral_v<T>) {
    if (n <= config.sort_boost_threshold) {
      std::sort(begin, end);
    } else if (n < config.sort_tbb_threshold) {
      boost::sort::spreadsort::integer_sort(begin, end);
    } else {
#if defined(ROBSCALE_HAS_SYSTEM_TBB) || defined(USE_DIRECT_TBB)
      tbb::parallel_sort(begin, end);
#else
      boost::sort::spreadsort::integer_sort(begin, end);
#endif
    }
  } else {
    if (n < config.sort_tbb_threshold) {
      std::sort(begin, end);
    } else {
#if defined(ROBSCALE_HAS_SYSTEM_TBB) || defined(USE_DIRECT_TBB)
      tbb::parallel_sort(begin, end);
#else
      std::sort(begin, end);
#endif
    }
  }
}

} // namespace robscale::qnsn

#endif // ROBSCALE_QNSN_SORT_UTILS_H
