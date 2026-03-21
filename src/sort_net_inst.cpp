// Explicit template instantiation for sort_net.h entry points.
// Defining ROBSCALE_SORT_NET_INST suppresses the extern template declarations
// at the bottom of sort_net.h, so this TU provides the one canonical compiled
// copy of small_sort and median_net for each required type.
#define ROBSCALE_SORT_NET_INST
#include "sort_net.h"

namespace robscale {
template void small_sort<double>(double*, size_t);
template void small_sort<int>(int*, size_t);
template double median_net<double>(double*, size_t);
} // namespace robscale
