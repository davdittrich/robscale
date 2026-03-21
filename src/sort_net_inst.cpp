// Explicit template instantiation for sort_net.h entry points.
// Defining ROBSCALE_SORT_NET_INST suppresses the extern template declarations
// at the bottom of sort_net.h, so this TU provides the one canonical compiled
// copy of small_sort and median_net for each required type.
//
// OPT-B (2026-03-21): visibility("hidden") on median_net<double> and
// small_sort<double> eliminates PLT indirect-call overhead on Linux -fPIC.
// Without this, every call through the extern template routes via the PLT
// (~3-5 ns per call on modern CPUs) because weak symbols in shared objects
// are routed via the dynamic linker even for intra-DSO calls.
// Effect: ~2 calls per rob_scale_core saved ~6-10 ns total.
#define ROBSCALE_SORT_NET_INST
#include "sort_net.h"

#pragma GCC visibility push(hidden)
namespace robscale {
template void small_sort<double>(double*, size_t);
template void small_sort<int>(int*, size_t);
template double median_net<double>(double*, size_t);
} // namespace robscale
#pragma GCC visibility pop
