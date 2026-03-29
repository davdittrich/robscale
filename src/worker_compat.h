#ifndef ROBSCALE_WORKER_COMPAT_H
#define ROBSCALE_WORKER_COMPAT_H

#include <RcppParallel.h>
#if defined(ROBSCALE_HAS_SYSTEM_TBB)
#include <oneapi/tbb/parallel_for.h>
#include <oneapi/tbb/parallel_reduce.h>
#include <oneapi/tbb/blocked_range.h>
struct WorkerBase {};
using SplitType = tbb::split;
#elif defined(USE_DIRECT_TBB)
#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>
#include <tbb/blocked_range.h>
struct WorkerBase {};
using SplitType = tbb::split;
#else
using WorkerBase = RcppParallel::Worker;
using SplitType = RcppParallel::Split;
#endif

#endif // ROBSCALE_WORKER_COMPAT_H
