#ifndef ROBSCALE_WORKER_COMPAT_H
#define ROBSCALE_WORKER_COMPAT_H

#include <RcppParallel.h>
#ifdef USE_DIRECT_TBB
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
