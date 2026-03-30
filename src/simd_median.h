// simd_median.h — AVX2 bitonic sort median kernels for small double arrays.
//
// Implements full bitonic sorts for n = 4, 8, 16, 32, 64 doubles using
// AVX2 _mm256_min_pd/_mm256_max_pd, then extracts the median.
// Guarded by ROBSCALE_HAS_AVX2_DISPATCH; falls back to scalar median_net
// for sizes without a dedicated SIMD kernel.
//
// Each function carries ROBSCALE_TARGET_AVX2 (__attribute__((target("avx2,fma"))))
// so no global -mavx2 is needed.  The binary remains portable; the caller
// activates SIMD paths via CPUID at runtime.

#ifndef ROBSCALE_SIMD_MEDIAN_H
#define ROBSCALE_SIMD_MEDIAN_H

#include "robscale_config.h"
#include "sort_net.h"  // fallback median_net, small_sort
#include <cstring>     // std::memcpy

#ifdef ROBSCALE_HAS_AVX2_DISPATCH
#if defined(__x86_64__) || defined(_M_X64)
#include <immintrin.h>
#endif

namespace robscale {
namespace simd {

// ============================================================================
// Building blocks: in-register sort of 4 doubles, reverse, bitonic merge of 4
// ============================================================================

// Sort 4 doubles in a single YMM register (ascending).
// Network: stage 1 (0,1)(2,3) → stage 2 (0,2)(1,3) → stage 3 (1,2).
ROBSCALE_TARGET_AVX2
static inline __m256d sort_4_reg(__m256d v) {
    __m256d p, lo, hi;
    // Stage 1: compare-swap adjacent pairs (0,1) and (2,3)
    p  = _mm256_permute_pd(v, 0x5);       // {v1,v0,v3,v2}
    lo = _mm256_min_pd(v, p);
    hi = _mm256_max_pd(v, p);
    v  = _mm256_blend_pd(lo, hi, 0xA);    // {min01,max01,min23,max23}
    // Stage 2: compare-swap cross-lane (0,2) and (1,3)
    p  = _mm256_permute4x64_pd(v, 0x4E);  // {v2,v3,v0,v1}
    lo = _mm256_min_pd(v, p);
    hi = _mm256_max_pd(v, p);
    v  = _mm256_blend_pd(lo, hi, 0xC);    // {min02,min13,max02,max13}
    // Stage 3: compare-swap middle pair (1,2) only
    p  = _mm256_permute4x64_pd(v, 0xD8);  // {v0,v2,v1,v3}
    lo = _mm256_min_pd(v, p);
    hi = _mm256_max_pd(v, p);
    v  = _mm256_blend_pd(lo, hi, 0x4);    // {v0,min12,max12,v3}
    return v;
}

// Reverse 4 elements in a YMM register: {v3,v2,v1,v0}.
ROBSCALE_TARGET_AVX2
static inline __m256d reverse_4_reg(__m256d v) {
    return _mm256_permute4x64_pd(v, 0x1B); // 00_01_10_11
}

// Bitonic merge of 4 elements (input must be a bitonic sequence).
// distance-2: (0,2)(1,3)  →  distance-1: (0,1)(2,3)
ROBSCALE_TARGET_AVX2
static inline __m256d bitonic_merge_4_reg(__m256d v) {
    __m256d p, lo, hi;
    p  = _mm256_permute4x64_pd(v, 0x4E);  // {v2,v3,v0,v1}
    lo = _mm256_min_pd(v, p);
    hi = _mm256_max_pd(v, p);
    v  = _mm256_blend_pd(lo, hi, 0xC);
    p  = _mm256_permute_pd(v, 0x5);       // {v1,v0,v3,v2}
    lo = _mm256_min_pd(v, p);
    hi = _mm256_max_pd(v, p);
    v  = _mm256_blend_pd(lo, hi, 0xA);
    return v;
}

// ============================================================================
// Bitonic merge for 8, 16, 32, 64 elements
// ============================================================================

ROBSCALE_TARGET_AVX2
static inline void bitonic_merge_8(__m256d& v0, __m256d& v1) {
    // distance-4: cross-register
    __m256d lo = _mm256_min_pd(v0, v1);
    __m256d hi = _mm256_max_pd(v0, v1);
    v0 = lo; v1 = hi;
    // distance-2 + distance-1 within each register
    v0 = bitonic_merge_4_reg(v0);
    v1 = bitonic_merge_4_reg(v1);
}

ROBSCALE_TARGET_AVX2
static inline void bitonic_merge_16(
    __m256d& v0, __m256d& v1, __m256d& v2, __m256d& v3)
{
    // distance-8
    __m256d lo, hi;
    lo = _mm256_min_pd(v0, v2); hi = _mm256_max_pd(v0, v2);
    v0 = lo; v2 = hi;
    lo = _mm256_min_pd(v1, v3); hi = _mm256_max_pd(v1, v3);
    v1 = lo; v3 = hi;
    bitonic_merge_8(v0, v1);
    bitonic_merge_8(v2, v3);
}

ROBSCALE_TARGET_AVX2
static inline void bitonic_merge_32(
    __m256d& v0, __m256d& v1, __m256d& v2, __m256d& v3,
    __m256d& v4, __m256d& v5, __m256d& v6, __m256d& v7)
{
    // distance-16
    __m256d lo, hi;
    lo = _mm256_min_pd(v0, v4); hi = _mm256_max_pd(v0, v4);
    v0 = lo; v4 = hi;
    lo = _mm256_min_pd(v1, v5); hi = _mm256_max_pd(v1, v5);
    v1 = lo; v5 = hi;
    lo = _mm256_min_pd(v2, v6); hi = _mm256_max_pd(v2, v6);
    v2 = lo; v6 = hi;
    lo = _mm256_min_pd(v3, v7); hi = _mm256_max_pd(v3, v7);
    v3 = lo; v7 = hi;
    bitonic_merge_16(v0, v1, v2, v3);
    bitonic_merge_16(v4, v5, v6, v7);
}

ROBSCALE_TARGET_AVX2
static inline void bitonic_merge_64(
    __m256d& v0,  __m256d& v1,  __m256d& v2,  __m256d& v3,
    __m256d& v4,  __m256d& v5,  __m256d& v6,  __m256d& v7,
    __m256d& v8,  __m256d& v9,  __m256d& v10, __m256d& v11,
    __m256d& v12, __m256d& v13, __m256d& v14, __m256d& v15)
{
    // distance-32
    __m256d lo, hi;
    lo = _mm256_min_pd(v0, v8);  hi = _mm256_max_pd(v0, v8);
    v0 = lo; v8  = hi;
    lo = _mm256_min_pd(v1, v9);  hi = _mm256_max_pd(v1, v9);
    v1 = lo; v9  = hi;
    lo = _mm256_min_pd(v2, v10); hi = _mm256_max_pd(v2, v10);
    v2 = lo; v10 = hi;
    lo = _mm256_min_pd(v3, v11); hi = _mm256_max_pd(v3, v11);
    v3 = lo; v11 = hi;
    lo = _mm256_min_pd(v4, v12); hi = _mm256_max_pd(v4, v12);
    v4 = lo; v12 = hi;
    lo = _mm256_min_pd(v5, v13); hi = _mm256_max_pd(v5, v13);
    v5 = lo; v13 = hi;
    lo = _mm256_min_pd(v6, v14); hi = _mm256_max_pd(v6, v14);
    v6 = lo; v14 = hi;
    lo = _mm256_min_pd(v7, v15); hi = _mm256_max_pd(v7, v15);
    v7 = lo; v15 = hi;
    bitonic_merge_32(v0, v1, v2, v3, v4, v5, v6, v7);
    bitonic_merge_32(v8, v9, v10, v11, v12, v13, v14, v15);
}

// ============================================================================
// Full bitonic sort for 8, 16, 32, 64 doubles
// ============================================================================

// Sort 8 doubles across 2 YMM registers (ascending).
ROBSCALE_TARGET_AVX2
static inline void sort_8(__m256d& v0, __m256d& v1) {
    v0 = sort_4_reg(v0);
    v1 = sort_4_reg(v1);
    // Form bitonic sequence: reverse second half
    v1 = reverse_4_reg(v1);
    bitonic_merge_8(v0, v1);
}

// Sort 16 doubles across 4 YMM registers (ascending).
ROBSCALE_TARGET_AVX2
static inline void sort_16(
    __m256d& v0, __m256d& v1, __m256d& v2, __m256d& v3)
{
    sort_8(v0, v1);
    sort_8(v2, v3);
    // Reverse second 8: swap register pair and reverse each
    __m256d tmp = reverse_4_reg(v3);
    v3 = reverse_4_reg(v2);
    v2 = tmp;
    bitonic_merge_16(v0, v1, v2, v3);
}

// Sort 32 doubles across 8 YMM registers (ascending).
ROBSCALE_TARGET_AVX2
static inline void sort_32(
    __m256d& v0, __m256d& v1, __m256d& v2, __m256d& v3,
    __m256d& v4, __m256d& v5, __m256d& v6, __m256d& v7)
{
    sort_16(v0, v1, v2, v3);
    sort_16(v4, v5, v6, v7);
    // Reverse second 16: swap outer↔inner pairs and reverse each
    __m256d tmp;
    tmp = v4; v4 = reverse_4_reg(v7); v7 = reverse_4_reg(tmp);
    tmp = v5; v5 = reverse_4_reg(v6); v6 = reverse_4_reg(tmp);
    bitonic_merge_32(v0, v1, v2, v3, v4, v5, v6, v7);
}

// Sort 64 doubles across 16 YMM registers (ascending).
// Uses all 16 YMM registers — register spills expected on AVX2.
ROBSCALE_TARGET_AVX2
static void sort_64(
    __m256d& v0,  __m256d& v1,  __m256d& v2,  __m256d& v3,
    __m256d& v4,  __m256d& v5,  __m256d& v6,  __m256d& v7,
    __m256d& v8,  __m256d& v9,  __m256d& v10, __m256d& v11,
    __m256d& v12, __m256d& v13, __m256d& v14, __m256d& v15)
{
    sort_32(v0, v1, v2, v3, v4, v5, v6, v7);
    sort_32(v8, v9, v10, v11, v12, v13, v14, v15);
    // Reverse second 32
    __m256d tmp;
    tmp = v8;  v8  = reverse_4_reg(v15); v15 = reverse_4_reg(tmp);
    tmp = v9;  v9  = reverse_4_reg(v14); v14 = reverse_4_reg(tmp);
    tmp = v10; v10 = reverse_4_reg(v13); v13 = reverse_4_reg(tmp);
    tmp = v11; v11 = reverse_4_reg(v12); v12 = reverse_4_reg(tmp);
    bitonic_merge_64(v0, v1, v2, v3, v4, v5, v6, v7,
                     v8, v9, v10, v11, v12, v13, v14, v15);
}

// ============================================================================
// Median extraction helpers
// ============================================================================

// Extract lane 3 (highest double) from a YMM register.
ROBSCALE_TARGET_AVX2
static inline double extract_hi(__m256d v) {
    __m128d hi = _mm256_extractf128_pd(v, 1);       // {lane2, lane3}
    return _mm_cvtsd_f64(_mm_unpackhi_pd(hi, hi));   // lane3
}

// Extract lane 0 (lowest double) from a YMM register.
ROBSCALE_TARGET_AVX2
static inline double extract_lo(__m256d v) {
    return _mm_cvtsd_f64(_mm256_castpd256_pd128(v));  // lane0
}

// ============================================================================
// SIMD full sort: sort n doubles in-place via bitonic network.
// Pads to the next supported size (4/8/16/32/64) with +inf, sorts, copies back.
// ============================================================================

ROBSCALE_TARGET_AVX2
static void simd_sort_4(double* x, size_t n) {
    alignas(32) double buf[4] = {
        (n > 0) ? x[0] : __builtin_inf(),
        (n > 1) ? x[1] : __builtin_inf(),
        (n > 2) ? x[2] : __builtin_inf(),
        (n > 3) ? x[3] : __builtin_inf()
    };
    __m256d v = _mm256_load_pd(buf);
    v = sort_4_reg(v);
    _mm256_store_pd(buf, v);
    for (size_t i = 0; i < n; ++i) x[i] = buf[i];
}

ROBSCALE_TARGET_AVX2
static void simd_sort_8(double* x, size_t n) {
    alignas(32) double buf[8];
    for (size_t i = 0; i < 8; ++i) buf[i] = (i < n) ? x[i] : __builtin_inf();
    __m256d v0 = _mm256_load_pd(buf);
    __m256d v1 = _mm256_load_pd(buf + 4);
    sort_8(v0, v1);
    _mm256_store_pd(buf, v0);
    _mm256_store_pd(buf + 4, v1);
    for (size_t i = 0; i < n; ++i) x[i] = buf[i];
}

ROBSCALE_TARGET_AVX2
static void simd_sort_16(double* x, size_t n) {
    alignas(32) double buf[16];
    for (size_t i = 0; i < 16; ++i) buf[i] = (i < n) ? x[i] : __builtin_inf();
    __m256d v0 = _mm256_load_pd(buf);
    __m256d v1 = _mm256_load_pd(buf + 4);
    __m256d v2 = _mm256_load_pd(buf + 8);
    __m256d v3 = _mm256_load_pd(buf + 12);
    sort_16(v0, v1, v2, v3);
    _mm256_store_pd(buf, v0);     _mm256_store_pd(buf + 4, v1);
    _mm256_store_pd(buf + 8, v2); _mm256_store_pd(buf + 12, v3);
    for (size_t i = 0; i < n; ++i) x[i] = buf[i];
}

ROBSCALE_TARGET_AVX2
static void simd_sort_32(double* x, size_t n) {
    alignas(32) double buf[32];
    for (size_t i = 0; i < 32; ++i) buf[i] = (i < n) ? x[i] : __builtin_inf();
    __m256d v0 = _mm256_load_pd(buf);      __m256d v1 = _mm256_load_pd(buf + 4);
    __m256d v2 = _mm256_load_pd(buf + 8);  __m256d v3 = _mm256_load_pd(buf + 12);
    __m256d v4 = _mm256_load_pd(buf + 16); __m256d v5 = _mm256_load_pd(buf + 20);
    __m256d v6 = _mm256_load_pd(buf + 24); __m256d v7 = _mm256_load_pd(buf + 28);
    sort_32(v0, v1, v2, v3, v4, v5, v6, v7);
    _mm256_store_pd(buf, v0);      _mm256_store_pd(buf + 4, v1);
    _mm256_store_pd(buf + 8, v2);  _mm256_store_pd(buf + 12, v3);
    _mm256_store_pd(buf + 16, v4); _mm256_store_pd(buf + 20, v5);
    _mm256_store_pd(buf + 24, v6); _mm256_store_pd(buf + 28, v7);
    for (size_t i = 0; i < n; ++i) x[i] = buf[i];
}

ROBSCALE_TARGET_AVX2
static void simd_sort_64(double* x, size_t n) {
    alignas(32) double buf[64];
    for (size_t i = 0; i < 64; ++i) buf[i] = (i < n) ? x[i] : __builtin_inf();
    __m256d v0  = _mm256_load_pd(buf);      __m256d v1  = _mm256_load_pd(buf + 4);
    __m256d v2  = _mm256_load_pd(buf + 8);  __m256d v3  = _mm256_load_pd(buf + 12);
    __m256d v4  = _mm256_load_pd(buf + 16); __m256d v5  = _mm256_load_pd(buf + 20);
    __m256d v6  = _mm256_load_pd(buf + 24); __m256d v7  = _mm256_load_pd(buf + 28);
    __m256d v8  = _mm256_load_pd(buf + 32); __m256d v9  = _mm256_load_pd(buf + 36);
    __m256d v10 = _mm256_load_pd(buf + 40); __m256d v11 = _mm256_load_pd(buf + 44);
    __m256d v12 = _mm256_load_pd(buf + 48); __m256d v13 = _mm256_load_pd(buf + 52);
    __m256d v14 = _mm256_load_pd(buf + 56); __m256d v15 = _mm256_load_pd(buf + 60);
    sort_64(v0, v1, v2, v3, v4, v5, v6, v7,
            v8, v9, v10, v11, v12, v13, v14, v15);
    _mm256_store_pd(buf, v0);      _mm256_store_pd(buf + 4, v1);
    _mm256_store_pd(buf + 8, v2);  _mm256_store_pd(buf + 12, v3);
    _mm256_store_pd(buf + 16, v4); _mm256_store_pd(buf + 20, v5);
    _mm256_store_pd(buf + 24, v6); _mm256_store_pd(buf + 28, v7);
    _mm256_store_pd(buf + 32, v8); _mm256_store_pd(buf + 36, v9);
    _mm256_store_pd(buf + 40, v10);_mm256_store_pd(buf + 44, v11);
    _mm256_store_pd(buf + 48, v12);_mm256_store_pd(buf + 52, v13);
    _mm256_store_pd(buf + 56, v14);_mm256_store_pd(buf + 60, v15);
    for (size_t i = 0; i < n; ++i) x[i] = buf[i];
}

// Dispatch: pick the smallest SIMD sort that fits n, pad with +inf.
ROBSCALE_TARGET_AVX2
static void simd_sort_dispatch_avx2(double* x, size_t n) {
    if      (n <= 4)  simd_sort_4(x, n);
    else if (n <= 8)  simd_sort_8(x, n);
    else if (n <= 16) simd_sort_16(x, n);
    else if (n <= 32) simd_sort_32(x, n);
    else              simd_sort_64(x, n);  // n <= 64 assumed
}

// ============================================================================
// SIMD median functions for even-n arrays
// ============================================================================

// n=4: median = (sorted[1] + sorted[2]) / 2
ROBSCALE_TARGET_AVX2
static double simd_median_4(double* x) {
    __m256d v = _mm256_loadu_pd(x);
    v = sort_4_reg(v);
    __m128d lo128 = _mm256_castpd256_pd128(v);     // {lane0, lane1}
    __m128d hi128 = _mm256_extractf128_pd(v, 1);   // {lane2, lane3}
    double a = _mm_cvtsd_f64(_mm_unpackhi_pd(lo128, lo128)); // lane1
    double b = _mm_cvtsd_f64(hi128);                          // lane2
    return (a + b) * 0.5;
}

// n=8: median = (sorted[3] + sorted[4]) / 2
ROBSCALE_TARGET_AVX2
static double simd_median_8(double* x) {
    __m256d v0 = _mm256_loadu_pd(x);
    __m256d v1 = _mm256_loadu_pd(x + 4);
    sort_8(v0, v1);
    return (extract_hi(v0) + extract_lo(v1)) * 0.5;
}

// n=16: median = (sorted[7] + sorted[8]) / 2
ROBSCALE_TARGET_AVX2
static double simd_median_16(double* x) {
    __m256d v0 = _mm256_loadu_pd(x);
    __m256d v1 = _mm256_loadu_pd(x + 4);
    __m256d v2 = _mm256_loadu_pd(x + 8);
    __m256d v3 = _mm256_loadu_pd(x + 12);
    sort_16(v0, v1, v2, v3);
    return (extract_hi(v1) + extract_lo(v2)) * 0.5;
}

// n=32: median = (sorted[15] + sorted[16]) / 2
ROBSCALE_TARGET_AVX2
static double simd_median_32(double* x) {
    __m256d v0 = _mm256_loadu_pd(x);
    __m256d v1 = _mm256_loadu_pd(x + 4);
    __m256d v2 = _mm256_loadu_pd(x + 8);
    __m256d v3 = _mm256_loadu_pd(x + 12);
    __m256d v4 = _mm256_loadu_pd(x + 16);
    __m256d v5 = _mm256_loadu_pd(x + 20);
    __m256d v6 = _mm256_loadu_pd(x + 24);
    __m256d v7 = _mm256_loadu_pd(x + 28);
    sort_32(v0, v1, v2, v3, v4, v5, v6, v7);
    return (extract_hi(v3) + extract_lo(v4)) * 0.5;
}

// n=64: median = (sorted[31] + sorted[32]) / 2
// All 16 YMM registers consumed — register spills expected.
ROBSCALE_TARGET_AVX2
static double simd_median_64(double* x) {
    __m256d v0  = _mm256_loadu_pd(x);
    __m256d v1  = _mm256_loadu_pd(x + 4);
    __m256d v2  = _mm256_loadu_pd(x + 8);
    __m256d v3  = _mm256_loadu_pd(x + 12);
    __m256d v4  = _mm256_loadu_pd(x + 16);
    __m256d v5  = _mm256_loadu_pd(x + 20);
    __m256d v6  = _mm256_loadu_pd(x + 24);
    __m256d v7  = _mm256_loadu_pd(x + 28);
    __m256d v8  = _mm256_loadu_pd(x + 32);
    __m256d v9  = _mm256_loadu_pd(x + 36);
    __m256d v10 = _mm256_loadu_pd(x + 40);
    __m256d v11 = _mm256_loadu_pd(x + 44);
    __m256d v12 = _mm256_loadu_pd(x + 48);
    __m256d v13 = _mm256_loadu_pd(x + 52);
    __m256d v14 = _mm256_loadu_pd(x + 56);
    __m256d v15 = _mm256_loadu_pd(x + 60);
    sort_64(v0, v1, v2, v3, v4, v5, v6, v7,
            v8, v9, v10, v11, v12, v13, v14, v15);
    return (extract_hi(v7) + extract_lo(v8)) * 0.5;
}

// ============================================================================
// SIMD selection networks: hybrid SIMD (regular early stages) + scalar tail.
//
// These vectorize the EXISTING median selection networks from sort_net.h,
// using SIMD only for stages where independent comparators pack into
// vector lanes. The irregular tail runs scalar on an aligned stack buffer.
// ============================================================================

#define SIMD_SWAP(a, b) do {                 \
    auto mn_ = (std::min)((a), (b));         \
    auto mx_ = (std::max)((a), (b));         \
    (a) = mn_; (b) = mx_;                    \
} while (0)

// n=8 selection network (16 comparators, depth 5).
// Stages A+B (8 comparators): SIMD.  Stages C-E (8 comparators): scalar.
ROBSCALE_TARGET_AVX2
static double simd_median_sel_8(double* x) {
    __m256d v0 = _mm256_loadu_pd(x);
    __m256d v1 = _mm256_loadu_pd(x + 4);

    // Stage A: (0,2),(1,3),(4,6),(5,7) — within-register distance-2
    {
        __m256d p, lo, hi;
        p  = _mm256_permute4x64_pd(v0, 0x4E);
        lo = _mm256_min_pd(v0, p);
        hi = _mm256_max_pd(v0, p);
        v0 = _mm256_blend_pd(lo, hi, 0xC);
        p  = _mm256_permute4x64_pd(v1, 0x4E);
        lo = _mm256_min_pd(v1, p);
        hi = _mm256_max_pd(v1, p);
        v1 = _mm256_blend_pd(lo, hi, 0xC);
    }
    // Stage B: (0,4),(1,5),(2,6),(3,7) — cross-register element-wise
    {
        __m256d lo = _mm256_min_pd(v0, v1);
        __m256d hi = _mm256_max_pd(v0, v1);
        v0 = lo; v1 = hi;
    }
    // Store to aligned buffer for scalar tail
    alignas(32) double b[8];
    _mm256_store_pd(b, v0);
    _mm256_store_pd(b + 4, v1);
    // Stage C: (0,1),(2,4),(3,5),(6,7)
    SIMD_SWAP(b[0], b[1]); SIMD_SWAP(b[2], b[4]);
    SIMD_SWAP(b[3], b[5]); SIMD_SWAP(b[6], b[7]);
    // Stage D: (2,3),(4,5)
    SIMD_SWAP(b[2], b[3]); SIMD_SWAP(b[4], b[5]);
    // Stage E: (1,4),(3,6)
    SIMD_SWAP(b[1], b[4]); SIMD_SWAP(b[3], b[6]);
    return (b[3] + b[4]) * 0.5;
}

// n=16 selection network (46 comparators, depth 10).
// Stage 1 (8 comparators): adjacent-pair SIMD.  Rest (38 comparators): scalar.
ROBSCALE_TARGET_AVX2
static double simd_median_sel_16(double* x) {
    __m256d v0 = _mm256_loadu_pd(x);
    __m256d v1 = _mm256_loadu_pd(x + 4);
    __m256d v2 = _mm256_loadu_pd(x + 8);
    __m256d v3 = _mm256_loadu_pd(x + 12);

    // Stage 1: (0,1),(2,3),(4,5),(6,7),(8,9),(10,11),(12,13),(14,15)
    // Adjacent-pair swap in all 4 registers
    {
        __m256d p, lo, hi;
        p = _mm256_permute_pd(v0, 0x5); lo = _mm256_min_pd(v0, p); hi = _mm256_max_pd(v0, p); v0 = _mm256_blend_pd(lo, hi, 0xA);
        p = _mm256_permute_pd(v1, 0x5); lo = _mm256_min_pd(v1, p); hi = _mm256_max_pd(v1, p); v1 = _mm256_blend_pd(lo, hi, 0xA);
        p = _mm256_permute_pd(v2, 0x5); lo = _mm256_min_pd(v2, p); hi = _mm256_max_pd(v2, p); v2 = _mm256_blend_pd(lo, hi, 0xA);
        p = _mm256_permute_pd(v3, 0x5); lo = _mm256_min_pd(v3, p); hi = _mm256_max_pd(v3, p); v3 = _mm256_blend_pd(lo, hi, 0xA);
    }
    // Store to aligned buffer for scalar tail (38 remaining comparators)
    alignas(32) double b[16];
    _mm256_store_pd(b, v0);
    _mm256_store_pd(b + 4, v1);
    _mm256_store_pd(b + 8, v2);
    _mm256_store_pd(b + 12, v3);
    // Remaining comparators from median_net_16 (lines 3279–3316)
    SIMD_SWAP(b[0],b[6]);   SIMD_SWAP(b[2],b[4]);
    SIMD_SWAP(b[9],b[15]);  SIMD_SWAP(b[11],b[13]);
    SIMD_SWAP(b[4],b[9]);   SIMD_SWAP(b[6],b[11]);
    SIMD_SWAP(b[1],b[9]);   SIMD_SWAP(b[3],b[11]);
    SIMD_SWAP(b[4],b[12]);  SIMD_SWAP(b[6],b[14]);
    SIMD_SWAP(b[0],b[4]);   SIMD_SWAP(b[1],b[10]);
    SIMD_SWAP(b[2],b[6]);   SIMD_SWAP(b[3],b[8]);
    SIMD_SWAP(b[5],b[14]);  SIMD_SWAP(b[7],b[12]);
    SIMD_SWAP(b[9],b[13]);  SIMD_SWAP(b[11],b[15]);
    SIMD_SWAP(b[1],b[12]);  SIMD_SWAP(b[3],b[14]);
    SIMD_SWAP(b[4],b[7]);   SIMD_SWAP(b[5],b[6]);
    SIMD_SWAP(b[8],b[11]);  SIMD_SWAP(b[9],b[10]);
    SIMD_SWAP(b[1],b[3]);   SIMD_SWAP(b[4],b[5]);
    SIMD_SWAP(b[6],b[7]);   SIMD_SWAP(b[8],b[9]);
    SIMD_SWAP(b[10],b[11]); SIMD_SWAP(b[12],b[14]);
    SIMD_SWAP(b[3],b[5]);   SIMD_SWAP(b[6],b[8]);
    SIMD_SWAP(b[7],b[9]);   SIMD_SWAP(b[10],b[12]);
    SIMD_SWAP(b[5],b[8]);   SIMD_SWAP(b[7],b[10]);
    SIMD_SWAP(b[5],b[7]);   SIMD_SWAP(b[8],b[10]);
    return (b[7] + b[8]) * 0.5;
}

// n=32 selection network (128 comparators, depth 15).
// Stages 1-3 (48 comparators): SIMD. Rest (80 comparators): scalar.
ROBSCALE_TARGET_AVX2
static double simd_median_sel_32(double* x) {
    __m256d v0 = _mm256_loadu_pd(x);      __m256d v1 = _mm256_loadu_pd(x + 4);
    __m256d v2 = _mm256_loadu_pd(x + 8);  __m256d v3 = _mm256_loadu_pd(x + 12);
    __m256d v4 = _mm256_loadu_pd(x + 16); __m256d v5 = _mm256_loadu_pd(x + 20);
    __m256d v6 = _mm256_loadu_pd(x + 24); __m256d v7 = _mm256_loadu_pd(x + 28);

    // Stage 1: 16 adjacent-pair comparators
    {
        __m256d p, lo, hi;
        p = _mm256_permute_pd(v0, 0x5); lo = _mm256_min_pd(v0, p); hi = _mm256_max_pd(v0, p); v0 = _mm256_blend_pd(lo, hi, 0xA);
        p = _mm256_permute_pd(v1, 0x5); lo = _mm256_min_pd(v1, p); hi = _mm256_max_pd(v1, p); v1 = _mm256_blend_pd(lo, hi, 0xA);
        p = _mm256_permute_pd(v2, 0x5); lo = _mm256_min_pd(v2, p); hi = _mm256_max_pd(v2, p); v2 = _mm256_blend_pd(lo, hi, 0xA);
        p = _mm256_permute_pd(v3, 0x5); lo = _mm256_min_pd(v3, p); hi = _mm256_max_pd(v3, p); v3 = _mm256_blend_pd(lo, hi, 0xA);
        p = _mm256_permute_pd(v4, 0x5); lo = _mm256_min_pd(v4, p); hi = _mm256_max_pd(v4, p); v4 = _mm256_blend_pd(lo, hi, 0xA);
        p = _mm256_permute_pd(v5, 0x5); lo = _mm256_min_pd(v5, p); hi = _mm256_max_pd(v5, p); v5 = _mm256_blend_pd(lo, hi, 0xA);
        p = _mm256_permute_pd(v6, 0x5); lo = _mm256_min_pd(v6, p); hi = _mm256_max_pd(v6, p); v6 = _mm256_blend_pd(lo, hi, 0xA);
        p = _mm256_permute_pd(v7, 0x5); lo = _mm256_min_pd(v7, p); hi = _mm256_max_pd(v7, p); v7 = _mm256_blend_pd(lo, hi, 0xA);
    }
    // Stage 2: 16 distance-2 comparators — within-register cross-lane
    {
        __m256d p, lo, hi;
        p = _mm256_permute4x64_pd(v0, 0x4E); lo = _mm256_min_pd(v0, p); hi = _mm256_max_pd(v0, p); v0 = _mm256_blend_pd(lo, hi, 0xC);
        p = _mm256_permute4x64_pd(v1, 0x4E); lo = _mm256_min_pd(v1, p); hi = _mm256_max_pd(v1, p); v1 = _mm256_blend_pd(lo, hi, 0xC);
        p = _mm256_permute4x64_pd(v2, 0x4E); lo = _mm256_min_pd(v2, p); hi = _mm256_max_pd(v2, p); v2 = _mm256_blend_pd(lo, hi, 0xC);
        p = _mm256_permute4x64_pd(v3, 0x4E); lo = _mm256_min_pd(v3, p); hi = _mm256_max_pd(v3, p); v3 = _mm256_blend_pd(lo, hi, 0xC);
        p = _mm256_permute4x64_pd(v4, 0x4E); lo = _mm256_min_pd(v4, p); hi = _mm256_max_pd(v4, p); v4 = _mm256_blend_pd(lo, hi, 0xC);
        p = _mm256_permute4x64_pd(v5, 0x4E); lo = _mm256_min_pd(v5, p); hi = _mm256_max_pd(v5, p); v5 = _mm256_blend_pd(lo, hi, 0xC);
        p = _mm256_permute4x64_pd(v6, 0x4E); lo = _mm256_min_pd(v6, p); hi = _mm256_max_pd(v6, p); v6 = _mm256_blend_pd(lo, hi, 0xC);
        p = _mm256_permute4x64_pd(v7, 0x4E); lo = _mm256_min_pd(v7, p); hi = _mm256_max_pd(v7, p); v7 = _mm256_blend_pd(lo, hi, 0xC);
    }
    // Stage 3: 16 distance-4 comparators — cross-register element-wise
    {
        __m256d lo, hi;
        lo = _mm256_min_pd(v0, v1); hi = _mm256_max_pd(v0, v1); v0 = lo; v1 = hi;
        lo = _mm256_min_pd(v2, v3); hi = _mm256_max_pd(v2, v3); v2 = lo; v3 = hi;
        lo = _mm256_min_pd(v4, v5); hi = _mm256_max_pd(v4, v5); v4 = lo; v5 = hi;
        lo = _mm256_min_pd(v6, v7); hi = _mm256_max_pd(v6, v7); v6 = lo; v7 = hi;
    }
    // Store to aligned buffer for scalar tail (80 remaining comparators)
    alignas(32) double b[32];
    _mm256_store_pd(b,      v0); _mm256_store_pd(b + 4,  v1);
    _mm256_store_pd(b + 8,  v2); _mm256_store_pd(b + 12, v3);
    _mm256_store_pd(b + 16, v4); _mm256_store_pd(b + 20, v5);
    _mm256_store_pd(b + 24, v6); _mm256_store_pd(b + 28, v7);
    // Remaining comparators from median_net_32 (lines 4711–4790)
    SIMD_SWAP(b[1],b[6]);   SIMD_SWAP(b[4],b[22]);  SIMD_SWAP(b[5],b[18]);
    SIMD_SWAP(b[9],b[27]);  SIMD_SWAP(b[11],b[19]); SIMD_SWAP(b[12],b[20]);
    SIMD_SWAP(b[13],b[26]); SIMD_SWAP(b[25],b[30]);
    SIMD_SWAP(b[1],b[25]);  SIMD_SWAP(b[2],b[13]);  SIMD_SWAP(b[3],b[27]);
    SIMD_SWAP(b[4],b[28]);  SIMD_SWAP(b[5],b[10]);  SIMD_SWAP(b[6],b[30]);
    SIMD_SWAP(b[9],b[17]);  SIMD_SWAP(b[14],b[22]); SIMD_SWAP(b[18],b[29]);
    SIMD_SWAP(b[21],b[26]);
    SIMD_SWAP(b[0],b[9]);   SIMD_SWAP(b[1],b[12]);  SIMD_SWAP(b[2],b[16]);
    SIMD_SWAP(b[3],b[17]);  SIMD_SWAP(b[4],b[8]);   SIMD_SWAP(b[5],b[24]);
    SIMD_SWAP(b[6],b[11]);  SIMD_SWAP(b[7],b[26]);  SIMD_SWAP(b[10],b[18]);
    SIMD_SWAP(b[13],b[21]); SIMD_SWAP(b[14],b[28]); SIMD_SWAP(b[15],b[29]);
    SIMD_SWAP(b[19],b[30]); SIMD_SWAP(b[20],b[25]); SIMD_SWAP(b[22],b[31]);
    SIMD_SWAP(b[23],b[27]);
    SIMD_SWAP(b[3],b[14]);  SIMD_SWAP(b[6],b[20]);  SIMD_SWAP(b[7],b[23]);
    SIMD_SWAP(b[8],b[24]);  SIMD_SWAP(b[9],b[16]);  SIMD_SWAP(b[10],b[13]);
    SIMD_SWAP(b[11],b[25]); SIMD_SWAP(b[15],b[22]); SIMD_SWAP(b[17],b[28]);
    SIMD_SWAP(b[18],b[21]);
    SIMD_SWAP(b[6],b[9]);   SIMD_SWAP(b[7],b[15]);  SIMD_SWAP(b[8],b[10]);
    SIMD_SWAP(b[11],b[13]); SIMD_SWAP(b[12],b[14]); SIMD_SWAP(b[16],b[24]);
    SIMD_SWAP(b[17],b[19]); SIMD_SWAP(b[18],b[20]); SIMD_SWAP(b[21],b[23]);
    SIMD_SWAP(b[22],b[25]);
    SIMD_SWAP(b[3],b[18]);  SIMD_SWAP(b[7],b[17]);  SIMD_SWAP(b[9],b[16]);
    SIMD_SWAP(b[10],b[12]); SIMD_SWAP(b[13],b[28]); SIMD_SWAP(b[14],b[24]);
    SIMD_SWAP(b[15],b[22]); SIMD_SWAP(b[19],b[21]);
    SIMD_SWAP(b[7],b[20]);  SIMD_SWAP(b[11],b[24]); SIMD_SWAP(b[12],b[16]);
    SIMD_SWAP(b[13],b[17]); SIMD_SWAP(b[14],b[18]); SIMD_SWAP(b[15],b[19]);
    SIMD_SWAP(b[7],b[11]);  SIMD_SWAP(b[13],b[15]); SIMD_SWAP(b[16],b[18]);
    SIMD_SWAP(b[20],b[24]);
    SIMD_SWAP(b[11],b[18]); SIMD_SWAP(b[13],b[20]); SIMD_SWAP(b[14],b[16]);
    SIMD_SWAP(b[15],b[17]);
    SIMD_SWAP(b[11],b[13]); SIMD_SWAP(b[18],b[20]);
    SIMD_SWAP(b[13],b[16]); SIMD_SWAP(b[15],b[18]);
    return (b[15] + b[16]) * 0.5;
}

// ============================================================================
// SIMD hybrid FULL sorting networks: SIMD for regular early stages,
// scalar for the irregular tail.  Sorts the array IN PLACE.
// ============================================================================

// sort_net_8 hybrid: stages 1-2 in SIMD (8 of 19 comparators), rest scalar.
ROBSCALE_TARGET_AVX2
static void simd_sort_net_8(double* x) {
    __m256d v0 = _mm256_loadu_pd(x);
    __m256d v1 = _mm256_loadu_pd(x + 4);
    // Stage 1: adjacent pairs (0,1),(2,3),(4,5),(6,7)
    {
        __m256d p, lo, hi;
        p = _mm256_permute_pd(v0, 0x5); lo = _mm256_min_pd(v0, p); hi = _mm256_max_pd(v0, p); v0 = _mm256_blend_pd(lo, hi, 0xA);
        p = _mm256_permute_pd(v1, 0x5); lo = _mm256_min_pd(v1, p); hi = _mm256_max_pd(v1, p); v1 = _mm256_blend_pd(lo, hi, 0xA);
    }
    // Stage 2: distance-2 (0,2),(1,3),(4,6),(5,7)
    {
        __m256d p, lo, hi;
        p = _mm256_permute4x64_pd(v0, 0x4E); lo = _mm256_min_pd(v0, p); hi = _mm256_max_pd(v0, p); v0 = _mm256_blend_pd(lo, hi, 0xC);
        p = _mm256_permute4x64_pd(v1, 0x4E); lo = _mm256_min_pd(v1, p); hi = _mm256_max_pd(v1, p); v1 = _mm256_blend_pd(lo, hi, 0xC);
    }
    alignas(32) double b[8];
    _mm256_store_pd(b, v0);
    _mm256_store_pd(b + 4, v1);
    // Stages 3-6 scalar (11 remaining comparators)
    SIMD_SWAP(b[1],b[2]); SIMD_SWAP(b[5],b[6]);                     // stage 3
    SIMD_SWAP(b[0],b[4]); SIMD_SWAP(b[1],b[5]);                     // stage 4
    SIMD_SWAP(b[2],b[6]); SIMD_SWAP(b[3],b[7]);
    SIMD_SWAP(b[2],b[4]); SIMD_SWAP(b[3],b[5]);                     // stage 5
    SIMD_SWAP(b[1],b[2]); SIMD_SWAP(b[3],b[4]); SIMD_SWAP(b[5],b[6]); // stage 6
    std::memcpy(x, b, 8 * sizeof(double));
}

// sort_net_32 hybrid: stages 1-4 in SIMD (64 of 185 comparators), rest scalar.
ROBSCALE_TARGET_AVX2
static void simd_sort_net_32(double* x) {
    __m256d v0 = _mm256_loadu_pd(x);      __m256d v1 = _mm256_loadu_pd(x + 4);
    __m256d v2 = _mm256_loadu_pd(x + 8);  __m256d v3 = _mm256_loadu_pd(x + 12);
    __m256d v4 = _mm256_loadu_pd(x + 16); __m256d v5 = _mm256_loadu_pd(x + 20);
    __m256d v6 = _mm256_loadu_pd(x + 24); __m256d v7 = _mm256_loadu_pd(x + 28);
    // Stage 1: 16 adjacent-pair comparators
    {
        __m256d p, lo, hi;
        p = _mm256_permute_pd(v0, 0x5); lo = _mm256_min_pd(v0, p); hi = _mm256_max_pd(v0, p); v0 = _mm256_blend_pd(lo, hi, 0xA);
        p = _mm256_permute_pd(v1, 0x5); lo = _mm256_min_pd(v1, p); hi = _mm256_max_pd(v1, p); v1 = _mm256_blend_pd(lo, hi, 0xA);
        p = _mm256_permute_pd(v2, 0x5); lo = _mm256_min_pd(v2, p); hi = _mm256_max_pd(v2, p); v2 = _mm256_blend_pd(lo, hi, 0xA);
        p = _mm256_permute_pd(v3, 0x5); lo = _mm256_min_pd(v3, p); hi = _mm256_max_pd(v3, p); v3 = _mm256_blend_pd(lo, hi, 0xA);
        p = _mm256_permute_pd(v4, 0x5); lo = _mm256_min_pd(v4, p); hi = _mm256_max_pd(v4, p); v4 = _mm256_blend_pd(lo, hi, 0xA);
        p = _mm256_permute_pd(v5, 0x5); lo = _mm256_min_pd(v5, p); hi = _mm256_max_pd(v5, p); v5 = _mm256_blend_pd(lo, hi, 0xA);
        p = _mm256_permute_pd(v6, 0x5); lo = _mm256_min_pd(v6, p); hi = _mm256_max_pd(v6, p); v6 = _mm256_blend_pd(lo, hi, 0xA);
        p = _mm256_permute_pd(v7, 0x5); lo = _mm256_min_pd(v7, p); hi = _mm256_max_pd(v7, p); v7 = _mm256_blend_pd(lo, hi, 0xA);
    }
    // Stage 2: 16 distance-2 comparators
    {
        __m256d p, lo, hi;
        p = _mm256_permute4x64_pd(v0, 0x4E); lo = _mm256_min_pd(v0, p); hi = _mm256_max_pd(v0, p); v0 = _mm256_blend_pd(lo, hi, 0xC);
        p = _mm256_permute4x64_pd(v1, 0x4E); lo = _mm256_min_pd(v1, p); hi = _mm256_max_pd(v1, p); v1 = _mm256_blend_pd(lo, hi, 0xC);
        p = _mm256_permute4x64_pd(v2, 0x4E); lo = _mm256_min_pd(v2, p); hi = _mm256_max_pd(v2, p); v2 = _mm256_blend_pd(lo, hi, 0xC);
        p = _mm256_permute4x64_pd(v3, 0x4E); lo = _mm256_min_pd(v3, p); hi = _mm256_max_pd(v3, p); v3 = _mm256_blend_pd(lo, hi, 0xC);
        p = _mm256_permute4x64_pd(v4, 0x4E); lo = _mm256_min_pd(v4, p); hi = _mm256_max_pd(v4, p); v4 = _mm256_blend_pd(lo, hi, 0xC);
        p = _mm256_permute4x64_pd(v5, 0x4E); lo = _mm256_min_pd(v5, p); hi = _mm256_max_pd(v5, p); v5 = _mm256_blend_pd(lo, hi, 0xC);
        p = _mm256_permute4x64_pd(v6, 0x4E); lo = _mm256_min_pd(v6, p); hi = _mm256_max_pd(v6, p); v6 = _mm256_blend_pd(lo, hi, 0xC);
        p = _mm256_permute4x64_pd(v7, 0x4E); lo = _mm256_min_pd(v7, p); hi = _mm256_max_pd(v7, p); v7 = _mm256_blend_pd(lo, hi, 0xC);
    }
    // Stage 3: 16 distance-4 element-wise cross-register
    {
        __m256d lo, hi;
        lo = _mm256_min_pd(v0, v1); hi = _mm256_max_pd(v0, v1); v0 = lo; v1 = hi;
        lo = _mm256_min_pd(v2, v3); hi = _mm256_max_pd(v2, v3); v2 = lo; v3 = hi;
        lo = _mm256_min_pd(v4, v5); hi = _mm256_max_pd(v4, v5); v4 = lo; v5 = hi;
        lo = _mm256_min_pd(v6, v7); hi = _mm256_max_pd(v6, v7); v6 = lo; v7 = hi;
    }
    // Stage 4: 16 distance-8 element-wise cross-register
    {
        __m256d lo, hi;
        lo = _mm256_min_pd(v0, v2); hi = _mm256_max_pd(v0, v2); v0 = lo; v2 = hi;
        lo = _mm256_min_pd(v1, v3); hi = _mm256_max_pd(v1, v3); v1 = lo; v3 = hi;
        lo = _mm256_min_pd(v4, v6); hi = _mm256_max_pd(v4, v6); v4 = lo; v6 = hi;
        lo = _mm256_min_pd(v5, v7); hi = _mm256_max_pd(v5, v7); v5 = lo; v7 = hi;
    }
    // Store to aligned buffer for scalar tail
    alignas(32) double b[32];
    _mm256_store_pd(b,      v0); _mm256_store_pd(b + 4,  v1);
    _mm256_store_pd(b + 8,  v2); _mm256_store_pd(b + 12, v3);
    _mm256_store_pd(b + 16, v4); _mm256_store_pd(b + 20, v5);
    _mm256_store_pd(b + 24, v6); _mm256_store_pd(b + 28, v7);
    // Stages 5+ scalar (121 remaining comparators from sort_net_32 lines 608-628)
    SIMD_SWAP(b[0],b[16]); SIMD_SWAP(b[1],b[8]);  SIMD_SWAP(b[2],b[4]);
    SIMD_SWAP(b[3],b[12]); SIMD_SWAP(b[5],b[10]); SIMD_SWAP(b[6],b[9]);
    SIMD_SWAP(b[7],b[14]); SIMD_SWAP(b[11],b[13]); SIMD_SWAP(b[15],b[31]);
    SIMD_SWAP(b[17],b[24]); SIMD_SWAP(b[18],b[20]); SIMD_SWAP(b[19],b[28]);
    SIMD_SWAP(b[21],b[26]); SIMD_SWAP(b[22],b[25]); SIMD_SWAP(b[23],b[30]);
    SIMD_SWAP(b[27],b[29]);
    SIMD_SWAP(b[1],b[2]);  SIMD_SWAP(b[3],b[5]);  SIMD_SWAP(b[4],b[8]);
    SIMD_SWAP(b[6],b[22]); SIMD_SWAP(b[7],b[11]); SIMD_SWAP(b[9],b[25]);
    SIMD_SWAP(b[10],b[12]); SIMD_SWAP(b[13],b[14]); SIMD_SWAP(b[17],b[18]);
    SIMD_SWAP(b[19],b[21]); SIMD_SWAP(b[20],b[24]); SIMD_SWAP(b[23],b[27]);
    SIMD_SWAP(b[26],b[28]); SIMD_SWAP(b[29],b[30]);
    SIMD_SWAP(b[1],b[17]); SIMD_SWAP(b[2],b[18]); SIMD_SWAP(b[3],b[19]);
    SIMD_SWAP(b[4],b[20]); SIMD_SWAP(b[5],b[10]); SIMD_SWAP(b[7],b[23]);
    SIMD_SWAP(b[8],b[24]); SIMD_SWAP(b[11],b[27]); SIMD_SWAP(b[12],b[28]);
    SIMD_SWAP(b[13],b[29]); SIMD_SWAP(b[14],b[30]); SIMD_SWAP(b[21],b[26]);
    SIMD_SWAP(b[3],b[17]); SIMD_SWAP(b[4],b[16]); SIMD_SWAP(b[5],b[21]);
    SIMD_SWAP(b[6],b[18]); SIMD_SWAP(b[7],b[9]);  SIMD_SWAP(b[8],b[20]);
    SIMD_SWAP(b[10],b[26]); SIMD_SWAP(b[11],b[23]); SIMD_SWAP(b[13],b[25]);
    SIMD_SWAP(b[14],b[28]); SIMD_SWAP(b[15],b[27]); SIMD_SWAP(b[22],b[24]);
    SIMD_SWAP(b[1],b[4]);  SIMD_SWAP(b[3],b[8]);  SIMD_SWAP(b[5],b[16]);
    SIMD_SWAP(b[7],b[17]); SIMD_SWAP(b[9],b[21]); SIMD_SWAP(b[10],b[22]);
    SIMD_SWAP(b[11],b[19]); SIMD_SWAP(b[12],b[20]); SIMD_SWAP(b[14],b[24]);
    SIMD_SWAP(b[15],b[26]); SIMD_SWAP(b[23],b[28]); SIMD_SWAP(b[27],b[30]);
    SIMD_SWAP(b[2],b[5]);  SIMD_SWAP(b[7],b[8]);  SIMD_SWAP(b[9],b[18]);
    SIMD_SWAP(b[11],b[17]); SIMD_SWAP(b[12],b[16]); SIMD_SWAP(b[13],b[22]);
    SIMD_SWAP(b[14],b[20]); SIMD_SWAP(b[15],b[19]); SIMD_SWAP(b[23],b[24]);
    SIMD_SWAP(b[26],b[29]);
    SIMD_SWAP(b[2],b[4]);  SIMD_SWAP(b[6],b[12]); SIMD_SWAP(b[9],b[16]);
    SIMD_SWAP(b[10],b[11]); SIMD_SWAP(b[13],b[17]); SIMD_SWAP(b[14],b[18]);
    SIMD_SWAP(b[15],b[22]); SIMD_SWAP(b[19],b[25]); SIMD_SWAP(b[20],b[21]);
    SIMD_SWAP(b[27],b[29]);
    SIMD_SWAP(b[5],b[6]);  SIMD_SWAP(b[8],b[12]); SIMD_SWAP(b[9],b[10]);
    SIMD_SWAP(b[11],b[13]); SIMD_SWAP(b[14],b[16]); SIMD_SWAP(b[15],b[17]);
    SIMD_SWAP(b[18],b[20]); SIMD_SWAP(b[19],b[23]); SIMD_SWAP(b[21],b[22]);
    SIMD_SWAP(b[25],b[26]);
    SIMD_SWAP(b[3],b[5]);  SIMD_SWAP(b[6],b[7]);  SIMD_SWAP(b[8],b[9]);
    SIMD_SWAP(b[10],b[12]); SIMD_SWAP(b[11],b[14]); SIMD_SWAP(b[13],b[16]);
    SIMD_SWAP(b[15],b[18]); SIMD_SWAP(b[17],b[20]); SIMD_SWAP(b[19],b[21]);
    SIMD_SWAP(b[22],b[23]); SIMD_SWAP(b[24],b[25]); SIMD_SWAP(b[26],b[28]);
    SIMD_SWAP(b[3],b[4]);  SIMD_SWAP(b[5],b[6]);  SIMD_SWAP(b[7],b[8]);
    SIMD_SWAP(b[9],b[10]); SIMD_SWAP(b[11],b[12]); SIMD_SWAP(b[13],b[14]);
    SIMD_SWAP(b[15],b[16]); SIMD_SWAP(b[17],b[18]); SIMD_SWAP(b[19],b[20]);
    SIMD_SWAP(b[21],b[22]); SIMD_SWAP(b[23],b[24]); SIMD_SWAP(b[25],b[26]);
    SIMD_SWAP(b[27],b[28]);
    std::memcpy(x, b, 32 * sizeof(double));
}

// Dispatch: hybrid SIMD sorting network for supported sizes.
ROBSCALE_TARGET_AVX2
static void simd_sort_net_dispatch_avx2(double* x, size_t n) {
    switch (n) {
        case 8:  simd_sort_net_8(x); return;
        case 32: simd_sort_net_32(x); return;
        default: robscale::small_sort(x, n); return;
    }
}

// ============================================================================
// Selection network dispatch
// ============================================================================

ROBSCALE_TARGET_AVX2
static double simd_median_sel_dispatch_avx2(double* x, size_t n) {
    switch (n) {
        case 8:  return simd_median_sel_8(x);
        case 16: return simd_median_sel_16(x);
        case 32: return simd_median_sel_32(x);
        default: return robscale::median_net(x, n);
    }
}

// ============================================================================
// Full-sort dispatch (bitonic sort + extract median)
// ============================================================================

ROBSCALE_TARGET_AVX2
static double simd_median_dispatch_avx2(double* x, size_t n) {
    switch (n) {
        case 4:  return simd_median_4(x);
        case 8:  return simd_median_8(x);
        case 16: return simd_median_16(x);
        case 32: return simd_median_32(x);
        case 64: return simd_median_64(x);
        default: return robscale::median_net(x, n);
    }
}

#undef SIMD_SWAP

} // namespace simd
} // namespace robscale

#endif // ROBSCALE_HAS_AVX2_DISPATCH
#endif // ROBSCALE_SIMD_MEDIAN_H
