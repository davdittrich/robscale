#ifndef ROBSCALE_SORT_ALTERNATIVES_H
#define ROBSCALE_SORT_ALTERNATIVES_H

#include <algorithm>
#include <cstring>
#include "sort_net.h"

// Platform detection for SIMD
#if defined(__aarch64__) || defined(_M_ARM64)
  #define ROBSCALE_HAS_NEON 1
  #include <arm_neon.h>
#elif defined(__x86_64__) || defined(_M_X64)
  #if defined(__AVX2__)
    #define ROBSCALE_HAS_AVX2 1
    #include <immintrin.h>
  #elif defined(__SSE2__)
    #define ROBSCALE_HAS_SSE2 1
    #include <emmintrin.h>
  #endif
#endif

namespace robscale {
namespace sort_alt {

// =========================================================================
// 1. Knuth sorting networks (delegate to existing sort_net.h)
// =========================================================================
template <typename T>
inline void knuth_network(T* x, size_t n) {
  small_sort(x, n);
}

// =========================================================================
// 2. std::sort wrapper
// =========================================================================
template <typename T>
inline void std_sort(T* x, size_t n) {
  std::sort(x, x + n);
}

// =========================================================================
// 3. Insertion sort with sentinel
// =========================================================================
template <typename T>
inline void insertion_sort(T* x, size_t n) {
  if (n <= 1) return;
  // Find minimum and place it at x[0] as sentinel
  size_t min_idx = 0;
  for (size_t i = 1; i < n; ++i) {
    if (x[i] < x[min_idx]) min_idx = i;
  }
  if (min_idx != 0) std::swap(x[0], x[min_idx]);
  // Sentinel-guarded insertion sort (no bounds check needed)
  for (size_t i = 2; i < n; ++i) {
    T key = x[i];
    size_t j = i;
    while (x[j - 1] > key) {
      x[j] = x[j - 1];
      --j;
    }
    x[j] = key;
  }
}

// =========================================================================
// 4. SIMD register sort
// =========================================================================

// --- NEON helpers (ARM64) ------------------------------------------------
#if ROBSCALE_HAS_NEON

// Sort 2 doubles via NEON
inline void simd_sort_neon_2(double* x) {
  float64x2_t a = vld1q_f64(x);
  // a = {x[0], x[1]}; extract and compare
  double v0 = vgetq_lane_f64(a, 0);
  double v1 = vgetq_lane_f64(a, 1);
  if (v0 > v1) {
    vst1q_f64(x, vcombine_f64(vcreate_f64(0), vcreate_f64(0)));
    x[0] = v1;
    x[1] = v0;
  }
}

// Sort 4 doubles via NEON using sorting network
inline void simd_sort_neon_4(double* x) {
  // Load pairs
  float64x2_t ab = vld1q_f64(x);      // {x[0], x[1]}
  float64x2_t cd = vld1q_f64(x + 2);  // {x[2], x[3]}

  // Stage 1: compare-swap (0,1) and (2,3) — 2 comparators in parallel
  float64x2_t mn1 = vminq_f64(ab, vextq_f64(ab, ab, 1));
  float64x2_t mx1 = vmaxq_f64(ab, vextq_f64(ab, ab, 1));
  double a0 = vgetq_lane_f64(mn1, 0);
  double a1 = vgetq_lane_f64(mx1, 0);

  float64x2_t mn2 = vminq_f64(cd, vextq_f64(cd, cd, 1));
  float64x2_t mx2 = vmaxq_f64(cd, vextq_f64(cd, cd, 1));
  double a2 = vgetq_lane_f64(mn2, 0);
  double a3 = vgetq_lane_f64(mx2, 0);

  // Stage 2: compare-swap (0,2) and (1,3)
  double b0 = std::min(a0, a2);
  double b2 = std::max(a0, a2);
  double b1 = std::min(a1, a3);
  double b3 = std::max(a1, a3);

  // Stage 3: compare-swap (1,2)
  x[0] = b0;
  x[1] = std::min(b1, b2);
  x[2] = std::max(b1, b2);
  x[3] = b3;
}

#endif // ROBSCALE_HAS_NEON

// --- AVX2 helpers (x86_64) -----------------------------------------------
#if ROBSCALE_HAS_AVX2

// Sort 4 doubles via AVX2 sorting network
inline void simd_sort_avx2_4(double* x) {
  __m256d v = _mm256_loadu_pd(x);

  // Stage 1: compare-swap (0,1) and (2,3)
  __m256d shuf1 = _mm256_permute4x64_pd(v, 0b10110001); // swap within pairs: {1,0,3,2}
  __m256d lo1 = _mm256_min_pd(v, shuf1);
  __m256d hi1 = _mm256_max_pd(v, shuf1);
  // Take min from even positions, max from odd: {min(0,1), max(0,1), min(2,3), max(2,3)}
  v = _mm256_blend_pd(lo1, hi1, 0b1010);

  // Stage 2: compare-swap (0,2) and (1,3)
  __m256d shuf2 = _mm256_permute4x64_pd(v, 0b01001110); // swap pairs: {2,3,0,1}
  __m256d lo2 = _mm256_min_pd(v, shuf2);
  __m256d hi2 = _mm256_max_pd(v, shuf2);
  v = _mm256_blend_pd(lo2, hi2, 0b1100);

  // Stage 3: compare-swap (1,2)
  __m256d shuf3 = _mm256_permute4x64_pd(v, 0b11011000); // {0,2,1,3}
  __m256d lo3 = _mm256_min_pd(v, shuf3);
  __m256d hi3 = _mm256_max_pd(v, shuf3);
  v = _mm256_blend_pd(lo3, hi3, 0b0110);

  _mm256_storeu_pd(x, v);
}

// Sort 8 doubles via AVX2 bitonic merge network
inline void simd_sort_avx2_8(double* x) {
  __m256d a = _mm256_loadu_pd(x);
  __m256d b = _mm256_loadu_pd(x + 4);

  // --- Sort each register of 4 (using the n=4 network inline) ---
  // Sort 'a'
  {
    __m256d s1 = _mm256_permute4x64_pd(a, 0b10110001);
    __m256d lo = _mm256_min_pd(a, s1);
    __m256d hi = _mm256_max_pd(a, s1);
    a = _mm256_blend_pd(lo, hi, 0b1010);

    __m256d s2 = _mm256_permute4x64_pd(a, 0b01001110);
    lo = _mm256_min_pd(a, s2);
    hi = _mm256_max_pd(a, s2);
    a = _mm256_blend_pd(lo, hi, 0b1100);

    __m256d s3 = _mm256_permute4x64_pd(a, 0b11011000);
    lo = _mm256_min_pd(a, s3);
    hi = _mm256_max_pd(a, s3);
    a = _mm256_blend_pd(lo, hi, 0b0110);
  }
  // Sort 'b'
  {
    __m256d s1 = _mm256_permute4x64_pd(b, 0b10110001);
    __m256d lo = _mm256_min_pd(b, s1);
    __m256d hi = _mm256_max_pd(b, s1);
    b = _mm256_blend_pd(lo, hi, 0b1010);

    __m256d s2 = _mm256_permute4x64_pd(b, 0b01001110);
    lo = _mm256_min_pd(b, s2);
    hi = _mm256_max_pd(b, s2);
    b = _mm256_blend_pd(lo, hi, 0b1100);

    __m256d s3 = _mm256_permute4x64_pd(b, 0b11011000);
    lo = _mm256_min_pd(b, s3);
    hi = _mm256_max_pd(b, s3);
    b = _mm256_blend_pd(lo, hi, 0b0110);
  }

  // --- Bitonic merge: reverse b, then merge ---
  b = _mm256_permute4x64_pd(b, 0b00011011); // reverse: {3,2,1,0}

  // Compare-swap across registers
  __m256d lo = _mm256_min_pd(a, b);
  __m256d hi = _mm256_max_pd(a, b);
  a = lo;
  b = hi;

  // Merge within each register (2 stages)
  // Stage 1: swap halves and compare
  {
    __m256d sa = _mm256_permute4x64_pd(a, 0b01001110);
    __m256d la = _mm256_min_pd(a, sa);
    __m256d ha = _mm256_max_pd(a, sa);
    a = _mm256_blend_pd(la, ha, 0b1100);

    __m256d sb = _mm256_permute4x64_pd(b, 0b01001110);
    __m256d lb = _mm256_min_pd(b, sb);
    __m256d hb = _mm256_max_pd(b, sb);
    b = _mm256_blend_pd(lb, hb, 0b1100);
  }
  // Stage 2: swap adjacent and compare
  {
    __m256d sa = _mm256_permute4x64_pd(a, 0b10110001);
    __m256d la = _mm256_min_pd(a, sa);
    __m256d ha = _mm256_max_pd(a, sa);
    a = _mm256_blend_pd(la, ha, 0b1010);

    __m256d sb = _mm256_permute4x64_pd(b, 0b10110001);
    __m256d lb = _mm256_min_pd(b, sb);
    __m256d hb = _mm256_max_pd(b, sb);
    b = _mm256_blend_pd(lb, hb, 0b1010);
  }

  _mm256_storeu_pd(x, a);
  _mm256_storeu_pd(x + 4, b);
}

#endif // ROBSCALE_HAS_AVX2

// --- SIMD dispatcher -----------------------------------------------------
template <typename T>
inline void simd_sort(T* x, size_t n) {
  // SIMD paths only for double
  if constexpr (std::is_same_v<T, double>) {
#if ROBSCALE_HAS_NEON
    switch (n) {
      case 0: case 1: return;
      case 2: simd_sort_neon_2(x); return;
      case 4: simd_sort_neon_4(x); return;
      default: break; // fall through to scalar network
    }
#elif ROBSCALE_HAS_AVX2
    switch (n) {
      case 0: case 1: return;
      case 4: simd_sort_avx2_4(x); return;
      case 8: simd_sort_avx2_8(x); return;
      default: break; // fall through to scalar network
    }
#endif
  }
  // Fallback: Knuth network for sizes without a SIMD path
  small_sort(x, n);
}

} // namespace sort_alt
} // namespace robscale

#endif // ROBSCALE_SORT_ALTERNATIVES_H
