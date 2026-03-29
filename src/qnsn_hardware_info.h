#ifndef ROBSCALE_QNSN_HARDWARE_INFO_H
#define ROBSCALE_QNSN_HARDWARE_INFO_H

#include <fstream>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

#if defined(__linux__)
#include <unistd.h>
#elif defined(__APPLE__)
#include <sys/sysctl.h>
#elif defined(_WIN32)
#include <limits>
#include <windows.h>
#endif

// R4: guard <cpuid.h> with GCC/Clang — MSVC uses <intrin.h> and __cpuid()
#if (defined(__x86_64__) || defined(__i386__)) && \
    (defined(__GNUC__) || defined(__clang__))
#include <cpuid.h>
#endif

namespace robscale::qnsn {

static constexpr size_t L2_FALLBACK_KB = 256; // 256 KB

enum class SIMDLevel { None, SSE4_2, AVX2, AVX512, Neon };

struct HardwareInfo {
  size_t l1_cache_size = 32768;   // Default 32KB
  size_t l2_cache_size = 262144;  // Default 256KB
  size_t l2_per_core = 262144;    // Default 256KB — the key threshold input
  size_t l3_cache_size = 4194304; // Default 4MB
  size_t cache_line_size = 64;
  size_t num_physical_cores = 1;
  size_t num_logical_cores = 1;
  SIMDLevel simd_level = SIMDLevel::None;

  HardwareInfo() = default;

  void discover() {
    num_logical_cores = std::max(std::thread::hardware_concurrency(), 1u);

#if defined(__linux__)
    discover_linux();
#elif defined(__APPLE__)
    discover_macos();
#elif defined(_WIN32)
    discover_windows();
#endif

    discover_simd();
  }

private:
  // Reads L2 cache size (per-core) and cache-line size from sysfs
  // (/sys/devices/system/cpu/cpu0/cache/); infers physical core count from
  // the SMT thread_siblings_list.  Falls back to L2_FALLBACK_KB and
  // num_logical_cores=1 on any parse error.
  void discover_linux() {
    try {
      auto read_sysfs = [](const std::string &path) -> size_t {
        std::ifstream f(path);
        if (!f.is_open())
          return 0;
        size_t val;
        f >> val;
        return val;
      };

      // sysfs L2 size is ALREADY per-core
      size_t l2_kb = read_sysfs("/sys/devices/system/cpu/cpu0/cache/index2/size");
      if (l2_kb > 0) {
        l2_per_core = l2_kb * 1024;
        l2_cache_size = l2_per_core;
      }

      size_t cls = read_sysfs(
          "/sys/devices/system/cpu/cpu0/cache/index0/coherency_line_size");
      if (cls > 0)
        cache_line_size = cls;

#ifdef _SC_LEVEL2_CACHE_SIZE
      if (l2_kb == 0) {
        long l2_sc = sysconf(_SC_LEVEL2_CACHE_SIZE);
        if (l2_sc > 0) {
          l2_cache_size = (size_t)l2_sc;
          // S10: sysconf may return total L2 on some kernels; divide as conservative estimate
          l2_per_core = l2_cache_size / std::max(num_logical_cores, size_t(1));
        }
      }
#endif

      // Detect SMT factor from thread_siblings_list to get physical core count
      size_t smt_factor = 1;
      std::ifstream siblings("/sys/devices/system/cpu/cpu0/topology/thread_siblings_list");
      if (siblings.is_open()) {
        std::string content;
        std::getline(siblings, content);
        size_t count = 0;
        std::istringstream ss(content);
        std::string token;
        while (std::getline(ss, token, ',')) {
          auto dash = token.find('-');
          if (dash != std::string::npos)
            count += std::stoul(token.substr(dash + 1)) - std::stoul(token.substr(0, dash)) + 1;
          else
            count += 1;
        }
        if (count > 0) smt_factor = count;
      }
      num_physical_cores = std::max(num_logical_cores / smt_factor, size_t(1));

      // Guard: any l2_per_core below 64KB is physically implausible for L2
      // and indicates a detection error. Fall back to the default (256KB).
      if (l2_per_core < 65536) {
        l2_per_core = L2_FALLBACK_KB * 1024;
        l2_cache_size = l2_per_core;
      }
    } catch (...) {
      l2_per_core = L2_FALLBACK_KB * 1024;
      l2_cache_size = l2_per_core;
      // R5: preserve num_logical_cores from hardware_concurrency() (line 43)
    }
  }

  // Reads total L2 cache size and physical core count via sysctl
  // (hw.l2cachesize, hw.physicalcpu); derives l2_per_core as their ratio.
  // Falls back to the HardwareInfo default (256 KB) if sysctl is unavailable.
  void discover_macos() {
#ifdef __APPLE__
    size_t len;
    len = sizeof(l2_cache_size);
    if (sysctlbyname("hw.l2cachesize", &l2_cache_size, &len, NULL, 0) != 0)
      l2_cache_size = L2_FALLBACK_KB * 1024;
    len = sizeof(cache_line_size);
    sysctlbyname("hw.cachelinesize", &cache_line_size, &len, NULL, 0);
    len = sizeof(num_physical_cores);
    sysctlbyname("hw.physicalcpu", &num_physical_cores, &len, NULL, 0);

    // hw.l2cachesize returns total L2; derive per-core
    if (num_physical_cores > 0)
      l2_per_core = l2_cache_size / num_physical_cores;
    else
      l2_per_core = l2_cache_size;

    // Guard: any l2_per_core below 64KB is physically implausible for L2
    if (l2_per_core < 65536) {
      l2_per_core = L2_FALLBACK_KB * 1024;
      l2_cache_size = l2_per_core;
    }
#endif
  }

  // Reads L2 cache size and physical core count via
  // GetLogicalProcessorInformation(); derives l2_per_core as their ratio.
  // Falls back to the HardwareInfo default (256 KB) if the API call fails.
  void discover_windows() {
#ifdef _WIN32
    DWORD len = 0;
    GetLogicalProcessorInformation(nullptr, &len);
    if (len > 0) {
      std::vector<SYSTEM_LOGICAL_PROCESSOR_INFORMATION>
        buf(len / sizeof(SYSTEM_LOGICAL_PROCESSOR_INFORMATION));
      if (GetLogicalProcessorInformation(buf.data(), &len)) {
        size_t phys = 0;
        for (size_t i = 0; i < buf.size(); ++i) {
          if (buf[i].Relationship == RelationCache &&
              buf[i].Cache.Level == 2) {
            l2_cache_size = buf[i].Cache.Size;
          }
          if (buf[i].Relationship == RelationProcessorCore) {
            phys++;
          }
        }
        if (phys > 0) num_physical_cores = phys;
        // Each RelationCache entry covers one sharing unit (typically per-core
        // L2), so Cache.Size is already the per-core value.  Do not divide by
        // num_physical_cores — that would give a fraction of the per-core size.
        l2_per_core = l2_cache_size;
      }
    }
    // Guard: any l2_per_core below 64KB is physically implausible for L2
    if (l2_per_core < 65536) {
      l2_per_core = L2_FALLBACK_KB * 1024;
      l2_cache_size = l2_per_core;
    }
#endif
  }

  void discover_simd() {
#if defined(__x86_64__) || defined(_M_X64)
#if defined(__GNUC__) || defined(__clang__)
    if (__builtin_cpu_supports("avx512f"))
      simd_level = SIMDLevel::AVX512;
    else if (__builtin_cpu_supports("avx2"))
      simd_level = SIMDLevel::AVX2;
    else if (__builtin_cpu_supports("sse4.2"))
      simd_level = SIMDLevel::SSE4_2;
#else
    simd_level = SIMDLevel::None;
#endif
#elif defined(__aarch64__)
    simd_level = SIMDLevel::Neon;
#endif
  }
};

} // namespace robscale::qnsn

#endif // ROBSCALE_QNSN_HARDWARE_INFO_H
