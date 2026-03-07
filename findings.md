# Findings & Decisions

## Requirements

- Integrate SLEEF for vectorized `tanh`.
- Improve performance of `robscale` package functions.
- Support Linux (AVX2, AVX512) and macOS (NEON).
- Maintain existing fallbacks (Accelerate, OpenMP).
- Follow `clean-code` and `scientific-method` principles.

## Research Findings

- `robscale-art/experiments/run-decomposition-v5.R` shows SLEEF provides significant speedup for large `n`.
- Existing implementation in `revss_temp/src/robust_core.h` uses `bulk_tanh` with Accelerate or OpenMP SIMD.
- SLEEF `tanh` functions: `Sleef_tanhd4_u10avx2` (AVX2), `Sleef_tanhd8_u10avx512f` (AVX512), `Sleef_tanhd2_u10` (NEON).
- CPU SIMD detection is needed in `configure` to set correct compiler flags (`-mavx2 -mfma` etc.).

## Technical Decisions

| Decision | Rationale |
|----------|-----------|
| Use `configure` + `Makevars.in` | Required for portable system dependency detection in R. |
| Prioritize Accelerate over SLEEF | On macOS, Accelerate (vvtanh) is faster and preferred; disabling SLEEF on Darwin. |
| Target SLEEF 1.0 accuracy (`_u10`) | Matches standard `std::tanh` accuracy closely while being fast (Linux). |
| Modularize `bulk_tanh` | Follow `clean-code` by separating SIMD logic into distinct inline functions or blocks. |

## Issues Encountered

| Issue | Resolution |
|-------|------------|
| Canceled tool call in configure | Verified configure state via view_file; it was applied. |

## Resources

- SLEEF Documentation: https://sleef.org/
- `robscale` source: `revss_temp/`
- Experiment script: `robscale-art/experiments/run-decomposition-v5.R`

## Visual/Browser Findings

- Experiment 5 results show `sleef/scalar` ratio improved as `n` increased, with significant gains at `n=100`.
