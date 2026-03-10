# Findings & Decisions

## Requirements
- Convert `README.md` to a Quarto notebook (`README.qmd`).
- Dynamically compute all benchmarks, statistics, figures, and tables.
- Compare `robscale` (ROBSCALE_FAST=1), `robscale` (ROBSCALE_FAST=0), `revss`, and `robustbase`.
- Ensure the generated `README.md` is "honest" (no hardcoded numbers, fully reproducible).
- Can be used to generate a new README before each release.
- Provide a comprehensive plan executable by someone else, explaining *why* things are done.

## Research Findings
- The current `README.md` has 589 lines of detailed mathematical and performance explanations.
- It contains hardcoded claims like "**11–39×** speedups", table data, and a statically linked image `benchmarks/cran_mode/full_speedup_panel.png`.
- The package exposes a C++ backend whose optimization level depends on compile-time environment variables (`ROBSCALE_FAST=1`).
- The project does not currently use a `_targets` pipeline or `renv` for the README generation.

## Technical Decisions
| Decision | Rationale |
|----------|-----------|
| Decouple Compute / Render | Using `targets` prevents massive latency when re-rendering the README simply to fix a typo. The benchmarks will be run occasionally; the rendering can happen frequently. |
| Temporary R Libraries | To accurately benchmark unoptimized vs optimized, the target workers must install the package into a temporary `.libPaths()` isolation to prevent state pollution. |
| `format: gfm` | Quarto must output `gfm` (GitHub Flavored Markdown) to seamlessly integrate back into the standard GitHub interface and CRAN checks. |
| `planning-with-files` | Storing these Manus files locally helps keep high-level goals in focus during long execution. |

## Issues Encountered
| Issue | Resolution |
|-------|------------|
|       |            |

## Visual/Browser Findings
- N/A

---
*Update this file after every 2 view/browser/search operations*
*This prevents visual information from being lost*
