# README Conversion and Automation Plan

## Goal
Convert the static `README.md` into a reproducible Quarto notebook (`README.qmd`) that dynamically computes all benchmarks, figures, and tables, comparing optimized vs unoptimized `robscale` against legacy alternatives, ensuring 100% honesty in performance claims.

## Current Phase
Phase 1: Planning

## Phases

### Phase 1: Planning & Structure
- [x] Analyze current README structure and hardcoded claims
- [x] Define pipeline architecture (targets + Quarto)
- [x] Document the comprehensive plan
- **Status:** complete

### Phase 2: Pipeline & Benchmarking Infrastructure
- [ ] Initialize `renv` to lock down dependencies (quarto, targets, bench, etc.)
- [ ] Create `R/benchmarks.R` with robust timing functions
- [ ] Setup a target that compiles `robscale` with `ROBSCALE_FAST=0` and benchmarks it
- [ ] Setup a target that compiles `robscale` with `ROBSCALE_FAST=1` and benchmarks it
- [ ] Save benchmark results and session info (CPU, BLAS, etc.) to target store
- **Status:** pending

### Phase 3: Quarto Conversion
- [ ] Rename `README.md` to `README.qmd`
- [ ] Configure YAML for `format: gfm`
- [ ] Replace hardcoded metric statements (e.g., "11-39x speedups") with inline R code (`tar_read()`)
- [ ] Replace markdown tables with dynamic generation (e.g., `knitr::kable()`)
- [ ] Generate Figure 1 dynamically from benchmark data
- **Status:** pending

### Phase 4: Release Automation & Verification
- [ ] Create `scripts/update-readme.sh` to run the pipeline and render
- [ ] Update `.Rbuildignore` to exclude pipeline files
- [ ] Validate final `README.md` output against expected layout
- **Status:** pending

## Key Questions
1. How exactly will we cleanly install the package twice with different `Sys.setenv` parameters in the `targets` pipeline without polluting the user's main R library?
   *Decision: We will use a temporary library path (`.libPaths(temp_lib)`) during the benchmarking targets for the active installations.*

## Decisions Made
| Decision | Rationale |
|----------|-----------|
| Use `targets` for benchmarks | Benchmarking takes time. We must separate computation from rendering so tweaking README prose doesn't trigger a 20-minute benchmark run. |
| Use Quarto + inline R | Guarantees "honesty". Text claims like "10x faster" will be programmatically linked to the actual data, eliminating drift. |
| Two-pass compilation | Running both `ROBSCALE_FAST=0` and `1` requires compiling the package from source twice during the benchmark phase to ensure apples-to-apples comparison on the same machine. |

## Errors Encountered
| Error | Attempt | Resolution |
|-------|---------|------------|
|       | 1       |            |

## Notes
- Update phase status as you progress: pending → in_progress → complete
- Re-read this plan before major decisions (attention manipulation)
- Log ALL errors - they help avoid repetition
