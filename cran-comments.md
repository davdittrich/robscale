# CRAN submission comments for robscale 0.2.2

## Changes in 0.2.2

This patch fixes ERROR-level check failures on r-devel-linux-x86_64-fedora-clang,
r-devel-linux-x86_64-fedora-gcc, and r-oldrel-windows-x86_64 caused by
`revss` v3.0.0 changing the defaults of `robScale()` and `robLoc()` to use
bias-corrected "AA" constants.

- Tests: Skip `robScale()` and `robLoc()` cross-validation against `revss`
  when `revss >= 3.0.0` (breaking upstream change). The `adm()` cross-check
  remains active for all versions.
- Tests: Added frozen golden reference tests (90 values each) for `robScale()`
  and `robLoc()`, verified against `revss` v2.0.0. These catch regressions in
  our own algorithm independent of upstream version changes.
- No changes to the robscale algorithm or any source code.

## Test environments

- macOS (ARM64), R 4.5.3, Apple Clang
- Fedora Linux (x86_64), R-devel, GCC 15.2.1
- Fedora Linux (x86_64), R-devel, GCC 15.2.1 (Intel)

## R CMD check results

0 errors | 0 warnings | 1 note

- **GNU make**: `SystemRequirements` field lists GNU make (needed for
  `$(shell)` in Makevars). Already declared in DESCRIPTION.
