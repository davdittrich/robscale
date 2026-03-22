#!/usr/bin/env bash
# Run a gate check script with CPU governor=performance to collapse bimodal
# frequency-scaling noise at n<=16 (~2µs operations).
#
# Optionally applies FIFO-99 real-time scheduling if sudo chrt is available
# without a password prompt. Falls back to nice -n -20.
#
# Usage: bash bench/run_gate.sh [gate_script.R]
#   Default: bench/iqr_gate_check.R
#
# Requires: cpupower NOPASSWD in sudoers (already configured on this machine).

set -euo pipefail

SCRIPT="${1:-bench/iqr_gate_check.R}"

# ---- 1. Set performance governor ----
GOV_FILE="/sys/devices/system/cpu/cpu0/cpufreq/scaling_governor"
PREV_GOV=""

if [[ -f "$GOV_FILE" ]]; then
  PREV_GOV=$(cat "$GOV_FILE")
  if [[ "$PREV_GOV" == "performance" ]]; then
    echo "CPU governor: already performance"
    PREV_GOV=""  # nothing to restore
  else
    echo "CPU governor: $PREV_GOV → performance"
    if ! sudo -n cpupower frequency-set -g performance >/dev/null 2>&1; then
      echo "WARNING: cpupower failed — governor unchanged. Results at n<=16 may be noisy." >&2
      PREV_GOV=""
    fi
  fi
fi

cleanup() {
  if [[ -n "$PREV_GOV" ]]; then
    echo "CPU governor: performance → $PREV_GOV"
    sudo -n cpupower frequency-set -g "$PREV_GOV" >/dev/null 2>&1 || true
  fi
}
trap cleanup EXIT

# ---- 2. Set scheduling priority for this shell (becomes child's priority) ----
# Try to set FIFO-99 on this PID via sudo (without exec, so trap fires).
# If unavailable, use nice -20 as fallback.
if sudo -n chrt -f -p 99 $$ 2>/dev/null; then
  echo "(scheduling: FIFO-99 on PID $$)"
else
  echo "(scheduling: FIFO unavailable; using nice -n -20)"
  renice -n -20 $$ 2>/dev/null || true
fi

echo "Running $SCRIPT..."
echo
Rscript "$SCRIPT"
