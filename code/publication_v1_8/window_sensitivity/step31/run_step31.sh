#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
THREADS="${STEP31_THREADS:-2}"
RAW_ROOT="${STEP31_RAW_ROOT:-auto}"
EXTRA=()
[[ "${STEP31_RESUME:-0}" == "1" ]] && EXTRA+=(--resume)
[[ "${STEP31_FORCE:-0}" == "1" ]] && EXTRA+=(--force)
[[ "${STEP31_GEOMETRY_ONLY:-0}" == "1" ]] && EXTRA+=(--geometry-only)
[[ -n "${STEP31_STEP2:-}" ]] && EXTRA+=(--step2 "$STEP31_STEP2")
[[ -n "${STEP31_STEP10:-}" ]] && EXTRA+=(--step10 "$STEP31_STEP10")
[[ -n "${STEP31_STEP24:-}" ]] && EXTRA+=(--step24 "$STEP31_STEP24")
[[ -n "${STEP31_STEP29:-}" ]] && EXTRA+=(--step29 "$STEP31_STEP29")
export NUMBA_NUM_THREADS="$THREADS"
python "$HERE/step31_window_robustness_final.py" --threads "$THREADS" --raw-root "$RAW_ROOT" "${EXTRA[@]}" "$@"
OUT="/root/nmda2/step_31/results_step_31_window_robustness_final"
if [[ -d "$OUT" ]]; then
  (cd "$(dirname "$OUT")" && zip -qr "results_step_31_window_robustness_final.zip" "$(basename "$OUT")")
fi
