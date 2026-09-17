#!/usr/bin/env bash
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE"

PYTHON_BIN="${PYTHON_BIN:-python}"
OUTDIR="${STEP28_OUTDIR:-/root/nmda2/step_28/results_step_28_residual_connectivity_stress}"
THREADS="${STEP28_THREADS:-2}"
SEARCH_POINTS="${STEP28_SEARCH_POINTS:-17}"
VERIFY_POINTS="${STEP28_VERIFY_POINTS:-61}"
POP_SIZE="${STEP28_POP_SIZE:-32}"
GENERATIONS="${STEP28_GENERATIONS:-20}"
RESTARTS="${STEP28_RESTARTS:-2}"
RESPONSE_TOP="${STEP28_RESPONSE_TOP:-10}"
PAIR_ATTEMPTS="${STEP28_PAIR_ATTEMPTS:-6}"
DEEP_PAIR_ATTEMPTS="${STEP28_DEEP_PAIR_ATTEMPTS:-2}"
MAX_COMPONENT_PAIRS="${STEP28_MAX_COMPONENT_PAIRS:-6}"
MAX_WAYPOINTS="${STEP28_MAX_WAYPOINTS:-3}"
INITIAL_SIGMA="${STEP28_INITIAL_SIGMA:-0.9}"
VERIFY_TRIGGER="${STEP28_VERIFY_TRIGGER:-1.12}"
FLEX_TRIGGER="${STEP28_FLEX_TRIGGER:-1.50}"
SEED="${STEP28_SEED:-20260914}"

export NUMBA_NUM_THREADS="$THREADS"
export OMP_NUM_THREADS="$THREADS"
export OPENBLAS_NUM_THREADS="$THREADS"
export MKL_NUM_THREADS="$THREADS"

EXTRA=()
if [[ "${STEP28_FORCE:-0}" == "1" ]]; then
  EXTRA+=(--force)
elif [[ "${STEP28_RESUME:-0}" == "1" ]]; then
  EXTRA+=(--resume)
fi
if [[ "${STEP28_DISABLE_FLEX_AUX:-0}" == "1" ]]; then
  EXTRA+=(--disable-flex-aux)
fi
if [[ "${STEP28_ALLOW_ACCEPTED_BOUNDS_FALLBACK:-0}" == "1" ]]; then
  EXTRA+=(--allow-accepted-bounds-fallback)
fi

"$PYTHON_BIN" step28_residual_connectivity_stress.py \
  --step25 auto \
  --step27 auto \
  --output "$OUTDIR" \
  --threads "$THREADS" \
  --search-points "$SEARCH_POINTS" \
  --verify-points "$VERIFY_POINTS" \
  --population "$POP_SIZE" \
  --generations "$GENERATIONS" \
  --restarts "$RESTARTS" \
  --response-top "$RESPONSE_TOP" \
  --pair-attempts "$PAIR_ATTEMPTS" \
  --deep-pair-attempts "$DEEP_PAIR_ATTEMPTS" \
  --max-component-pairs "$MAX_COMPONENT_PAIRS" \
  --max-waypoints "$MAX_WAYPOINTS" \
  --initial-sigma "$INITIAL_SIGMA" \
  --verify-trigger "$VERIFY_TRIGGER" \
  --flex-trigger "$FLEX_TRIGGER" \
  --seed "$SEED" \
  "${EXTRA[@]}"

RESULT_ZIP="${OUTDIR}.zip"
rm -f "$RESULT_ZIP" "${RESULT_ZIP}.sha256"
(
  cd "$(dirname "$OUTDIR")"
  zip -qr "$RESULT_ZIP" "$(basename "$OUTDIR")"
)
sha256sum "$RESULT_ZIP" > "${RESULT_ZIP}.sha256"
printf 'Results: %s\n' "$OUTDIR"
printf 'Archive: %s\n' "$RESULT_ZIP"
printf 'SHA256: '
cat "${RESULT_ZIP}.sha256"
