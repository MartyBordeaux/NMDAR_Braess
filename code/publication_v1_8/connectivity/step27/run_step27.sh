#!/usr/bin/env bash
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE"

PYTHON_BIN="${PYTHON_BIN:-python}"
OUTDIR="${STEP27_OUTDIR:-/root/nmda2/step_27/results_step_27_curved_bridge_optimization}"
THREADS="${STEP27_THREADS:-2}"
SEARCH_POINTS="${STEP27_SEARCH_POINTS:-13}"
VERIFY_POINTS="${STEP27_VERIFY_POINTS:-41}"
N_PCS="${STEP27_N_PCS:-5}"
N_MODES="${STEP27_N_MODES:-2}"
POP_SIZE="${STEP27_POP_SIZE:-24}"
GENERATIONS="${STEP27_GENERATIONS:-18}"
JOINT_GENERATIONS="${STEP27_JOINT_GENERATIONS:-10}"
MAX_COMPONENT_PAIRS="${STEP27_MAX_COMPONENT_PAIRS:-12}"
PAIR_ATTEMPTS="${STEP27_PAIR_ATTEMPTS:-3}"
SEED="${STEP27_SEED:-20260914}"

export NUMBA_NUM_THREADS="$THREADS"
export OMP_NUM_THREADS="$THREADS"
export OPENBLAS_NUM_THREADS="$THREADS"
export MKL_NUM_THREADS="$THREADS"

EXTRA=()
if [[ "${STEP27_FORCE:-0}" == "1" ]]; then
  EXTRA+=(--force)
elif [[ "${STEP27_RESUME:-0}" == "1" ]]; then
  EXTRA+=(--resume)
fi
if [[ "${STEP27_ALLOW_ACCEPTED_BOUNDS_FALLBACK:-0}" == "1" ]]; then
  EXTRA+=(--allow-accepted-bounds-fallback)
fi
if [[ "${STEP27_DISABLE_FLEX_AUX:-0}" == "1" ]]; then
  EXTRA+=(--disable-flex-aux)
fi

"$PYTHON_BIN" step27_curved_bridge_optimization.py \
  --step25 auto \
  --step26 auto \
  --output "$OUTDIR" \
  --threads "$THREADS" \
  --search-points "$SEARCH_POINTS" \
  --verify-points "$VERIFY_POINTS" \
  --n-pcs "$N_PCS" \
  --n-modes "$N_MODES" \
  --population "$POP_SIZE" \
  --generations "$GENERATIONS" \
  --joint-generations "$JOINT_GENERATIONS" \
  --max-component-pairs "$MAX_COMPONENT_PAIRS" \
  --pair-attempts "$PAIR_ATTEMPTS" \
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
