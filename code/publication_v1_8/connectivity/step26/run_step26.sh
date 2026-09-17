#!/usr/bin/env bash
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE"

PYTHON_BIN="${PYTHON_BIN:-python}"
OUTDIR="${STEP26_OUTDIR:-/root/nmda2/step_26/results_step_26_strong_connectivity}"
THREADS="${STEP26_THREADS:-2}"
COARSE_POINTS="${STEP26_COARSE_POINTS:-21}"
REFINE_POINTS="${STEP26_REFINE_POINTS:-41}"
MAX_EDGES="${STEP26_MAX_EDGES:-180}"
KNN_K="${STEP26_KNN_K:-4}"

export NUMBA_NUM_THREADS="$THREADS"
export OMP_NUM_THREADS="$THREADS"
export OPENBLAS_NUM_THREADS="$THREADS"
export MKL_NUM_THREADS="$THREADS"

EXTRA=()
if [[ "${STEP26_FORCE:-0}" == "1" ]]; then
  EXTRA+=(--force)
fi
if [[ "${STEP26_ALLOW_STEP25_STATUS:-0}" == "1" ]]; then
  EXTRA+=(--allow-step25-status)
fi

"$PYTHON_BIN" step26_strong_connectivity.py \
  --step25 auto \
  --output "$OUTDIR" \
  --threads "$THREADS" \
  --coarse-points "$COARSE_POINTS" \
  --refine-points "$REFINE_POINTS" \
  --max-edges "$MAX_EDGES" \
  --knn-k "$KNN_K" \
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
