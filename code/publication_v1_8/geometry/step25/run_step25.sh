#!/usr/bin/env bash
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE"

PYTHON_BIN="${PYTHON_BIN:-python}"
OUTDIR="${STEP25_OUTDIR:-/root/nmda2/step_25/results_step_25_corrected_strong_geometry}"
LOCAL_N="${STEP25_LOCAL_N:-1000}"
THREADS="${STEP25_THREADS:-2}"
SEED="${STEP25_SEED:-20260913}"

export NUMBA_NUM_THREADS="$THREADS"
export OMP_NUM_THREADS="$THREADS"
export OPENBLAS_NUM_THREADS="$THREADS"
export MKL_NUM_THREADS="$THREADS"

EXTRA=()
if [[ "${STEP25_FORCE:-0}" == "1" ]]; then
  EXTRA+=(--force)
fi
if [[ "${STEP25_ALLOW_REPLAY_MISMATCH:-0}" == "1" ]]; then
  EXTRA+=(--allow-replay-mismatch)
fi

"$PYTHON_BIN" step25_corrected_strong_geometry.py \
  --step24 auto \
  --step10 auto \
  --output "$OUTDIR" \
  --local-n "$LOCAL_N" \
  --threads "$THREADS" \
  --seed "$SEED" \
  "${EXTRA[@]}"

RESULT_ZIP="${OUTDIR}.zip"
rm -f "$RESULT_ZIP"
(
  cd "$(dirname "$OUTDIR")"
  zip -qr "$RESULT_ZIP" "$(basename "$OUTDIR")"
)
sha256sum "$RESULT_ZIP" > "${RESULT_ZIP}.sha256"
printf 'Results: %s\n' "$OUTDIR"
printf 'Archive: %s\n' "$RESULT_ZIP"
printf 'SHA256: '
cat "${RESULT_ZIP}.sha256"
