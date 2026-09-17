#!/usr/bin/env bash
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE"

THREADS="${STEP29_THREADS:-2}"
LOCAL_N="${STEP29_LOCAL_N:-1000}"
ARGS=(--threads "$THREADS" --local-n "$LOCAL_N")

if [[ "${STEP29_FORCE:-0}" == "1" ]]; then
  ARGS+=(--force)
fi
if [[ "${STEP29_RESUME:-0}" == "1" ]]; then
  ARGS+=(--resume)
fi
if [[ "${STEP29_SKIP_LOCAL:-0}" == "1" ]]; then
  ARGS+=(--skip-local)
fi

python step29_corrected_publication_consolidation.py "${ARGS[@]}"
