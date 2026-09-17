#!/usr/bin/env bash
set -euo pipefail

REPO_URL="https://github.com/MartyBordeaux/NMDAR_Braess.git"
BRANCH="publication-v1.8-staging"
WORKDIR="${WORKDIR:-/root/NMDAR_Braess_publication_v1_8}"
RAW_ROOT="${RAW_ROOT:-/root/nmda/IV_NMDA}"
NMDA2_ROOT="${NMDA2_ROOT:-/root/nmda2}"
NMDA_ROOT="${NMDA_ROOT:-/root/nmda}"

command -v git >/dev/null || { echo "git is required" >&2; exit 1; }
command -v git-lfs >/dev/null || { echo "git-lfs is required (install git-lfs first)" >&2; exit 1; }

if [[ ! -d "$RAW_ROOT" ]]; then
  echo "Raw ABF root not found: $RAW_ROOT" >&2
  exit 1
fi

rm -rf "$WORKDIR"
git clone "$REPO_URL" "$WORKDIR"
cd "$WORKDIR"
git checkout "$BRANCH"
git lfs install
mkdir -p data/raw code/publication_v1_8/{calibration,model,connectivity,consolidation,window_sensitivity} data/source_data_v1_8

# -----------------------------------------------------------------------------
# Raw electrophysiology
# -----------------------------------------------------------------------------
rm -rf data/raw/IV_NMDA
cp -a "$RAW_ROOT" data/raw/IV_NMDA
find data/raw/IV_NMDA -type f -iname '*.abf' -print0 | sort -z | xargs -0 sha256sum > data/raw/SHA256SUMS.txt
find data/raw/IV_NMDA -type f -iname '*.abf' | sort > data/raw/ABF_FILELIST.txt
printf 'Raw source copied from: %s\nGenerated: %s\n' "$RAW_ROOT" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" > data/raw/PROVENANCE.txt

# -----------------------------------------------------------------------------
# Final computational code used by manuscript v1.8
# -----------------------------------------------------------------------------
copy_pipeline () {
  local src="$1"
  local dst="$2"
  [[ -d "$src" ]] || { echo "Missing pipeline directory: $src" >&2; exit 1; }
  mkdir -p "$dst"
  find "$src" -maxdepth 1 -type f \( -name '*.py' -o -name '*.sh' -o -name 'requirements.txt' -o -name 'README.md' -o -name 'VALIDATION.md' -o -name 'SHA256SUMS.txt' \) -exec cp -f {} "$dst"/ \;
}

# Step 24: corrected calibration-window scoring
copy_pipeline "$NMDA2_ROOT/NMDAR_Braess_step24_calibration_window_correction_v1_0" \
              code/publication_v1_8/calibration/step24

# Step 26 includes the frozen Step25 model engine snapshot used downstream.
copy_pipeline "$NMDA2_ROOT/NMDAR_Braess_step26_strong_connectivity_v1_0" \
              code/publication_v1_8/connectivity/step26
if [[ -f "$NMDA2_ROOT/NMDAR_Braess_step26_strong_connectivity_v1_0/step25_engine_snapshot.py" ]]; then
  cp -f "$NMDA2_ROOT/NMDAR_Braess_step26_strong_connectivity_v1_0/step25_engine_snapshot.py" \
        code/publication_v1_8/model/step25_engine.py
fi
copy_pipeline "$NMDA2_ROOT/NMDAR_Braess_step27_curved_bridge_optimization_v1_0" \
              code/publication_v1_8/connectivity/step27
copy_pipeline "$NMDA2_ROOT/NMDAR_Braess_step28_residual_connectivity_stress_v1_0" \
              code/publication_v1_8/connectivity/step28
copy_pipeline "$NMDA2_ROOT/NMDAR_Braess_step29_corrected_publication_consolidation_v1_0" \
              code/publication_v1_8/consolidation/step29
copy_pipeline "$NMDA2_ROOT/NMDAR_Braess_step31_window_robustness_final_v1_0_2" \
              code/publication_v1_8/window_sensitivity/step31

# Preserve historical candidate generator / Base model implementations already
# used upstream. These locations may already exist in the repository; refresh
# from the server only when source files are available.
if [[ -f "$NMDA_ROOT/10_matched_rerouting.py" ]]; then
  cp -f "$NMDA_ROOT/10_matched_rerouting.py" code/publication_v1_8/model/
fi

# -----------------------------------------------------------------------------
# Final machine-readable source data
# -----------------------------------------------------------------------------
STEP29_RESULTS="$NMDA2_ROOT/step_29/results_step_29_corrected_publication_consolidation"
STEP31_RESULTS="$NMDA2_ROOT/step_31/results_step_31_window_robustness_final"
STEP24_RESULTS="$NMDA2_ROOT/NMDAR_Braess_step24_calibration_window_correction_v1_0/results_step_24_calibration_window_correction"

if [[ -d "$STEP29_RESULTS" ]]; then
  mkdir -p data/source_data_v1_8/step29
  find "$STEP29_RESULTS" -maxdepth 1 -type f \( -name '*.csv' -o -name '*.csv.gz' -o -name '*.json' -o -name '*.md' \) -exec cp -f {} data/source_data_v1_8/step29/ \;
else
  echo "Missing Step29 results: $STEP29_RESULTS" >&2
  exit 1
fi

if [[ -d "$STEP31_RESULTS" ]]; then
  mkdir -p data/source_data_v1_8/step31
  find "$STEP31_RESULTS" -maxdepth 1 -type f \( -name '*.csv' -o -name '*.csv.gz' -o -name '*.json' -o -name '*.md' \) -exec cp -f {} data/source_data_v1_8/step31/ \;
else
  echo "Missing Step31 results: $STEP31_RESULTS" >&2
  exit 1
fi

if [[ -d "$STEP24_RESULTS" ]]; then
  mkdir -p data/source_data_v1_8/step24
  find "$STEP24_RESULTS" -maxdepth 1 -type f \( -name '*.csv' -o -name '*.csv.gz' -o -name '*.json' -o -name '*.md' \) -exec cp -f {} data/source_data_v1_8/step24/ \;
fi

# Frozen full Step10 candidate table, needed to reproduce corrected calibration.
if [[ -f "$NMDA_ROOT/results_step10.zip" ]]; then
  mkdir -p data/model_inputs/step10_full
  cp -f "$NMDA_ROOT/results_step10.zip" data/model_inputs/step10_full/results_step10.zip
fi

# -----------------------------------------------------------------------------
# Audit and publish
# -----------------------------------------------------------------------------
find code/publication_v1_8 -type f -print | sort > code/publication_v1_8/FILELIST.txt
find data/source_data_v1_8 -type f -print | sort > data/source_data_v1_8/FILELIST.txt

ABF_N=$(find data/raw/IV_NMDA -type f -iname '*.abf' | wc -l)
echo "ABF files staged: $ABF_N"
[[ "$ABF_N" -gt 0 ]] || { echo "No ABF files found" >&2; exit 1; }

# Check that the final headline files exist before pushing.
required=(
  data/source_data_v1_8/step29/02_final_strong40_topology.csv
  data/source_data_v1_8/step29/03_final_verified_network_edges.csv
  data/source_data_v1_8/step29/05_publication_representatives.csv
  data/source_data_v1_8/step31/06_final40_window_retention.csv
  data/source_data_v1_8/step31/10_scientific_summary.json
  code/publication_v1_8/model/step25_engine.py
)
for f in "${required[@]}"; do
  [[ -f "$f" ]] || { echo "Required publication file missing: $f" >&2; exit 1; }
done

git add .gitattributes README.md REPRODUCIBILITY.md PUBLICATION_VERSION.md code data tools
if git diff --cached --quiet; then
  echo "Nothing new to commit."
  exit 0
fi

git commit -m "Freeze publication v1.8 code, source data, and raw ABF archive"
git push origin "$BRANCH"

echo "Publication staging branch updated:"
echo "https://github.com/MartyBordeaux/NMDAR_Braess/tree/$BRANCH"
