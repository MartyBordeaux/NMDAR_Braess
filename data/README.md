# Data layout

The publication release contains both raw electrophysiology and manuscript source data.

- `raw/` — original Axon Binary Format (`.abf`) recordings used by the publication analyses, preserving the experimental archive structure. `raw/SHA256SUMS.txt` is the authoritative file manifest.
- `source_data_v1_8/` — machine-readable tables underlying the final manuscript summaries, figures, calibration, experiment-scale amplification, connectivity analysis, local perturbations and 27-scheme window-sensitivity analysis.
- `derived/`, `model_inputs/`, and `model_outputs/` — earlier analysis products retained for provenance.

Raw ABFs are required for a de novo rerun of trace extraction, baseline subtraction, file-level QC and the experimental window-sensitivity replay. The final source tables are provided so that manuscript-level numerical results can also be audited without rerunning every raw trace.

The final source-data tables correspond to the frozen v1.8 manuscript analysis. Older historical tables remain for provenance but should not override `source_data_v1_8/` when values differ.
