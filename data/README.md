# Data included in this repository

The repository intentionally does **not** redistribute the raw Axon Binary Format (`.abf`) electrophysiology files.

The experimental data distributed here are the derived median inter-pulse plateau values used by the final analysis:

- `derived/plateau_by_cell_voltage.csv` — cell-series × voltage table of median same-day reference and Ro25 plateau levels, ratios, and QC fields.
- `derived/plateau_minus40_primary.csv` — frozen primary endpoint at -40 mV.
- `derived/potentiating_targets_minus40.csv` — five plateau-potentiating ratios used as continuous experimental magnitude targets.
- `derived/negative_voltage_robustness.csv` — robustness across negative holding potentials.

These are derived data, not raw recordings. They support reproduction of the reported endpoint summaries and manuscript figures without releasing the original ABF archive.

`model_outputs/` contains source tables for the Base, B-slow, and dual-timescale analyses reported in the manuscript.

## Important reproducibility note

The code for all analysis stages is included. Steps 1-3 require the original ABF recordings and are therefore provided for transparency rather than immediate public rerun. Downstream experimental summaries can start from `derived/plateau_by_cell_voltage.csv`.

A full de novo rerun of the model ensemble additionally requires the frozen Step10 accepted parameter ensemble (`matched_parameter_accepted.csv`) or the exact Step10 generator output. That frozen input was not present in the publication export used to initialize this repository. The reported model source tables and final numerical outputs used in the manuscript are included here; the Step10 ensemble should be added before claiming one-command model regeneration from scratch.
