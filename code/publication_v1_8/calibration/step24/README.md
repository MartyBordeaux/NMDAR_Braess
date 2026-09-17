# NMDAR Braess Step 24 — calibration-window correction v1.0

This pipeline audits and corrects a window mismatch discovered in the historical Step10 control calibration.

- Experimental charge target: detrended current integrated over **140-500 ms**.
- Historical Step10 model `control_charge_PO_ms`: open-state occupancy accumulated over **0-600 ms**.
- Corrected model integral: **140-500 ms**, with all other calibration targets and weights frozen.

The main correction is computationally cheap. Because Step10 already stores mean model occupancy over 140-200 ms (`E`) and 200-500 ms (`L`), the historical rectangular-sum integral over 140-500 ms is exactly

`J_corrected = 60*E + 300*L` ms.

The pipeline therefore does **not** rerun 100,000 control trajectories or the 80-node Braess grid. It reuses the frozen candidate observables and previously confirmed response screens.

## Server launch

Copy this directory to the server, for example:

```bash
cd /root/nmda2
unzip NMDAR_Braess_step24_calibration_window_correction_v1_0.zip
cd NMDAR_Braess_step24_calibration_window_correction_v1_0
python -m pip install -r requirements.txt
./run_step24.sh
```

With the existing project layout, automatic discovery first tries:

- `/root/nmda/results_step10.zip::results/10_matched_rerouting/matched_parameter_all_samples.csv`
- `/root/nmda2/step_18/results_step_18_prior_geometry`
- `/root/nmda2/step_20/results_step_20_exact_witness_compatibility`
- `/root/nmda2/step_21/results_step_21_basin_population`
- `/root/nmda2/step_23/results_step_23_publication_consolidation`

Paths can be overridden explicitly:

```bash
./run_step24.sh \
  --step10 '/root/nmda/results_step10.zip::results/10_matched_rerouting/matched_parameter_all_samples.csv' \
  --step18 /root/nmda2/step_18/results_step_18_prior_geometry \
  --step20 /root/nmda2/step_20/results_step_20_exact_witness_compatibility \
  --step21 /root/nmda2/step_21/results_step_21_basin_population \
  --step23 /root/nmda2/step_23/results_step_23_publication_consolidation
```

Expected runtime is short relative to Steps 18-21 because no ODE grid is rerun; on the existing two-core server it should normally be on the order of minutes, dominated by reading/writing compressed tables.

## Main outputs

- `00_input_audit.json`
- `01_historical_score_replay.json` — hard gate that the historical score is reproduced.
- `02_corrected_candidate_scores.csv.gz` — old/new score, rank and retention status.
- `03_acceptance_shift_by_prior.csv`
- `04_primary_strong_support_shift.csv`
- `04A_threshold_strong_support_shift.csv`
- `05_step20_local_cloud_corrected_summary.csv`
- `06_step21_expanded_basin_corrected_summary.csv`
- `06A_step21_multistart_corrected_summary.json`
- `07_headline_comparison.csv`
- `08_scientific_decision.json`
- `README_RESULTS.md`
- Figures 24A-24C in PDF and PNG.

If the full Step10 table is unavailable, the pipeline can fall back to the Step19 merged broad table and run a broad-only audit. The final publication decision should use the full Step10 two-prior table.
