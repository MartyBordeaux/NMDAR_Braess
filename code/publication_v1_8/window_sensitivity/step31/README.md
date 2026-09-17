# NMDAR Braess Step 31 — final analysis-window robustness

This pipeline finalizes the sensitivity analysis of the operational response windows used for control calibration and regenerates the manuscript-facing parameter-space geometry.

## What is fixed in this version

The calibration input is restricted to untreated `SlowEPSC` reference recordings. Files from `+Ro25` and `+Mem`/memantine conditions are hard-excluded at raw-ABF resolution. The output mapping records these condition checks explicitly and the run stops if a drug-condition file enters calibration.

The five frozen baseline series-level ratios remain the absolute experimental calibration reference. Raw ABFs are used only to measure the *within-series relative change* caused by moving an operational window boundary. Thus, for each series and observable,

`anchored alternative ratio = frozen baseline ratio × (raw alternative ratio / raw baseline ratio)`.

At the baseline scheme, all five series are required to reproduce their frozen baseline ratios after anchoring. This is checked series by series, not only at the median-target level.

The same paired-ratio principle is retained for model observables: exact frozen baseline values define the absolute reference, while the Step 31 integrator supplies only the relative change under moved window boundaries.

## Prespecified window grid

The post-train start remains fixed at 140 ms. The other boundaries are varied as:

- split: 180, 200, 220 ms
- primary-window end: 230, 250, 270 ms
- late/integral end: 450, 500, 550 ms

This gives 27 schemes. The baseline is 140–200 / 140–250 / 200–500 / 140–500 ms.

## Main gates

The run stops if any of the following fails:

1. a calibration ABF resolves to Ro25 or memantine;
2. a calibration ABF is not a SlowEPSC recording;
3. any of the five anchored baseline series ratios differs from its frozen value;
4. the model prefix integrator differs grossly from the frozen baseline observables;
5. the anchored baseline model score does not reproduce the frozen calibration score;
6. the baseline top-5000 accepted sets do not reproduce the frozen sets.

## Main outputs

- `00_input_audit.json`
- `01_experimental_file_mapping.csv`
- `02A_integration_method_replay.csv`
- `02_baseline_series_raw_replay.csv`
- `02B_baseline_series_anchored_replay.csv`
- `03_window_schemes.csv`
- `04_experimental_targets_by_window.csv`
- `04A_experimental_series_metrics_by_window.csv.gz`
- `05_model_window_sensitivity_by_prior.csv`
- `06_final40_window_retention.csv`
- `07_parameter_space_coordinates.csv.gz`
- `08_geometry_edges.csv`
- `09_representative_path_profile.csv`
- `10_scientific_summary.json`
- `Fig31A_window_sensitivity.pdf/png`
- `Fig31B_parameter_space_geometry.pdf/png`
- `Fig31C_rate_filtering_marginals.pdf/png`
- `README_RESULTS.md`

`Fig31C` uses density curves only for the 50,000 sampled broad models and the 5,000 control-compatible models. The 40 experiment-scale solutions are shown as rug marks rather than a density estimate.

## Run

```bash
cd /root/nmda2
unzip NMDAR_Braess_step31_window_robustness_final_v1_0_0.zip
cd NMDAR_Braess_step31_window_robustness_final_v1_0_0

STEP31_THREADS=2 \
STEP31_RAW_ROOT=/root/nmda/IV_NMDA \
STEP31_FORCE=1 \
nohup ./run_step31.sh > out.log 2>&1 &
```

Monitor:

```bash
tail -f out.log
```

Resume an interrupted run:

```bash
STEP31_THREADS=2 \
STEP31_RAW_ROOT=/root/nmda/IV_NMDA \
STEP31_RESUME=1 \
nohup ./run_step31.sh > out.log 2>&1 &
```

Results are written to:

`/root/nmda2/step_31/results_step_31_window_robustness_final`

and archived as:

`/root/nmda2/step_31/results_step_31_window_robustness_final.zip`
