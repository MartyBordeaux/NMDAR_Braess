# NMDAR Braess Step 29 — corrected publication consolidation v1.0

## Purpose

Step 24 corrected the control-calibration integral from the historical model window `0–600 ms` to the experiment-matched window `140–500 ms`. Steps 25–28 then validated the resulting 40 strong/control-compatible vectors and tested their sampled connectivity.

Step 29 freezes that corrected result into one manuscript-facing package. It is intentionally **not** another topology search.

The pipeline does five things:

1. freezes the corrected accepted cohorts and the 40 numerically confirmed strong vectors;
2. attaches the final Step-28 topology labels;
3. selects publication representatives without centering the analysis on historical witness `37318`;
4. recomputes local simultaneous 12-rate robustness for the new representative set under the corrected calibration criterion;
5. generates final tables, figures, and a compact manuscript-number ledger.

No new experimental target, calibration tolerance, response threshold, denominator cutoff, or post-hoc exclusion is introduced.

## Frozen upstream state

The production run expects:

```text
/root/nmda2/step_24/results_step_24_calibration_window_correction
/root/nmda2/step_25/results_step_25_corrected_strong_geometry
/root/nmda2/step_28/results_step_28_residual_connectivity_stress
```

Alternative locations used by earlier package layouts are auto-discovered.

Hard gates require:

```text
Step24 status = MATERIAL_CHANGE
Step25 status = KINETIC_REGION_STRUCTURE_UNRESOLVED
Step28 status = PARTIAL_RESIDUAL_CONNECTIVITY_AFTER_FULLDIM_STRESS_TEST
confirmed strong vectors = 40/40
Step28 strict components = 36 + 3 + 1
Step28 native-only strict components = 35 + 3 + 1 + 1
historical witness 37318 present
Step28 engine replay gate = PASS
```

The strong threshold and strict corrected-calibration cutoff must match exactly between Steps 25 and 28.

## Final regime labels

The 40 confirmed strong/control-compatible vectors are assigned to four disjoint manuscript-facing groups:

```text
dominant_native_connected       35
dominant_flex_only_accession     1
secondary_native_connected       3
historical_witness_singleton     1
```

The first two together are the 36-vector dominant strict connected network. The distinction is retained because one member requires endpoint-preserving flexibility in the two auxiliary calibration coordinates to join that network.

## Publication representatives

Selection is deterministic.

`primary_typical_dominant_native`
: within the 35-vector dominant native-connected region, choose the model response closest in log-ratio distance to the median of the five experimental potentiating ratios.

`best_calibrated_dominant_native`
: minimum corrected calibration score within the same dominant native-connected region.

`minimum_experimental_match_dominant_native`
: closest dominant-native model response to the smallest experimental potentiating ratio.

`extreme_experimental_match_dominant_native`
: closest dominant-native model response to the largest experimental potentiating ratio.

`largest_response_dominant_native`
: largest confirmed response within the dominant native-connected region. Its denominator is retained explicitly; no denominator cutoff is added.

`secondary_region_representative`
: best calibrated member of the secondary native-connected three-vector region.

`flex_aux_only_accession`
: the vector that joins the dominant 36-vector strict network only under endpoint-preserving auxiliary-coordinate flexibility.

`historical_witness_sensitivity`
: sample `37318`, retained as a historical/sensitivity solution rather than as the primary representative.

With the frozen Step-25/28 inputs, the primary representative should be sample `30902`.

## Local robustness recomputation

By default Step 29 recomputes local robustness for four centres:

```text
primary typical dominant-native representative
best calibrated dominant-native representative
secondary-region representative
historical witness 37318
```

For each centre, all 12 kinetic rates are perturbed simultaneously and independently by uniform multiplicative factors at:

```text
±5%
±10%
```

Default sample size:

```text
1000 points per centre per perturbation level
```

For every perturbed point:

- the centre's frozen forcing node is used for the Base/B-deleted response calculation;
- response integration uses `dt = 0.05 ms`;
- corrected control calibration uses `dt = 0.05 ms`;
- the inherited historical calibration-drive coordinates for that centre are kept fixed;
- strong response requires the frozen experimental threshold;
- compatibility requires the strict corrected cutoff;
- no denominator cutoff is applied.

This recomputation replaces the old witness-centric robustness presentation with a direct comparison of the new primary region, the secondary region, and the historical witness.

## Run

```bash
cd /root/nmda2
unzip NMDAR_Braess_step29_corrected_publication_consolidation_v1_0.zip
cd NMDAR_Braess_step29_corrected_publication_consolidation_v1_0

STEP29_THREADS=2 ./run_step29.sh
```

Background run:

```bash
STEP29_THREADS=2 nohup ./run_step29.sh > out.log 2>&1 &
```

Monitor:

```bash
tail -f out.log
```

If an interrupted run already created the output directory:

```bash
STEP29_THREADS=2 STEP29_RESUME=1 nohup ./run_step29.sh > out.log 2>&1 &
```

For a clean replacement:

```bash
STEP29_THREADS=2 STEP29_FORCE=1 nohup ./run_step29.sh > out.log 2>&1 &
```

A fast consolidation without new local perturbations is available only for diagnostics:

```bash
STEP29_SKIP_LOCAL=1 ./run_step29.sh
```

The publication run should not use `STEP29_SKIP_LOCAL=1`.

## Output

Default directory:

```text
/root/nmda2/step_29/results_step_29_corrected_publication_consolidation
```

Archive:

```text
/root/nmda2/step_29/results_step_29_corrected_publication_consolidation.zip
```

Principal files:

```text
00_input_audit.json
00_source_manifest.csv
01_frozen_corrected_accepted_cohorts.csv.gz
02_final_strong40_topology.csv
03_final_verified_network_edges.csv
03A_final_residual_barriers.csv
04_experimental_model_matches.csv
05_publication_representatives.csv
06_representative_local_perturbation_points.csv.gz
06A_representative_local_robustness_summary.csv
07_manuscript_headline_numbers.csv
08_scientific_decision.json
MANUSCRIPT_NUMBERS.md
SCIENTIFIC_INTERPRETATION.md
README_RESULTS.md
run_summary.json
```

Figures:

```text
Fig29A_corrected_support.pdf
Fig29B_response_denominator_topology.pdf
Fig29C_experimental_ratio_matches.pdf
Fig29D_connectivity_composition.pdf
Fig29E_primary_response_map.pdf
Fig29F_representative_local_robustness.pdf
```

PNG copies are generated for inspection.

## Interpretation boundary

Step 29 is a consolidation of the finite sampled evidence. It does not convert the residual Step-28 connectivity failures into a proof of global topological disconnection. The correct manuscript-level result is that multiple corrected strong-compatible solutions exist; most sampled solutions belong to one verified connected network, while a small secondary region and the historical witness remained separate under the tested finite path searches.
