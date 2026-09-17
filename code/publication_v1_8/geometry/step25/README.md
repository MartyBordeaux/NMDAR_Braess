# NMDAR Braess Step 25 v1.0

## Purpose

Step 24 showed that correcting the control-calibration integral from the historical model window `0-600 ms` to the experiment-matched `140-500 ms` window materially changes the retained ensemble. In the broad prior, the corrected retained set contains 40 parameter vectors at or above the frozen experimental strong-response benchmark.

Step 25 validates those corrected strong vectors and asks a narrower question:

> Do corrected strong/control-compatible Base solutions occupy one sampled kinetic region, or several distinguishable kinetic regimes?

The pipeline does not refit the experiment and does not change experimental labels, the Base topology, the strong-response benchmark, or the forcing grid.

## Frozen inputs and conventions

The script auto-discovers the Step-24 result directory and the full Step-10 archive. The expected full Step-10 source is:

```text
/root/nmda/results_step10.zip::results/10_matched_rerouting/matched_parameter_all_samples.csv
```

The corrected Step-24 cohort is read from:

```text
02_corrected_candidate_scores.csv.gz
```

The frozen strong benchmark is:

```text
r_plateau >= 1.6359295439505144
```

The response grid is unchanged:

```text
G_peak = 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50, 0.65, 0.80, 0.95
tau_G  = 3, 5, 7.5, 10, 15, 20, 30, 40 ms
```

The Base response endpoint is the frozen nested median over 19 inter-pulse windows after discarding the first five intervals. Complete B-route removal disables both reversible pairs incident on `P1B`.

Historical control compatibility is replayed with the original Step-10 explicit-midpoint control simulator. The corrected integrated occupancy is

```text
J_140_500 = 60 * E_140_200 + 300 * L_200_500
```

and the exact corrected targets are reconstructed algebraically from the Step-24 score table. This gives an internal hard replay gate instead of duplicating rounded target constants.

## Analyses

1. Freeze corrected accepted broad and reference cohorts (5000 each).
2. Recover the exact corrected calibration targets and reproduce the corrected Step-24 scores.
3. Select the 40 corrected broad strong candidates.
4. Replay their historical control calibration with the source Step-10 simulator.
5. Recompute all 40 vectors over the full 80-node Base forcing grid at `dt=0.1 ms`.
6. Compare the recomputed maxima with the Step-24 stored maxima as a Base-model implementation gate.
7. Confirm every candidate at its recomputed maximum node using `dt=0.05 ms`.
8. Report raw control and blocked plateaus; no new denominator cutoff is applied.
9. Match the confirmed model responses to each of the five frozen potentiating experimental ratios.
10. Analyze 12-dimensional log-rate geometry using standardized distances, PCA, nearest neighbours, a minimum-spanning tree, and OPTICS sensitivity without fixing the number of clusters in advance.
11. Select distinct representative centres: best calibrated, experimental-typical match, extreme-response match, largest response, minimum-response match, and OPTICS cluster medoids where applicable.
12. Around each representative centre, perturb all 12 kinetic rates simultaneously by +/-5% and +/-10%. For every perturbed vector, recompute both the strong response and the corrected historical control-calibration score.
13. Issue a conservative sampled-region verdict. The output explicitly does not claim a proof of global basin connectivity.

## Run on the server

```bash
cd /root/nmda2
unzip NMDAR_Braess_step25_corrected_strong_geometry_v1_0.zip
cd NMDAR_Braess_step25_corrected_strong_geometry_v1_0
./run_step25.sh
```

Default output:

```text
/root/nmda2/step_25/results_step_25_corrected_strong_geometry
/root/nmda2/step_25/results_step_25_corrected_strong_geometry.zip
```

Default local sampling is 1000 points per representative centre at each perturbation level. On a small server it can be changed explicitly, for example:

```bash
STEP25_LOCAL_N=500 STEP25_THREADS=2 ./run_step25.sh
```

Use `STEP25_FORCE=1` only to discard an existing output directory and recompute it.

## Safety/replay gate

The final Base simulator is independently implemented from the frozen model equations, so Step 25 compares its full-grid maxima against the stored Step-24 maxima before interpreting downstream results. The default gate requires at least 95% of the 40 candidates to agree within either 0.01 absolute ratio units or 1% relative error, and requires sample 37318 to pass. A failed gate stops the pipeline.

The override

```bash
STEP25_ALLOW_REPLAY_MISMATCH=1 ./run_step25.sh
```

exists only for diagnostic inspection after a failed gate. Results from an overridden failed replay gate should not be used scientifically until the discrepancy is resolved.

## Principal outputs

```text
00_input_audit.json
01_corrected_target_reconstruction.json
02_corrected_accepted_cohorts.csv.gz
03_strong_candidates_step24.csv
03A_strong_calibration_replay_audit.csv
04_full_80_node_grid.csv.gz
04A_step24_response_replay_audit.csv
05_strong_candidates_confirmed.csv
06_experimental_ratio_matches.csv
07_denominator_diagnostics.csv
08_geometry_coordinates.csv.gz
09_geometry_pairwise_distances.csv
10_geometry_nearest_neighbors.csv
11_optics_sensitivity.csv
11A_optics_assignments.csv
11B_optics_pairwise_ari.csv
11C_strong_mst_edges.csv
11D_random_subset_compactness.csv
12_representative_centres.csv
13_local_perturbation_points.csv.gz
14_local_perturbation_summary.csv
15_scientific_decision.json
README_RESULTS.md
run_summary.json
Fig25A_response_denominator.{pdf,png}
Fig25B_pca_geometry.{pdf,png}
Fig25C_local_robustness.{pdf,png}
Fig25D_score_response.{pdf,png}
```

## Interpretation rule

`MULTIPLE_ROBUST_KINETIC_REGIONS_INDICATED` requires at least two canonical OPTICS groups, local +/-5% joint strong-and-compatible support in at least two of those groups, and reasonable stability across the OPTICS sensitivity grid. `SINGLE_DOMINANT_SAMPLED_REGION_COMPATIBLE` means the sampled geometry is consistent with one dominant group. Otherwise the result is `KINETIC_REGION_STRUCTURE_UNRESOLVED`.

None of these labels proves continuous global connectivity or identifies a physiological population distribution.
