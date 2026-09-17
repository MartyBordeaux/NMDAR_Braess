# NMDAR Braess Step 28 v1.0

## Purpose

Step 27 established a verified connected network containing 36 of the 40 corrected strong/control-compatible vectors. Three residual blocks remained:

```text
{27498, 29959}
{6113}
{37318}
```

Step 28 is a targeted final connectivity stress test for those residual blocks. It removes the main geometric restriction of Step 27: paths are no longer confined to a low-dimensional PCA sine basis.

The primary question is:

> Can the residual Step-27 components be connected to the verified network, or to each other, by endpoint-preserving multi-waypoint paths in the full 12-dimensional kinetic log-rate space while remaining simultaneously strong and corrected-control compatible?

A negative result remains a finite-search result, not a proof of global disconnection.

## Frozen inputs

The pipeline auto-discovers:

```text
/root/nmda2/step_25/results_step_25_corrected_strong_geometry
/root/nmda2/step_27/results_step_27_curved_bridge_optimization
/root/nmda/results_step10.zip
```

Hard gates require:

```text
Step25 status = KINETIC_REGION_STRUCTURE_UNRESOLVED
Step27 status = PARTIAL_CURVED_CONNECTIVITY_ONLY
confirmed strong vectors = 40/40
Step27 strict component sizes = 36 + 2 + 1 + 1
residual blocks = {27498,29959}, {6113}, {37318}
all Step27 network edges pass the strict corrected cutoff
```

The original full broad Step-10 candidate table is used only to recover the original broad-prior box. The 5000 corrected accepted broad vectors define standardized log-rate coordinates and provide calibration-manifold search seeds.

## Path families

### Primary: `native_aux`

The 12 kinetic rates are represented in standardized log space. A candidate path contains 2 or 3 freely optimized internal waypoints by default. Consecutive waypoints are joined by log-linear segments.

Thus, for each segment, the kinetic path is continuous and all endpoints are preserved exactly. Internal waypoints are constrained to the original broad-prior box.

The historical calibration-drive coordinates

```text
pulse_amplitude_au
glutamate_tau_ms
```

remain on their endpoint-preserving geometric interpolation.

This is the primary kinetic-connectivity stress test.

### Secondary: `flex_aux`

The 12 kinetic waypoints are optimized as above, but the two auxiliary calibration coordinates also receive freely optimized internal waypoints. Their endpoint values remain fixed exactly.

This is a continuous 14-dimensional sensitivity analysis. A bridge requiring `flex_aux` is not described as a purely kinetic connector.

## Why this is stronger than Step 27

Step 27 used a low-dimensional curved path:

```text
5 PCA directions x 2 sine modes = 10 kinetic curve coefficients
```

Step 28 uses full-dimensional internal waypoints:

```text
2 waypoints -> 24 free kinetic coordinates
3 waypoints -> 36 free kinetic coordinates
```

and, for the optional `flex_aux` family:

```text
+4 or +6 auxiliary waypoint coordinates
```

No PCA-direction restriction is imposed on the path itself.

## Search initialization

The optimizer is deterministic for a fixed seed and uses several predefined initialization families:

1. the direct log-linear path;
2. waypoint seeds taken from nearby known strong vectors;
3. waypoint seeds taken from the corrected accepted broad ensemble near the direct path;
4. stochastic exploration within the original broad-prior box.

These seeds only guide optimization. They do not alter the acceptance criteria.

## Search objective

For each proposal, corrected control calibration is evaluated along the path. A subset of the most promising proposals is also screened for strong response on the frozen forcing nodes represented among the 40 confirmed strong vectors.

The minimax search quantity is approximately

```text
max(
    max_path_calibration_score / strict_cutoff,
    strong_threshold / min_path_response
)
```

with only a very small path-length regularizer to avoid unnecessary detours.

After optimization, the best candidate is screened once on the complete frozen 80-node forcing grid before deciding whether final verification is warranted.

## Final verification

Search success is never accepted directly.

Every retained bridge is re-evaluated on 61 path points. At each point:

1. corrected control calibration is evaluated at `dt = 0.05 ms`;
2. the full frozen 80-node Base-response grid is scanned at `dt = 0.1 ms`;
3. the best forcing node is reconfirmed at `dt = 0.05 ms`;
4. if the reconfirmed response falls below threshold, the complete 80-node grid is rescanned at `dt = 0.05 ms`.

The primary bridge criterion is the **strict** corrected calibration cutoff:

```text
calibration_score <= strict corrected cutoff
r_plateau >= 1.6359295439505144
```

at every verification point.

No denominator cutoff or post-hoc denominator exclusion is introduced.

## Default search scope

All six pairings among the four Step-27 strict components are eligible. They are prioritized by the best Step-27 near-feasible joint barrier when available, then by standardized log-rate distance.

For each component pair, endpoint pairs are selected using both:

- the best Step-27 failed attempts;
- nearest full-dimensional endpoint distances.

Default budget:

```text
17 search points
61 verification points
32 proposals per generation
20 generations
2 restarts
10 response-evaluated proposals per generation
6 endpoint pairs per component pair
2 endpoint pairs receive the deeper 3-waypoint search
all 6 Step-27 component pairs may be tested
```

## Run

```bash
cd /root/nmda2
unzip NMDAR_Braess_step28_residual_connectivity_stress_v1_0.zip
cd NMDAR_Braess_step28_residual_connectivity_stress_v1_0

STEP28_THREADS=2 ./run_step28.sh
```

Background run:

```bash
STEP28_THREADS=2 nohup ./run_step28.sh > out.log 2>&1 &
```

Monitor:

```bash
tail -f out.log
```

Resume an interrupted run:

```bash
STEP28_THREADS=2 STEP28_RESUME=1 nohup ./run_step28.sh > out.log 2>&1 &
```

Discard existing Step-28 output and restart:

```bash
STEP28_THREADS=2 STEP28_FORCE=1 nohup ./run_step28.sh > out.log 2>&1 &
```

## Useful budget controls

A slightly larger search without changing scientific definitions:

```bash
STEP28_POP_SIZE=40 \
STEP28_GENERATIONS=28 \
STEP28_RESTARTS=3 \
STEP28_PAIR_ATTEMPTS=8 \
STEP28_RESPONSE_TOP=12 \
STEP28_THREADS=2 \
./run_step28.sh
```

A 4-waypoint escalation is available but is intentionally not the default:

```bash
STEP28_MAX_WAYPOINTS=4 STEP28_THREADS=2 ./run_step28.sh
```

## Output

Default directory:

```text
/root/nmda2/step_28/results_step_28_residual_connectivity_stress
```

Archive:

```text
/root/nmda2/step_28/results_step_28_residual_connectivity_stress.zip
```

Principal files:

```text
00_input_audit.json
01_step27_residual_components.csv
02_component_pair_priority.csv
03_search_attempts.csv
04_search_history.csv.gz
05_verified_bridges.csv
06_verified_path_points.csv.gz
07_final_network_edges.csv
08_final_components.csv
09_residual_barrier_summary.csv
10_scientific_decision.json
README_RESULTS.md
run_summary.json
Fig28A_verified_network.{pdf,png}
Fig28B_search_barriers.{pdf,png}
Fig28C_verified_bridge_profiles.{pdf,png}
```

The `_checkpoints/` directory stores completed endpoint/mode/waypoint stages for resume.

## Decision labels

```text
FULL12D_MULTIWAYPOINT_CONNECTIVITY_ESTABLISHED_NATIVE_AUX
```

All 40 vectors are connected using the Step-27 network plus new bridges that require only full-12D kinetic multi-waypoint paths.

```text
FULL14D_MULTIWAYPOINT_CONNECTIVITY_ESTABLISHED_WITH_FLEX_AUX
```

All 40 vectors become connected, but at least one new bridge requires endpoint-preserving auxiliary-coordinate bending.

```text
PARTIAL_RESIDUAL_CONNECTIVITY_AFTER_FULLDIM_STRESS_TEST
```

At least one new residual bridge is verified, but the 40-vector network remains split.

```text
RESIDUAL_CONNECTIVITY_NOT_ESTABLISHED_AFTER_FULLDIM_STRESS_TEST
```

No new residual bridge is verified within the frozen search budget.

The last two outcomes are not mathematical proofs of disconnected basins.
