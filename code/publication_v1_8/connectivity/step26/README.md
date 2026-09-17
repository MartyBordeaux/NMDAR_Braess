# NMDAR Braess Step 26 v1.0

## Purpose

Step 25 confirmed **40/40** corrected broad-prior strong vectors at `dt=0.05 ms`, but the 12-dimensional point-cloud analysis did not resolve whether they form one kinetic region or several. The canonical OPTICS result contained one small cluster plus substantial noise and the clustering assignments were unstable across the sensitivity grid. Step 26 therefore replaces cluster labels with an explicit path-connectivity test.

The scientific question is:

> Can the 40 corrected strong/control-compatible solutions be linked by continuous tested paths that remain both strong and corrected-control compatible?

This is a finite-resolution sampled-path analysis. A connected spanning network is evidence for one sampled connected region along the tested paths. Failure to find such a network is **not** a proof of global disconnection.

## Frozen inputs

The pipeline auto-discovers the Step-25 result directory, normally:

```text
/root/nmda2/step_25/results_step_25_corrected_strong_geometry
```

It requires the Step-25 status:

```text
KINETIC_REGION_STRUCTURE_UNRESOLVED
```

and verifies the Step-25 numerical gate before doing new calculations.

Frozen quantities are inherited directly from Step 25:

```text
strong threshold = 1.6359295439505144
corrected broad midpoint cutoff = 0.0806508971425229
strict max-accepted cutoff = 0.0806451494470612
```

The corrected historical control calibration remains based on:

```text
J_140_500 = 60 * E_140_200 + 300 * L_200_500
```

The final Base response grid remains:

```text
G_peak = 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50, 0.65, 0.80, 0.95
tau_G  = 3, 5, 7.5, 10, 15, 20, 30, 40 ms
```

No new control-denominator cutoff is introduced.

## Path definition

For every tested edge between two Step-25 strong vectors, all 12 kinetic rates are interpolated geometrically:

```text
log k(lambda) = (1-lambda) log k_A + lambda log k_B
```

Three prespecified continuous schedules are used for the two historical calibration-drive coordinates (`pulse_amplitude_au`, `glutamate_tau_ms`):

1. `interp_aux` — geometric interpolation between endpoint auxiliary values;
2. `fixed_A_aux` — endpoint-A auxiliary values held fixed over the entire path;
3. `fixed_B_aux` — endpoint-B auxiliary values held fixed over the entire path.

Primary `native` connectivity uses `interp_aux` only. The secondary `any_schedule` criterion accepts an edge only if **one single schedule** passes at every path point. It does not permit pointwise switching among schedules.

## Edge search

The full 40-point graph contains 780 possible edges. Step 26 prioritizes:

1. the Step-25 minimum-spanning-tree edges;
2. the union of 4-nearest-neighbour edges;
3. direct pairs among Step-25 representative centres;
4. remaining edges in increasing standardized log-rate distance.

By default the first 180 edges are screened. This is configurable with `STEP26_MAX_EDGES`.

The coarse screen uses 21 path points. It first tests the forcing nodes that were maxima for the 40 Step-25 strong candidates; a path point not proven strong there is automatically rescanned over the remaining nodes of the full 80-node grid. Therefore a coarse `strong` failure has already received a complete frozen-grid search.

If the coarse passing graph contains a spanning tree, only the edges needed to certify connectivity are refined to 41 path points. Refined paths are evaluated on the **full 80-node grid at every point** and every path point is reconfirmed at `dt=0.05 ms` at its best `dt=0.1 ms` node, with a full `dt=0.05 ms` grid fallback if needed.

## Decision hierarchy

The primary result is one of:

```text
ONE_SAMPLED_CONNECTED_REGION_SUPPORTED_NATIVE_INTERPOLATION
```

All 40 vectors are connected by refined paths using interpolated historical calibration-drive coordinates.

```text
ONE_SAMPLED_KINETIC_REGION_SUPPORTED_WITH_PRESPECIFIED_AUX_SCHEDULES
```

Native interpolation alone is insufficient, but all 40 vectors can be connected when each edge may use one of the three prespecified whole-edge auxiliary schedules.

```text
SAMPLED_CONNECTIVITY_NOT_ESTABLISHED
```

No complete refined spanning network was found within the tested edge budget. This label does not claim global disconnection.

The strict-cutoff result is reported separately and does not replace the midpoint-cutoff primary criterion.

## Run on the server

```bash
cd /root/nmda2
unzip NMDAR_Braess_step26_strong_connectivity_v1_0.zip
cd NMDAR_Braess_step26_strong_connectivity_v1_0
./run_step26.sh
```

Default output:

```text
/root/nmda2/step_26/results_step_26_strong_connectivity
/root/nmda2/step_26/results_step_26_strong_connectivity.zip
```

For the two-core server, leave:

```bash
STEP26_THREADS=2 ./run_step26.sh
```

A larger edge budget can be requested explicitly:

```bash
STEP26_MAX_EDGES=300 STEP26_THREADS=2 ./run_step26.sh
```

Use `STEP26_FORCE=1` only to replace an existing output directory.

## Principal outputs

```text
00_input_audit.json
01_strong_candidates.csv
02_edge_pool.csv
03_coarse_edge_summary.csv
04_coarse_path_points.csv.gz
05_coarse_components.csv
06_refined_edge_summary.csv
07_refined_path_points.csv.gz
08_refined_components.csv
09_representative_connectivity.csv
10_bottleneck_summary.csv
11_scientific_decision.json
README_RESULTS.md
run_summary.json
Fig26A_connectivity.{pdf,png}
Fig26B_bridge_bottlenecks.{pdf,png}
Fig26C_path_denominator.{pdf,png}
Fig26D_representative_paths.{pdf,png}   # only when direct representative paths were refined
```

## Important interpretation constraint

The calibration score constrains normalized control-shape descriptors, not absolute open-state amplitude. Step 25 showed that several strong vectors have very small control plateaus. Step 26 therefore carries control plateaus through every refined path as a diagnostic but does not invent a post hoc denominator threshold.

The pipeline addresses sampled connectivity only. It does not infer a physiological distribution of kinetic parameters and does not establish a mathematically complete global basin topology.
