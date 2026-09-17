# NMDAR Braess Step 27 v1.0

## Scientific purpose

Step 26 found that direct log-linear paths did not establish global connectivity of the 40 corrected strong/control-compatible Base-model vectors. It also showed that this negative result cannot be interpreted as proof of multiple globally disconnected basins because only a finite set of direct paths was tested.

Step 27 therefore asks a narrower and more informative question:

> Can the Step-26 sampled components be connected by endpoint-preserving curved paths that remain simultaneously strong and corrected-control compatible?

The analysis deliberately separates two path families.

### Primary: `native_aux`

Only the 12 kinetic rates bend. The two historical control-calibration drive coordinates remain on their geometric interpolation between the two endpoint values.

This is the primary kinetic-connectivity test.

### Secondary: `flex_aux`

The 12 kinetic rates bend as above, and the two calibration-drive coordinates may also bend smoothly. Their values are still fixed exactly at both endpoints.

This is a continuous 14-dimensional parameter-space sensitivity. It is more permissive than `native_aux`, but unlike the Step-26 fixed-A/fixed-B schedules it still terminates at both stored endpoint parameter vectors.

Step 27 does **not** use fixed-A or fixed-B whole-edge schedules as topological connectors, because those schedules generally do not end at the stored auxiliary coordinates of both candidate vectors.

## Frozen inputs

The pipeline auto-discovers:

```text
/root/nmda2/step_25/results_step_25_corrected_strong_geometry
/root/nmda2/step_26/results_step_26_strong_connectivity
/root/nmda/results_step10.zip
```

The Step-10 full broad table is used only to enforce the original broad-prior parameter box during path optimization. The corrected accepted broad ensemble from Step 25 defines the standardized log-rate geometry and PCA directions.

Hard gates require:

```text
Step25 status = KINETIC_REGION_STRUCTURE_UNRESOLVED
Step26 status = SAMPLED_CONNECTIVITY_NOT_ESTABLISHED
confirmed strong vectors = 40/40
```

Frozen scientific thresholds are inherited from Step 25/26:

```text
strong threshold = 1.6359295439505144
corrected broad midpoint cutoff ~= 0.0806508971
strict max-accepted cutoff ~= 0.0806451494
```

No denominator cutoff is introduced.

## Curved path parameterization

Rates are transformed to standardized log coordinates using the 5000 corrected accepted broad vectors. A path between endpoints A and B begins with the straight line

```text
z(lambda) = (1-lambda) z_A + lambda z_B
```

and is bent by a sine basis aligned to the leading PCA directions of the corrected accepted broad ensemble:

```text
z(lambda) = z_linear(lambda)
          + sum_m sin(m*pi*lambda) * sum_p a_mp * PC_p
```

Because every sine term is zero at lambda=0 and lambda=1, the endpoint kinetic vectors are preserved exactly.

Default search basis:

```text
5 PCA directions
2 sine modes
10 kinetic curve coefficients
```

For `flex_aux`, four additional coefficients bend the two log auxiliary coordinates with the same two sine modes while preserving their endpoints.

## Search strategy

Step 27 starts from the **native** coarse components from Step 26. The current Step-26 result contains seven such sampled components.

First, a dt=0.05-ms native spanning forest is certified inside every Step-26 coarse component. Already refined Step-26 edges are reused when possible; only additional edges needed for a certified forest are recomputed.

Potential inter-component links are then prioritized by minimum standardized log-rate distance. For each selected component pair, the nearest endpoint pairs are tried in this order:

1. direct endpoint-preserving log-linear path;
2. curved `native_aux` kinetic path;
3. curved `flex_aux` path, if enabled.

Curve coefficients are optimized with a deterministic seeded evolutionary minimax search. Calibration is optimized first because it is the dominant Step-26 barrier. If a calibration-compatible candidate curve loses strong support, a smaller joint rescue search adds the strong-response constraint.

## Final verification

Search success is never accepted directly.

Every retained bridge is re-evaluated on 41 path points. At each point:

1. corrected control calibration is evaluated at `dt=0.05 ms`;
2. the full frozen 80-node Base response grid is scanned at `dt=0.1 ms`;
3. the best node is reconfirmed at `dt=0.05 ms`;
4. if it falls below the strong threshold, a complete 80-node `dt=0.05 ms` fallback is performed.

A path passes only if **every** one of its verification points is both:

```text
calibration_score <= corrected cutoff
r_plateau >= strong threshold
```

The minimum control plateau along the support nodes is retained as a diagnostic but is never used as an exclusion rule.

## Decision labels

```text
CURVED_KINETIC_CONNECTIVITY_ESTABLISHED_NATIVE_AUX
```

All 40 vectors are connected by a verified network whose bridges require no bending of auxiliary calibration coordinates.

```text
CURVED_FULL_PARAMETER_CONNECTIVITY_ESTABLISHED
```

All 40 vectors are connected, but at least one verified bridge requires a continuous endpoint-preserving bend of the auxiliary calibration coordinates.

```text
PARTIAL_CURVED_CONNECTIVITY_ONLY
```

One or more curved bridges are verified, but a single connected network is not obtained within the frozen search budget.

```text
CURVED_CONNECTIVITY_NOT_ESTABLISHED_WITHIN_BUDGET
```

No inter-component bridge is verified within the frozen search budget.

The last two labels do not prove global disconnection.

## Run

```bash
cd /root/nmda2
unzip NMDAR_Braess_step27_curved_bridge_optimization_v1_0.zip
cd NMDAR_Braess_step27_curved_bridge_optimization_v1_0
./run_step27.sh
```

Default output:

```text
/root/nmda2/step_27/results_step_27_curved_bridge_optimization
/root/nmda2/step_27/results_step_27_curved_bridge_optimization.zip
```

On the two-core server:

```bash
STEP27_THREADS=2 ./run_step27.sh
```

If an interrupted run has already created checkpoints:

```bash
STEP27_RESUME=1 ./run_step27.sh
```

To discard all existing Step-27 output and restart:

```bash
STEP27_FORCE=1 ./run_step27.sh
```

The default search budget is:

```text
13 search points per candidate path
41 final verification points
24 evolutionary proposals per generation
18 calibration-search generations
10 optional joint-rescue generations
12 component-pair links maximum
3 endpoint pairs per component-pair link
```

The most useful budget controls are:

```bash
STEP27_MAX_COMPONENT_PAIRS=18 STEP27_PAIR_ATTEMPTS=5 ./run_step27.sh
```

Only increase these after inspecting the default result.

## Principal outputs

```text
00_input_audit.json
01_step26_native_components.csv
02_component_pair_priority.csv
03_internal_native_certification.csv
04_bridge_attempts.csv
05_verified_bridges.csv
06_verified_path_points.csv.gz
07_verified_network_edges.csv
08_final_components.csv
09_network_bottlenecks.csv
10_scientific_decision.json
README_RESULTS.md
run_summary.json
Fig27A_verified_network.{pdf,png}
Fig27B_bridge_barriers.{pdf,png}
Fig27C_verified_bridge_profiles.{pdf,png}
```

The `_checkpoints/` directory stores verified paths and optimization histories so an interrupted run can be resumed.
