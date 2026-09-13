# NMDAR_Braess

Publication repository accompanying the manuscript **Paradoxical NMDA response amplification in a robust receptor-state regime** (manuscript v1.5, 13 September 2026).

## Scientific scope

The study combines cerebellar Purkinje-cell recordings with a minimal receptor-state model used as a mechanistic test model. The experimental endpoint is the baseline-subtracted negative inter-pulse current plateau during a 25-pulse train. At -40 mV, threshold-free sweep separation resolves 5 Ro25 series as potentiating, 5 as suppressing, and 1 as unresolved. The five potentiating plateau ratios are 1.636, 1.695, 1.877, 1.891, and 4.073.

The Base receptor-state topology permits Braess-like amplification after removal of the B entry route while control and blocked simulations share the same kinetic parameters and forcing. Strong experimental-scale amplification is rare under the original independent broad sampling measure, but exact local replay and expanded log-space searches identify a finite control-compatible strong basin around the accepted witness. Bayesian analysis therefore estimates allocation to a model-defined strong-compatible regime, not a unique 12-dimensional biochemical prior.

These model results establish a possible dynamical mechanism and its quantitative domain. They do **not** prove the molecular action of Ro25-6981, a graph-theoretic theorem for arbitrary networks, or a biological population distribution over kinetic parameters.

## Authoritative result hierarchy

The repository contains several historical stages. For the current manuscript use the following hierarchy:

1. `data/derived/final/step15/` — final threshold-free experimental classification and negative-voltage robustness.
2. `data/derived/final/step18/` — global prior geometry and full-grid strong-response confirmation.
3. `data/derived/final/step19/` — calibration/Braess frontier.
4. `data/derived/final/step20/` — exact historical calibration replay and local witness compatibility.
5. `data/derived/final/step21/` — expanded basin geometry and multi-start search.
6. `data/derived/final/step22/` — Bayesian physiological occupancy analysis.
7. `data/derived/final/step23/` — publication consolidation and threshold sensitivity.

Earlier `data/derived/plateau_minus40_primary.csv` and `negative_voltage_robustness.csv` are retained as historical plateau-freeze products. Their old ±2 pA neutral labeling is **not** the final classification used in manuscript v1.5.

## Experimental endpoint

For each sweep, the inter-pulse plateau is estimated after stimulus-artifact exclusion and baseline subtraction. At -40 mV, the final classification does not use a fixed ±2 pA neutral band. Instead, the three Ro25 sweep values are compared with the three same-day reference sweep values:

- potentiating: all Ro25 plateau values are more inward than all reference values;
- suppressing: all Ro25 values are less inward than all reference values;
- unresolved: the sweep distributions overlap.

The reference is same-day and same-voltage; it is not claimed to be a within-cell pre-drug baseline. Repeated sweeps are technical repeats, not independent biological replicates.

Only negative holding potentials are used for the voltage-robustness result in the current manuscript.

## Model and glutamate forcing

The model is a minimal receptor-state network with two ligand-entry routes feeding shared downstream gating/desensitization states. The downstream kinetic motif is related to published NMDAR state schemes, including Santucci & Raghavachari (PLoS Comput Biol. 2008;4:e1000208). The current implementation uses **effective first-order rates** in ms^-1.

Glutamate input is represented by a normalized dimensionless drive. In the historical Step10 implementation an internal accumulator `g` is mapped to `G = g/(1+g)`, and this dimensionless `G` multiplies the effective association rates. No unique conversion of the current `G` values to mM is defined; assigning such a concentration would require an additional concentration scale that is not part of this model.

## Historical control calibration

`code/historical_step10/10_matched_parameter_rerouting.py` is the original Step10 parameter-generation/scoring source now deposited for provenance. The control-calibration score uses three experimentally derived control-shape targets from `data/derived/final/historical_control_calibration/experimental_targets.json`:

- early/primary ratio: 1.100;
- late/primary ratio: 0.638;
- charge/primary ratio: 257.22 ms.

Fixed logarithmic tolerances are 0.35, 0.35, and 0.50, with weight 0.5 on the charge term. These tolerances are frozen calibration settings; the repository does not claim that they are independently measured physiological variances.

The original proposal used Latin-hypercube sampling. All coordinates except the Ro25 residual term were sampled logarithmically over their specified support. The `reference` and `broad` supports are encoded explicitly in the historical source.

## Acute memantine control

The ancillary acute slice experiment used **30 µM memantine in the bath**, consistent with the cited experimental source. The current manuscript does not describe chronic oral memantine exposure as part of this acute control.

## Repository layout

- `manuscript/` — v1.5 LaTeX source, bibliography, reporting scripts and figure sources.
- `code/historical_step10/` — original Step10 sampling/calibration implementation.
- `code/step_18/` ... `code/step_23/` — final analysis source snapshots used for the current global/local/Bayesian results.
- `data/derived/final/` — compact final tables and machine-readable decision summaries underlying manuscript claims.
- `data/model_inputs/step10/` — previously deposited retained Step10 ensemble.
- `code/model/` — compact Base/B-slow implementation retained from the earlier public repository.
- `provenance/` — result hierarchy and manuscript provenance notes.

Raw ABF recordings are not redistributed. The manuscript repository contains an ABF-derived representative display and sweep/record-level measurements for Figure 1; this is not equivalent to independent raw-waveform re-extraction.

## Reproduction levels

There are three distinct levels of reproducibility:

1. **Reporting reproduction** — rebuild manuscript tables and selected figures from frozen final outputs.
2. **Final-analysis rerun** — rerun Steps 18-23 from their required frozen upstream inputs.
3. **Raw experimental re-extraction** — requires the original ABF archive and is not possible from the public repository alone.

See `REPRODUCIBILITY.md` for commands and limitations.

## Manuscript status

Manuscript v1.5 is the current author draft. The repository is being aligned to that version; no journal acceptance or peer-review outcome is implied.

## License

Analysis software is released under the MIT License; see `LICENSE`. Experimental data remain subject to the provenance and redistribution limitations described above.
