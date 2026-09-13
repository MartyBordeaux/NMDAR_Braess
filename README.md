# NMDAR_Braess

Publication repository accompanying **Paradoxical NMDA response amplification in a robust receptor-state regime** (manuscript v1.5, 13 September 2026).

## Main scientific result

Archival Purkinje-cell recordings show heterogeneous Ro25-6981 responses at -40 mV: threshold-free sweep separation resolves 5 series as potentiating, 5 as suppressing, and 1 as unresolved. The five potentiating plateau ratios are 1.636, 1.695, 1.877, 1.891, and 4.073. Across the five negative holding potentials used for robustness analysis, the resolved primary direction is never reversed by complete sweep separation.

A minimal receptor-state Base model permits Braess-like amplification after removal of one ligand-entry route while the remaining kinetic parameters and input are held fixed. Experimental-scale amplification is rare under the original independent global sampling measure, but exact replay and local/expanded searches identify a finite control-compatible strong-response basin. The Bayesian analysis therefore estimates allocation to a model-defined strong-compatible regime; it does not recover a unique 12-dimensional biochemical prior.

The model result is a mechanistic possibility and robustness statement for the specified kinetic network. It is not proof of the molecular action of Ro25-6981, a graph theorem for arbitrary networks, or a direct biological prevalence estimate.

## Current repository hierarchy

- `data/derived/final/step15/` — final threshold-free experimental classification and negative-voltage robustness.
- `data/derived/final/step18/` ... `step23/` — compact final tables and decision summaries supporting the global/local/Bayesian claims.
- `data/derived/final/historical_control_calibration/` — historical control-shape targets used by Step10.
- `code/final_v1_5/` — frozen historical/final analysis source snapshots and a manuscript-level integrity verifier.
- `data/model_inputs/step10/` — previously deposited accepted Step10 ensemble.
- `data/model_outputs/` — retained Base/B-slow/dual supporting outputs from the earlier publication package.
- `manuscript/` — manuscript-version notes and provenance boundary for the current author draft.
- `provenance/` — historical claim ledgers and audit material.

Earlier files directly under `data/derived/` are retained for provenance. In particular, their old +/-2 pA `neutral` label is superseded by the threshold-free `unresolved` classification used in manuscript v1.5.

## Experimental endpoint

The primary endpoint is the baseline-subtracted negative inter-pulse current level during a 25-pulse train. At -40 mV, each drug/reference comparison contains three Ro25 sweeps and three same-day same-voltage reference sweeps. A series is potentiating when all three Ro25 plateau values are more inward than all three reference values, suppressing when all three are less inward, and unresolved when the sets overlap. Sweeps are technical repeats, not biological replicates.

The same-day reference is not claimed to be a verified within-cell pre-drug baseline. Only negative holding potentials are used for the voltage-robustness result in manuscript v1.5.

## Model, forcing and calibration

The Base model is a minimal receptor-state network with two ligand-entry routes feeding shared downstream gating/desensitization states. Its kinetic architecture is related to published NMDAR schemes, including Santucci & Raghavachari (PLoS Comput Biol. 2008;4:e1000208), but the current sampled association coordinates are effective first-order rates.

The current glutamate input is a dimensionless normalized drive. In the historical Step10 implementation an internal accumulator `g` is mapped to `G = g/(1+g)`. The repository defines no unique conversion of this `G` to mM.

Historical control calibration uses experimentally derived median shape targets: early/primary 1.100, late/primary 0.638, and charge/primary 257.22 ms. Fixed logarithmic tolerances are 0.35, 0.35 and 0.50, with weight 0.5 on the charge term. The tolerances are frozen calibration settings, not claimed physiological variances. Ro25 outcomes and responder labels were not used for historical control calibration.

## Acute memantine control

The ancillary acute slice experiment used **30 micromolar memantine in the bath**. Chronic oral memantine exposure is not part of this acute control.

## Quick verification

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r environment/requirements-v1.5.txt
python code/final_v1_5/verify_publication_v1_5.py
```

This checks the deposited headline results, including the 5/5/1 experimental classification, exact historical calibration replay, local basin fraction, primary Bayesian occupancy, and frozen-threshold global counts. It does not rerun the expensive historical ensembles.

See `REPRODUCIBILITY.md` for analysis scope and limitations.

## Data availability boundary

Raw ABF files are not redistributed. The repository provides the final derived measurements used for the manuscript-level conclusions. Raw-waveform extraction, stimulus-artifact masking and recording QC therefore cannot be repeated from this repository alone without the original ABF archive.

## License

Analysis software is released under the MIT License; see `LICENSE`.
