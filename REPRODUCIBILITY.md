# Reproducibility

This document defines the reproducibility scope of manuscript v1.5, **Paradoxical NMDA response amplification in a robust receptor-state regime**.

## 1. Fast manuscript-level integrity audit

From the repository root:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r environment/requirements-v1.5.txt
python code/final_v1_5/verify_publication_v1_5.py
```

The verifier checks the frozen headline outputs without rerunning the expensive model searches. It confirms:

- -40 mV threshold-free classification: 5 potentiating / 5 suppressing / 1 unresolved;
- Step20 exact historical calibration replay and zero accepted/rejected classification mismatches;
- Step21 factor-1.5 strong+compatible fraction (~0.682);
- Step22 primary mean strong-state allocation (~0.458);
- Step23 frozen-threshold global counts: 1261/50,000 pre-calibration broad candidates and 1/5,000 accepted broad candidates strong somewhere.

## 2. Experimental result hierarchy

The final -40 mV classifier is `data/derived/final/step15/minus40_sweep_separation.csv`, not the earlier +/-2 pA neutral rule. Each series compares three Ro25 sweeps with three same-day same-voltage reference sweeps; sweeps are technical repeats. The final result is 5 potentiating, 5 suppressing and 1 unresolved series.

Only -80, -70, -60, -40 and -20 mV enter the voltage-robustness result. Among the ten series resolved at -40 mV, 43/50 within-series negative-voltage comparisons show complete separation in the same direction, seven are unresolved and none show opposite complete separation.

## 3. Historical Step10 proposal and calibration

The historical Step10 source snapshot under `code/final_v1_5/historical_step10/` defines the Latin-hypercube proposal, broad/reference supports, log/linear sampling flags, and exact control-calibration score. The associated frozen control targets are in `data/derived/final/historical_control_calibration/experimental_targets.json`.

Control calibration uses early/primary, late/primary and charge/primary control-shape ratios. Ro25 outcomes and responder labels are not calibration targets. The glutamate variable is a dimensionless normalized drive and has no defined direct conversion to mM in this model.

## 4. Final analysis chain

The publication-aligned chain is:

- Step18 — global prior geometry and full-grid strong-response confirmation;
- Step19 — calibration/Braess frontier;
- Step20 — exact historical calibration replay and local witness compatibility;
- Step21 — expanded basin geometry and multi-start search;
- Step22 — conditional Bayesian occupancy analysis;
- Step23 — threshold sensitivity and publication consolidation.

Compact final outputs required for manuscript claims are deposited under `data/derived/final/step18` ... `step23`. Source snapshots are kept under `code/final_v1_5/` and retain their original server-oriented discovery logic. They document the exact analyses, but a full heavy rerun still requires the corresponding frozen upstream tables and expected directory layout.

## 5. Figure 1 and raw-data provenance

The current manuscript uses an ABF-derived representative display at -40 mV to illustrate the measured inter-pulse level. The plotted example is explanatory; quantitative plateau values and all class assignments come from the independent sweep-level extraction pipeline rather than digitization of the figure.

Raw ABFs are not redistributed. Consequently, a clean clone can audit the deposited final numerical claims but cannot independently repeat raw waveform extraction, stimulus-artifact masking, baseline subtraction or raw-recording QC. A future raw-data archive should be linked as a separate deposit.

## 6. Acute memantine

The ancillary acute slice experiment used 30 micromolar memantine in the bath. It is distinct from chronic oral memantine treatment in the cited SCA1 study and from the primary Ro25 analysis.

## 7. Interpretation boundary

Fractions of sampled parameter sets or forcing contexts are measures of accessibility under specified computational sampling schemes, not biological prevalence. The Bayesian analysis estimates allocation conditional on the chosen model family and observed labels; it does not identify a unique biochemical prior.

## 8. Software environment

Core final-analysis dependencies are Python 3, NumPy, pandas, SciPy, Matplotlib and scikit-learn. Raw ABF extraction additionally requires pyABF. Historical source snapshots are preserved for provenance even where their original server path discovery must be adapted for a new machine.
