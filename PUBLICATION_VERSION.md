# Publication version map

Current manuscript: **v1.5, 13 September 2026**  
Title: **Paradoxical NMDA response amplification in a robust receptor-state regime**

## Final evidence map

| Manuscript claim | Authoritative repository source |
|---|---|
| -40 mV threshold-free 5/5/1 classification | `data/derived/final/step15/minus40_sweep_separation.csv` |
| Negative-voltage robustness | `data/derived/final/step15/cell_series_consistency_negative.csv` and `voltage_summary_negative.csv` |
| Historical proposal bounds and calibration implementation | `code/final_v1_5/historical_step10/` |
| Historical control targets | `data/derived/final/historical_control_calibration/experimental_targets.json` |
| Global strong-response accessibility | `data/derived/final/step18/` |
| Calibration/Braess frontier | `data/derived/final/step19/` |
| Exact local witness compatibility | `data/derived/final/step20/` |
| Expanded basin geometry | `data/derived/final/step21/` |
| Bayesian strong-state occupancy | `data/derived/final/step22/` |
| Threshold sensitivity and publication consolidation | `data/derived/final/step23/` |
| Headline numerical integrity audit | `code/final_v1_5/verify_publication_v1_5.py` |

## Superseded labels and stages

The earlier +/-2 pA plateau-neutral classification is retained only for provenance. Manuscript v1.5 uses threshold-free sweep separation and the term **unresolved** for overlapping sweep distributions.

The older July `results_step4.zip` is a historical control-calibration stage and is not the same analysis as the later plateau-only Step04 target freeze. Its older Ro25 responder fields are not used for current experimental classification.

Positive nominal holding potentials are outside the current manuscript's voltage-robustness result.

## Interpretation boundary

The Base model is a minimal mechanistic test model. Current results support the existence of a robust control-compatible Braess-like regime and strongly non-uniform occupancy relative to the original independent broad parameter measure. They do not establish a unique molecular mechanism of Ro25-6981 or identify a unique physiological 12-dimensional prior.

The current normalized glutamate drive is dimensionless; no direct conversion to mM is defined by the model.
