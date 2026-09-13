# Publication version map

Current manuscript: **v1.5, 13 September 2026**  
Title: **Paradoxical NMDA response amplification in a robust receptor-state regime**

## Final evidence map

| Manuscript claim | Authoritative source |
|---|---|
| Representative ABF-derived measurement and level definition | `manuscript/source_data/figure1/` |
| -40 mV 5/5/1 classification | `data/derived/final/step15/04_cell_voltage_sweep_separation.csv` |
| Negative-voltage robustness | `data/derived/final/step15/07_cell_series_consistency.csv` |
| Historical proposal bounds and calibration implementation | `code/historical_step10/10_matched_parameter_rerouting.py` |
| Historical control targets | `data/derived/final/historical_control_calibration/experimental_targets.json` |
| Global strong-response accessibility | `data/derived/final/step18/` |
| Calibration/Braess frontier | `data/derived/final/step19/` |
| Exact local witness compatibility | `data/derived/final/step20/` |
| Expanded basin geometry | `data/derived/final/step21/` |
| Bayesian strong-state occupancy | `data/derived/final/step22/` |
| Threshold sensitivity and publication consolidation | `data/derived/final/step23/` |

## Superseded labels

The earlier ±2 pA plateau-neutral classification is retained only for provenance. Manuscript v1.5 uses threshold-free sweep separation and the term **unresolved** for overlapping sweep distributions.

The older Step4 `results_step4.zip` is a historical control-calibration stage and is not the same analysis as the later plateau-only Step04 target freeze. Its older Ro25 responder fields are not used for current experimental classification.

## Interpretation boundary

The Base model is a minimal mechanistic test model. Current results support the existence of a robust control-compatible Braess-like regime and strongly non-uniform occupancy relative to the original independent broad parameter measure. They do not establish a unique molecular mechanism of Ro25-6981 or identify a unique physiological 12-dimensional prior.
