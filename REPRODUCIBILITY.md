# Reproducibility

This document defines the reproducibility scope of the final Neuropharmacology author manuscript (publication analysis v1.8, September 2026).

## 1. Experimental data

The raw electrophysiology source is the Axon Binary Format archive deposited under `data/raw/`. The publication release preserves the original relative folder structure and includes `data/raw/SHA256SUMS.txt` as the authoritative raw-file manifest.

The primary Ro25 result is based on complete separation of the technical sweeps at -40 mV: 5 series are potentiating, 5 suppressing and 1 unresolved. Negative-voltage robustness is evaluated at -80, -70, -60, -40 and -20 mV. Control-shape calibration metrics are extracted from untreated SlowEPSC reference recordings; Ro25 and memantine recordings are excluded from calibration.

The post-train experimental observables are computed after baseline subtraction as early mean current (140-200 ms), primary mean current (140-250 ms), late mean current (200-500 ms), and integrated current over 140-500 ms. The corresponding control targets are formed as within-series ratios and then aggregated across the five control series.

## 2. Model and calibration

The Base receptor-state model is implemented under `code/model/` and the final publication engine is deposited under `code/publication_v1_8/model/`. Route removal is implemented by making the B-entry route unavailable while keeping all other kinetic rates, the forcing protocol and initial conditions fixed.

The final calibration observables are

- `E = <P_O>_140-200`,
- `H = <P_O>_140-250`,
- `L = <P_O>_200-500`,
- `J = integral_140^500 P_O(t) dt`.

The control targets are 1.0998943768 for E/H, 0.6375504980 for L/H and 257.2199452 ms for J/H. The conservative minimum log-tolerances are 0.35, 0.35 and 0.50, respectively, with weight 0.5 on the integral term. Ro25 response labels are not used to select the control-calibrated ensemble.

The experiment-scale amplification benchmark is `r_plateau >= 1.63592954395`, the smallest observed potentiating ratio at -40 mV.

## 3. Final computational chain

The publication-aligned sequence is:

1. replay/generate the Base-model parameter candidates;
2. evaluate the final control-shape calibration score and retain the calibrated ensemble;
3. confirm amplification on the frozen 80-node forcing grid;
4. test sampled connectivity with Steps 26-28;
5. consolidate representatives, experimental-ratio matches and local perturbation summaries with Step 29;
6. evaluate sensitivity to 27 predefined post-train window schemes with Step 31.

The final scripts are grouped under `code/publication_v1_8/`. Earlier code snapshots remain in the repository only for provenance.

## 4. Connectivity interpretation

A verified parameter path is an interpolation between two complete kinetic parameter sets for which every evaluated intermediate model remains both control-compatible and experiment-scale amplifying. The final sampled endpoint-preserving graph contains a dominant connected component of 36 solutions, a smaller component of 3 and 1 isolated sampled solution.

This is a finite numerical connectivity result. It does not prove that the corresponding regions are globally disconnected in the continuous 12-rate parameter space. Likewise, local perturbation pass fractions are computational robustness summaries, not biological prevalence estimates.

## 5. Window sensitivity

The operational boundaries 200, 250 and 500 ms were tested by varying the split boundary over 180/200/220 ms, the primary-window endpoint over 230/250/270 ms and the late endpoint over 450/500/550 ms, yielding 27 schemes. The primary representative is retained under every scheme. The exact membership of the calibrated ensemble is quantitatively window-sensitive, whereas experiment-scale amplification remains available throughout the tested definitions.

## 6. Software environment

Core dependencies are Python 3, NumPy, pandas, SciPy, Matplotlib and Numba; raw ABF extraction additionally requires pyABF. The per-step `requirements.txt` files in `code/publication_v1_8/` document the packages used by each final pipeline.
