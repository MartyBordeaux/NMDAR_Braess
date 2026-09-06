# Reproducibility

## Source-table reproduction

The public derived tables permit recalculation of cell/reference plateau ratios and the descriptive 5 potentiating / 5 suppressing / 1 neutral classification. Published model summary tables support inspection of Base maps, constructive examples, B-slow amplification, dual sensitivity and numerical audits. `python code/reproduce_publication.py` checks distributed result tables and rebuilds summary plots. It does not independently regenerate every model trajectory behind those tables.

## Experimental scope

Raw ABF recordings are not distributed. Median plateau values and QC fields support downstream comparisons, but cannot independently validate waveform extraction, stimulus masking, baseline subtraction or raw-recording QC. A common same-day reference is not an established within-cell baseline. Reused reference sweeps must not be counted as independent experimental records.

The revised manuscript source package includes the baseline-subtracted samples, pointwise median waveforms and interval measurements used for its Figure 1, plus plotting code. Such waveform exports and a full raw ABF archive are different levels of reproducibility.

## Deposited Step10 input

The accepted input is now present at:

`data/model_inputs/step10/accepted_parameter_ensemble.csv.gz`

It contains 10,000 rows (5,000 broad and 5,000 reference), 17 columns and all parameter fields needed by the downstream kinetics. Its decompressed CSV SHA-256 is `9a6b6a1388392c6efb96ed3b3e3e199bb737bc5199d7fd4a8135214d16e8acc2`. Its compressed Git blob SHA is `c9d5d781f2bd0aee6620a3841f2ddcab15d05139`.

The source archive `results_step10.zip` has SHA-256 `bffddd17332d16c3048b7e23bcb2ab05a3a4c15213de1aa5eba816b23fd59d99`. The recovered run summary records 50,000 candidates per prior, 5,000 retained per prior, seed 20260728, integration step 0.05 ms, duration 600 ms and 25 pulses spaced by 5.1 ms. Calibration used control-shape information and excluded Ro25 outcomes and responder labels.

Availability of the accepted table does not reconstruct the original proposal distributions, proposal bounds or exact calibration-score implementation. Those upstream details cannot be inferred from retained-sample extrema. Obsolete post-train outcome columns in the original Step10 export are not inputs to the final plateau analysis.

## Model execution scope

The repository currently contains a compact Base/B-slow model implementation and a separate plateau observable. It does not contain the complete historical Step01-Step09 batch orchestration or a verified clean-clone, one-command regeneration of all Base, B-slow and dual results. Historical hash-identified model sources are additionally supplied in the manuscript source package; they are not the same thing as the compact public implementation.

The former missing-ensemble limitation is closed. Complete historical workflow deposition and end-to-end rerun verification are separate tasks and are not implied merely by the presence of that input.

## Interpretation

Fractions of accepted parameter sets or designed forcing contexts are model accessibility measures, not biological prevalence. The forcing grid is not a probability distribution of glutamate exposures. Software is distributed under the MIT License.
