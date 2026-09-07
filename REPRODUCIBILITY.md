# Reproducibility

## Source-table reproduction

The primary derived tables permit recalculation of cell/reference plateau ratios and the descriptive 5 potentiating / 5 suppressing / 1 neutral classification. Published model summary tables support inspection of Base maps, constructive examples, B-slow amplification, dual sensitivity and numerical audits. `python code/reproduce_publication.py` checks these tables and rebuilds summary plots; it does not independently regenerate every model trajectory.

## Step11 v1.2 ancillary controls

The exact frozen raw-extraction Python source is deposited in standard XZ compression under `code/step_11/`, with its SHA-256 checked by `run_frozen_v1_2.py`. Decompression exposes the complete original source, not a replacement model or rewritten extraction algorithm. Raw execution requires local ABFs, numpy, pandas, scipy, matplotlib and pyabf.

`python code/step_11/reproduce.py` requires only the saved measurements, numpy, pandas and matplotlib. It combines the two standard archive parts, verifies every member's hash, rechecks condition gates and reported ratios from saved sweep-level metrics, and generates Tables S10-S13 and Figures S4-S5. This source-table path has been executed successfully on the exact staged public files. Raw extraction was not rerun for the v0.9 integration.

Only server path prefixes were normalized in the public sweep tables. No numeric values or flags were changed. The manifest records both original-file and public-member identities. The author manuscript package retains the complete original Step11 results and raw-waveform QC PDF.

The correct negative-voltage count is 11 eligible comparisons in four cell series, from 40 candidates. Fifteen eligible comparisons occurs only over all voltages. At -40 mV two of eight candidates qualify. The ancillary gate was specified after review of v1.1, not preregistered; it was not applied retrospectively to the primary Ro25 cohort. Requiring a negative measurable post-drug plateau can exclude near-complete suppression. All candidates and exclusions are retained.

## Experimental scope

Raw ABF recordings are not distributed. Median plateau values and saved QC measurements cannot independently validate waveform extraction, stimulus masking, baseline subtraction or recording quality. Reused controls and repeated voltages are not independent biological replicates. Missing NBQX labels are undocumented background, not proven absence. PPF signal eligibility does not identify the site of drug action.

The manuscript source package includes baseline-subtracted samples, pointwise median waveforms and interval measurements for its approved Figure 1. Waveform exports and a full raw ABF archive are different levels of reproducibility.

## Deposited Step10 input

`data/model_inputs/step10/accepted_parameter_ensemble.csv.gz` contains 10,000 rows (5,000 broad, 5,000 reference), 17 columns and the downstream parameter fields. Its decompressed CSV SHA-256 is `9a6b6a1388392c6efb96ed3b3e3e199bb737bc5199d7fd4a8135214d16e8acc2`; its compressed Git blob SHA is `c9d5d781f2bd0aee6620a3841f2ddcab15d05139`.

The original archive SHA-256 is `bffddd17332d16c3048b7e23bcb2ab05a3a4c15213de1aa5eba816b23fd59d99`. Its run summary records 50,000 candidates per prior, 5,000 retained per prior, seed 20260728, 0.05-ms integration, 600-ms duration and 25 pulses spaced by 5.1 ms. Calibration used control-shape information, excluding Ro25 outcomes and responder labels.

The accepted table does not reconstruct original proposal distributions, bounds or the exact calibration-score implementation. These cannot be inferred from retained-sample extrema. Obsolete post-train outcome columns are not inputs to the final plateau analysis.

## Model execution scope

The repository contains a compact Base/B-slow implementation and plateau observable, not the complete historical Step01-Step09 batch orchestration or a verified one-command regeneration of all Base/B-slow/dual outputs. Historical hash-identified model sources are additionally supplied in the manuscript package. The deposited Step10 ensemble and the newly tested Step11 source-table workflow do not close this separate model-orchestration gap.

Software is distributed under the MIT License. Fractions of sampled parameters or designed forcing contexts are accessibility measures, not biological prevalence or glutamate-exposure probabilities.
