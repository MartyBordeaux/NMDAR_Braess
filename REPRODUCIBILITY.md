# Reproducibility

## Publicly reproducible from this repository

The following can be reproduced directly from the distributed derived/source data:

1. the cell-series plateau ratios and the frozen -40 mV endpoint;
2. the descriptive 5 potentiating / 5 suppressing / 1 neutral classification;
3. manuscript plots based on the distributed source-data tables;
4. the reported Base-topology maps, constructive witnesses, B-slow amplification summaries, dual sensitivity, and numerical-audit tables.

## Raw-data stages

Steps 1-3 operate on the original ABF recordings. The raw ABF archive is intentionally not redistributed in this public repository. The public experimental input is the derived plateau dataset in `data/derived/`, including one row per cell series and negative holding potential plus QC metadata. These values are sufficient to reproduce every experimental numerical statement used in the manuscript after the raw-to-plateau extraction stage.

## Recovered Step10 model input

The previously missing Step10 results archive has been recovered. The archive SHA-256 is:

`bffddd17332d16c3048b7e23bcb2ab05a3a4c15213de1aa5eba816b23fd59d99`

Verified metadata:

- pipeline: `10_matched_parameter_intrareceptor_rerouting`;
- candidate sets: 50,000 per prior;
- accepted sets: 5,000 `broad` + 5,000 `reference`;
- random seed: `20260728`;
- Step10 integration step: 0.05 ms;
- duration: 600 ms;
- pulse count: 25;
- pulse interval: 5.1 ms.

The Step10 run summary records control-only calibration using early/primary, late/primary, and charge/primary control-shape information. It explicitly states that Ro25 outcomes, responder-like labels, separate blocked-condition parameters, and PPF were not used as kinetic calibration constraints.

The final topology-first analyses (Steps 06-09) do not regenerate the original 50,000 candidates. They propagate the frozen accepted ensemble with equal weights and impose the common forcing grid used in the final manuscript.

## Model reruns

The exact downstream kinetic-network sources and plateau-rerun pipelines are preserved under the publication code/provenance layer. The accepted ensemble required by those reruns has now been recovered and verified. A compact 10,000-row publication copy contains `prior`, `sample_id`, the 12 kinetic-rate fields, inherited drive fields, and `calibration_score` — all fields required by the downstream validation and integration code.

The original Step10 archive contains additional legacy endpoint columns (including the obsolete post-train quantities from the earlier analysis). Those columns are not required to reproduce the final inter-pulse plateau calculations and are intentionally not treated as manuscript data.

## Scope

The repository is intended to support reproducibility of the final topology-first manuscript analysis from the frozen derived experimental endpoint and frozen model ensemble. It does not claim de-novo reconstruction of raw electrophysiology without the original ABF archive, nor does it treat the Step10 proposal distribution as a biological population distribution.
