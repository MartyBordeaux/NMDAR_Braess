# Frozen Step10 accepted ensemble

The exact Step10 results archive has been recovered and verified.

- Source archive SHA-256: `bffddd17332d16c3048b7e23bcb2ab05a3a4c15213de1aa5eba816b23fd59d99`
- Accepted rows: 10,000 total = 5,000 `broad` + 5,000 `reference`
- Compact publication copy SHA-256: `8abb26655348b8f7c5e48d1a51610ac3853b50e21e7d37e37889a235cec57613`
- Compact publication filename: `step10_accepted_parameter_ensemble.csv`

The compact copy contains all columns required by the downstream topology-first rerun code: `prior`, `sample_id`, the 12 kinetic-rate fields, inherited pulse-amplitude and glutamate-decay fields, and `calibration_score`.

The full source table also contains legacy post-train endpoint outputs from an earlier analysis. Those columns are not used by the final inter-pulse plateau analysis and are intentionally excluded from the compact publication input.

The binary/compressed accepted-ensemble file is included in the manuscript v0.6 source package and is intended to be attached to the immutable repository release at submission. The GitHub connector used to populate this repository cannot directly transfer that local binary asset, so this file records its exact identity and provenance until release deposition.
