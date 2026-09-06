# Frozen Step10 accepted ensemble

The compact accepted ensemble is deposited as `accepted_parameter_ensemble.csv.gz` in this directory. It is not awaiting upload or release deposition.

- Rows: 10,000 = 5,000 broad + 5,000 reference.
- Columns: 17, retaining prior, sample_id, 12 kinetic rates, two inherited forcing fields and calibration_score.
- Source archive SHA-256: `bffddd17332d16c3048b7e23bcb2ab05a3a4c15213de1aa5eba816b23fd59d99`.
- Deposited gzip SHA-256: `596723b5ce4959a42422f9862ea6561c20ed2fc5ee4e31bc6fd99d9d82cbc6cf`.
- Decompressed CSV SHA-256: `9a6b6a1388392c6efb96ed3b3e3e199bb737bc5199d7fd4a8135214d16e8acc2`.
- Git blob SHA-1: `c9d5d781f2bd0aee6620a3841f2ddcab15d05139`.

Row identifiers match the original accepted table. Comparing the retained numeric columns after parsing the original and compact CSVs gives a maximum absolute difference of approximately 7.11e-15 from serialization roundoff. The compressed public copy must therefore be checked against its own byte hashes, not against an earlier serialization.

The original Step10 output included obsolete post-train outcomes. Those extra columns are intentionally not used as final plateau-analysis targets. The run summary documents control-only calibration, excluding Ro25 outcomes and responder labels. It does not recover the complete upstream proposal distributions or the exact calibration-score implementation.

This file closes the missing accepted-input gap. It does not, by itself, certify that the public compact model and quick summary script reproduce the complete historical batch workflow; see the repository REPRODUCIBILITY.md.
