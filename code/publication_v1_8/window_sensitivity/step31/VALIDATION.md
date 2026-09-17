# Validation

Package version: 1.0.0

Static checks performed before packaging:

- Python syntax compilation passed.
- Built-in numerical self-test passed.
- Experimental baseline anchoring is checked independently for all five frozen calibration series.
- Raw calibration mapping has explicit hard gates against Ro25 and memantine condition paths.
- SlowEPSC protocol identity is required.
- Baseline model-score and accepted-set replay gates are retained.
- No new response threshold, denominator cutoff, acceptance rule, or Bayesian posterior is introduced.

The full raw-ABF run must be executed on the project server because the ABF archive and the full Step-10/24 inputs are not redistributed in this package.
