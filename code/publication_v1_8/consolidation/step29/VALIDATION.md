# Step 29 validation

Package version: `1.0.0`

Validation performed before packaging:

1. `python -m py_compile` passed for the Step29 script and frozen Step25 engine snapshot.
2. Built-in `--self-test` passed.
3. End-to-end smoke test passed against the actual Step24 and Step28 result schemas plus a synthetic Step25-compatible directory.
4. Figure generation passed for the full consolidation path.
5. The local-perturbation branch was exercised with small `n` and successfully ran the frozen Step25 Base-response and corrected-calibration engines at `dt=0.05 ms`.
6. The smoke test selected sample `30902` as `primary_typical_dominant_native`, consistent with the frozen 40-vector input table.
7. Output archiving completed successfully.

Production hard gates prevent interpretation if upstream statuses, strong count, component sizes, thresholds, witness membership, or the Step28 replay gate differ from the frozen state.

The smoke-test numerical fractions are not scientific results and are not distributed in the production package.
