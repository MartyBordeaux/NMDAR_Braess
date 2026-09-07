# Frozen Step11 v1.2: ancillary experimental controls

This is supporting evidence for manuscript v0.9, not a new model fit.

## Reproduce the publication supplement from saved measurements

From the repository root:

```bash
python -m pip install numpy pandas matplotlib
python code/step_11/reproduce.py
```

The command verifies the losslessly compressed source and data, recomputes
cell-condition eligibility and every reported ratio from saved sweep metrics,
and generates Tables S10-S13, Figures S4-S5, two negative-voltage CSV views and
a machine-readable audit in `reproduced/step11_publication/`.

## Inspect or rerun the exact extraction code

`step11_experimental_controls_v1_2.py.xz` is the original Python file in standard
XZ compression, not a new implementation. To inspect it, decompress with
`xz -dc code/step_11/step11_experimental_controls_v1_2.py.xz` or Python's `lzma`.
The original source SHA-256 is
`b994c232f779d1a2a04940530c359c9be17a9895b35617b03ddd2cfc34476419`.

```bash
python -m pip install numpy pandas scipy matplotlib pyabf
python code/step_11/run_frozen_v1_2.py --self-test
python code/step_11/run_frozen_v1_2.py --raw-root ~/nmda/IV_NMDA --output reproduced/step11_raw_v1_2
```

Raw ABFs are not distributed. The self-test does not validate real recordings.
`build_publication.py.xz` is the independently run source-table audit/figure
builder; its hash is checked by `reproduce.py`.

## Interpretation

Eight acute memantine pairs are candidates at -40 mV; only two pass the additional
endpoint gates. Eleven of forty planned negative-voltage pairs pass, nested in
four cell series, and all eleven are suppressive. Fifteen eligible rows occur
only across the entire voltage output, including zero/positive voltages.

The gate was specified after v1.1 inspection. It requires measurable negative
plateaus in both conditions and can exclude near-complete blockade. It was not
used to reselect the primary Ro25 cohort. No new p-values are produced.

One Ro25 PPF pair fails the drug-condition SNR gate (15/42 valid sweeps); the
other has one post-drug sweep. The memantine PPF pair has marked response loss.
PPR does not exclude presynaptic contributions. Dates/cell labels are archive
identifiers, not animal identifiers. Undocumented NBQX is not proven absence.
