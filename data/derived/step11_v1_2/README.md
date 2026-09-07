# Step11 v1.2 derived measurements

The native small CSV files expose all eight -40 mV eligibility decisions, the
two admitted ratios, negative-voltage robustness, and PPF condition/pair results.

The complete sweep- and cell-level input is stored losslessly in
`tables.tar.xz.001` and `tables.tar.xz.002`. Concatenate in that order to recover
a standard XZ-compressed tar archive. Its SHA-256 is
`7ed39e3e12dfd9c4c0455cef8a705956ba1a11f9f92fa7ed44525cfa9d846303`.
`python code/step_11/reproduce.py` does this automatically, verifies every member
against `PUBLIC_DATA_MANIFEST.json`, and rebuilds the supplement.

Only absolute server path prefixes were normalized to `IV_NMDA/` in the public
sweep tables. No numeric entries or eligibility flags were modified. Original
result-file hashes are recorded in the manifest; the author source package
retains the original complete server results and waveform-QC PDF.

There are 66 candidate voltage pairs overall, 40 at the five planned negative
voltages, 11 eligible negative-voltage pairs in four cell series, and two
eligible pairs out of eight at -40 mV. Fifteen is the eligible count across all
voltages, not the negative-voltage count. These are not independent animal
replicates and do not support a between-drug population comparison.

The saved data reproduce downstream aggregation/QC but cannot independently
validate stimulus detection or raw waveform quality without the ABFs.
