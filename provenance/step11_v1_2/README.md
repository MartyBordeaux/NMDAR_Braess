# Ancillary evidence audit for manuscript v0.9

The frozen result tables were independently checked against all saved sweep
measurements. The table-reproduction command emits `step11_publication_audit.json`
and verifies hashes before use.

Verified counts: 66 candidate voltage pairs; 15 eligible over all voltages;
40 candidates at the five planned negative voltages; 11 eligible negative-voltage
pairs nested in four series (2, 1, 5, 3); eight candidates and two eligible pairs
at -40 mV. All eleven eligible negative-voltage ratios are less than one.

This corrects a previous prose summary that called all 15 eligible comparisons
negative-voltage comparisons. No measurement values were changed.

The two -40 mV ratios are 0.14532864202803603 and 0.7455393925822574.
The 15 March control valid-sweep median is -382.3 pA, distinct from the all-sweep
median -405.8 pA. Excluded records are retained, not treated as no-effect cases.

There are 74 technical PPF sweeps across three pairs. Ro25 cell1 fails the
condition gate (15/42 drug sweeps), Ro25 cell2 has one drug sweep, and the
memantine pair retains 9.17% of the defined mean response amplitude. No
presynaptic exclusion, animal-level replication, or between-drug significance
claim is supported.

The original code's `within-cell acute control` label is retained in the CSV
for provenance. The manuscript qualifies this as archive-based matching,
especially for the unnumbered 29 June series. Exposure timing and undocumented
NBQX backgrounds are not established merely by these labels.
