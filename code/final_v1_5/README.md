# Final v1.5 analysis source snapshots

This directory is the publication-aligned source/provenance area for manuscript v1.5.

Deposited source snapshots include the historical Step10 proposal/calibration implementation and the final geometry/compatibility/Bayesian/consolidation analyses used for Steps18-23. Large historical snapshots may be gzip-compressed; they are byte-preserving standard gzip streams and can be inspected with `gzip -dc`.

`verify_publication_v1_5.py` is a fast integrity audit over the compact final tables. It checks the current 5/5/1 experimental classification, exact Step20 calibration replay, Step21 local-basin fraction, Step22 primary occupancy estimate and Step23 frozen-threshold global counts. It does not rerun the expensive global searches.

The source snapshots preserve their original server-oriented discovery logic and historical conventions. Re-running the heavy analyses on a new machine requires restoring the frozen upstream inputs and adapting paths where necessary. Raw ABF extraction remains outside the public-repository scope.
