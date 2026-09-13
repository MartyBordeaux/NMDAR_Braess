# Final v1.5 analysis provenance

This directory is the publication-aligned source/provenance area for manuscript v1.5.

The repository currently deposits byte-preserving source snapshots for the historical Step10 proposal/calibration implementation and Steps18-19, together with `verify_publication_v1_5.py`. Compact authoritative outputs for Steps20-23 are deposited under `data/derived/final/step20` ... `step23`; the complete Step20-23 historical source files remain part of the author analysis package and should be added to a later archival release if a full heavy rerun is required.

`verify_publication_v1_5.py` is a fast integrity audit over the compact final tables. It checks the current 5/5/1 experimental classification, exact Step20 calibration replay, Step21 local-basin fraction, Step22 primary occupancy estimate and Step23 frozen-threshold global counts. It does not rerun the expensive global searches.

The deposited historical source snapshots preserve their original server-oriented discovery logic and conventions. Re-running the heavy analyses on a new machine requires restoring the frozen upstream inputs and adapting paths where necessary. Raw ABF extraction remains outside the public-repository scope.
