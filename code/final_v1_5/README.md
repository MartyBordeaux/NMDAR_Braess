# Final v1.5 analysis source snapshots

This directory contains frozen source snapshots aligned to manuscript v1.5. Large source files are stored as standard gzip streams to keep the repository diff compact. Inspect with, for example:

```bash
gzip -dc code/final_v1_5/historical_step10/10_matched_parameter_rerouting.py.gz | less
gzip -dc code/final_v1_5/step18/step18_prior_geometry.py.gz | less
```

The historical Step10 source defines the original Latin-hypercube proposal and exact control-calibration score. Step18 is the global prior-geometry analysis. Additional final analysis archives/source snapshots for Steps19-23 are deposited alongside this directory as the publication branch is consolidated.

These snapshots preserve the original server-oriented discovery logic and historical conventions. They should be read together with `REPRODUCIBILITY.md` and `PUBLICATION_VERSION.md`; raw ABF extraction is outside the public-repository scope.
