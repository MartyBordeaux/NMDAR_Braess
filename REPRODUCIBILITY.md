# Reproducibility

## Publicly reproducible from this repository

The following can be reproduced directly from the distributed derived/source data:

1. the cell-series plateau ratios and the frozen -40 mV endpoint;
2. the descriptive 5 potentiating / 5 suppressing / 1 neutral classification;
3. manuscript plots based on the distributed source-data tables;
4. the reported Base-topology maps, constructive witnesses, B-slow amplification summaries, dual sensitivity, and numerical-audit tables.

## Raw-data stages

Steps 1-3 operate on the original ABF recordings. The scripts are included so that the extraction procedure is inspectable, but the raw ABF archive is not redistributed in this repository. The public derived experimental input is `data/derived/plateau_by_cell_voltage.csv`.

## Model rerun input

The exact model equations and rerun pipelines are included under `code/model/` and `code/step_05`-`code/step_09`.

The frozen Step10 accepted parameter ensemble is an upstream input to the final model reruns. It is not present in the current publication export. Until that ensemble is added, the repository reproduces the published model outputs and figures from source tables but does not yet provide a complete one-command regeneration of every trajectory from the original parameter sampling stage.

This limitation is explicit so that the public repository does not overstate reproducibility.
