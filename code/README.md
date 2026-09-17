# Analysis code

The publication branch is organized around the final analysis used for manuscript v1.8.

## Final publication pipelines

The directory `code/publication_v1_8/` contains the final computational chain used for the manuscript:

- control-shape calibration and frozen candidate scoring;
- confirmation of experiment-scale amplification on the fixed 80-node forcing grid;
- sampled connectivity searches (Steps 26-28);
- publication consolidation and representative/local-robustness analysis (Step 29);
- post-train analysis-window sensitivity with raw-ABF replay (Step 31).

The compact Base receptor-state implementation remains under `model/`. Earlier code under `final_v1_5/` and `step_11/` is retained only for provenance and upstream experimental utilities.

## Experimental raw-data processing

Raw ABF files are deposited under `data/raw/`. Final raw-file mapping and the post-train experimental metric implementation used for the window-sensitivity audit are contained in the Step 31 pipeline. The primary sweep-separation endpoint and voltage-robustness tables are distributed under `data/source_data_v1_8/`.

## Interpretation

The parameter candidate distributions and calibration-filtered ensembles are computational samples. They are not Bayesian priors/posteriors. Connectivity is defined by explicitly tested paths in the full kinetic parameter space, not by distances in a two-dimensional projection.
