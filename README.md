# NMDAR Braess-like amplification

Publication repository for the Neuropharmacology manuscript **“Paradoxical NMDA response amplification in a robust receptor-state regime.”**

The study combines archival Purkinje-cell electrophysiology with a minimal NMDAR-inspired receptor-state model. Ro25-6981-associated recordings contain both potentiating and suppressing inter-pulse responses. In the model, making one receptor-state entry route inaccessible can increase open-state occupancy while all other kinetic rates, stimulation and initial conditions remain unchanged. The final analysis asks whether experiment-scale amplification exists within models that independently reproduce the control-response shape, whether qualifying solutions can be connected through other qualifying kinetic models, and whether representative amplification survives local parameter and analysis-window perturbations.

## Final publication structure

- `code/publication_v1_8/` — final model calibration, response confirmation, geometry, connectivity, consolidation and window-sensitivity code.
- `code/model/`, `code/step_11/`, `code/final_v1_5/` — retained upstream model and experimental utilities needed for provenance and replay.
- `data/source_data_v1_8/` — final machine-readable source tables used in the manuscript.
- `data/model_inputs/step10_full/results_step10.zip` — frozen full Step 10 candidate table used by the downstream publication analysis.
- `data/raw/IV_NMDA/` — original ABF recordings, tracked with Git LFS while preserving the experimental archive structure.
- `data/raw/SHA256SUMS.txt` and `data/raw/ABF_FILELIST.txt` — integrity and inventory manifests for the raw archive.
- `REPRODUCIBILITY.md` — execution order, numerical definitions and interpretation limits.

## Headline results

At -40 mV, complete sweep separation classifies the 11 Ro25 series as 5 potentiating, 5 suppressing and 1 unresolved. Among the control-calibrated model ensemble, 40 solutions reach at least the smallest observed potentiating ratio. Finite parameter-path searches place 36 of these solutions in one dominant connected component, three in a smaller component and one as an isolated sampled solution under the final endpoint-preserving criterion. The primary representative has an amplification ratio of approximately 1.87 and remains compatible under all 27 predefined post-train window schemes.

## Repository status

This repository contains the final manuscript v1.8 computational code, machine-readable source data, the frozen Step 10 model-input archive and the raw electrophysiological ABF archive. Raw ABFs are tracked with Git LFS; users cloning the repository for raw-trace replay should have Git LFS installed and run `git lfs pull` if necessary.

The publication analysis is frozen to the definitions documented in `PUBLICATION_VERSION.md` and `REPRODUCIBILITY.md`. Historical code and earlier analysis products are retained only for provenance and should not supersede the v1.8 publication chain.

## Reuse and citation

Analysis software is released under the MIT License. See `PUBLICATION_VERSION.md` and `REPRODUCIBILITY.md` for the frozen publication definitions and scope.
