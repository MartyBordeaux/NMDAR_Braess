# NMDAR_Braess

Publication repository for **Route removal reveals a latent Braess-like amplification regime in an NMDA receptor-state network**.

## Main result

The specified Base receptor-state network permits a Braess-like increase in the open-probability plateau after removal of the B entry route. Forcing and kinetic parameters are fixed across each control/blocked comparison. A rare broad-ensemble realization reaches the observed experimental magnitude range. B-route slowing modifies accessibility and magnitude; it is not necessary for existence. This is a result about the specified kinetic model, not a graph-only theorem for arbitrary rates or proof of the biological action of Ro25-6981.

## Deposited material

- `data/derived/`: primary Ro25 median plateau levels, ratios, frozen targets and voltage robustness.
- `data/derived/step11_v1_2/`: ancillary acute memantine/PPF measurements, eligibility decisions and complete saved-sweep archive.
- `data/model_inputs/step10/accepted_parameter_ensemble.csv.gz`: 10,000 retained sets, 5,000 broad and 5,000 reference.
- `data/model_outputs/`: Base/B-slow/dual summaries and numerical-audit tables.
- `code/model/nmdar_braess_model.py`: compact Base/B-slow implementation.
- `code/plateau_observable.py`: inter-pulse plateau observable.
- `code/reproduce_publication.py`: source-table verification and summary plots.
- `code/step_11/`: exact frozen v1.2 extraction source and tested ancillary publication audit/figure builder.
- `provenance/`: result hierarchy, claim ledger, summaries and ancillary audit notes.

The manuscript is under author proofreading. Its v0.9 update adds ancillary controls without changing the receptor-model results or the approved Figures 1-2.

## Experimental scope

The endpoint is the baseline-subtracted negative current level between stimuli during the 25-pulse train. Artifact neighborhoods are excluded; the first five intervals are discarded and the remaining 19 interval medians are aggregated per sweep before condition aggregation. Artifact peak height is not the endpoint.

The primary Ro25 ratio uses a same-day reference, not an established within-cell baseline. Step11 acute memantine pairs are archive-matched and use additional post-review endpoint gates. At -40 mV two of eight candidates pass. Across planned negative voltages, eleven of forty pairs pass, nested in four cell series; all eleven are suppressive. Fifteen eligible pairs refers to all voltages, including zero/positive voltages. No between-drug population comparison is claimed. PPF cannot exclude presynaptic contributions.

Raw ABF files are not distributed. Saved measurements support downstream aggregation, not independent validation of raw waveform extraction or recording quality control.

## Frozen Step10 input

The source archive SHA-256 is `bffddd17332d16c3048b7e23bcb2ab05a3a4c15213de1aa5eba816b23fd59d99`. Its summary records 50,000 candidates and 5,000 retained sets per prior, seed `20260728`, and control-only shape calibration. Ro25 outcomes and responder labels were not calibration targets.

The deposited input's decompressed CSV SHA-256 is `9a6b6a1388392c6efb96ed3b3e3e199bb737bc5199d7fd4a8135214d16e8acc2`. Its compressed Git blob SHA is `c9d5d781f2bd0aee6620a3841f2ddcab15d05139`. Availability of this table does not reconstruct proposal distributions or the original calibration-score implementation.

## Source-table reproduction

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r environment/requirements.txt
python code/reproduce_publication.py
python code/step_11/reproduce.py
```

The second command verifies all saved ancillary measurements and generates Tables S10-S13, Figures S4-S5 and a machine-readable audit in `reproduced/step11_publication/`. Step11 source and data are standard losslessly compressed files; the wrapper verifies and extracts them automatically. See `code/step_11/README.md` for inspecting the source or rerunning extraction with local ABFs.

These commands do not rerun all historical receptor-model ensembles. Complete historical Step01-Step09 batch orchestration and a verified end-to-end model rerun remain separate tasks. See `REPRODUCIBILITY.md`.

## License

Analysis software is released under the MIT License; see `LICENSE`.
