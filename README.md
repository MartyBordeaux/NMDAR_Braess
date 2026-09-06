# NMDAR_Braess

Publication repository for **Route removal reveals a latent Braess-like amplification regime in an NMDA receptor-state network**.

## Main result

The specified Base receptor-state network permits a Braess-like increase in the open-probability plateau after removal of the B entry route. The same forcing and kinetic parameters are held fixed across each control/blocked comparison. Under the broad accepted ensemble, a rare numerically stable realization also reaches the observed experimental magnitude range. B-route slowing modifies accessibility and magnitude; it is not necessary for existence of the Base effect. This is a result about the specified kinetic model, not a graph-only theorem for arbitrary rates or proof of the biological action of Ro25-6981.

## Files actually provided

- `data/derived/`: median experimental inter-pulse plateau levels, cell/reference ratios, frozen targets and voltage robustness.
- `data/model_inputs/step10/accepted_parameter_ensemble.csv.gz`: the deposited 10,000-row, 17-column accepted input (5,000 broad and 5,000 reference rows).
- `data/model_inputs/step10/`: input metadata and Step10 provenance.
- `data/model_outputs/`: published Base/B-slow/dual summaries and numerical-audit tables.
- `code/model/nmdar_braess_model.py`: compact Base/B-slow implementation.
- `code/plateau_observable.py`: inter-pulse plateau observable.
- `code/reproduce_publication.py`: verification of distributed result tables and regeneration of summary plots.
- `provenance/`: result hierarchy, claim ledger and summary tables.

The manuscript remains under author proofreading. The repository does not currently contain the complete historical Step01-Step09 batch workflow or a tested clean-clone regeneration of every trajectory. The manuscript source package contains additional historical model sources; these should not be confused with the compact public implementation.

## Endpoint and experimental scope

The experimental endpoint is the baseline-subtracted negative current level between stimuli during the 25-pulse train. Artifact neighborhoods are excluded, the first five intervals are discarded, and the median of the remaining 19 interval medians is calculated per sweep before aggregation across sweeps. Artifact peak heights are not the endpoint.

The experimental ratio is `abs(Ro25 plateau) / abs(same-day reference plateau)`. The reference is not established as a within-cell baseline. Counts of cell series and fractions of sampled model parameters are not biological prevalence estimates.

Raw ABF files are intentionally absent. Derived tables permit downstream numerical comparisons but do not independently reproduce raw-waveform extraction or raw-recording quality control.

## Frozen model input

The original `results_step10.zip` SHA-256 is:

`bffddd17332d16c3048b7e23bcb2ab05a3a4c15213de1aa5eba816b23fd59d99`

Its run summary records 50,000 candidates and 5,000 retained sets per prior, seed `20260728`, and control-only shape calibration. Ro25 outcomes and responder labels were not calibration targets. The deposited accepted table is the downstream input, not a reconstruction of the original proposal distributions or exact score implementation.

The deposited compressed input is present now, not pending release. Its Git blob SHA is `c9d5d781f2bd0aee6620a3841f2ddcab15d05139`; the SHA-256 of its decompressed CSV is `9a6b6a1388392c6efb96ed3b3e3e199bb737bc5199d7fd4a8135214d16e8acc2`.

## Quick source-table reproduction

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r environment/requirements.txt
python code/reproduce_publication.py
```

This command verifies distributed summaries and regenerates diagnostic plots in `reproduced/`. It does not rerun all original model ensembles. See `REPRODUCIBILITY.md` for the distinction between data availability, source-table reproduction and complete computational reproduction.

## License

Analysis software is released under the MIT License; see `LICENSE`.
