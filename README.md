# NMDAR_Braess

Publication repository for the study **“Route removal reveals a latent Braess-like amplification regime in an NMDA receptor-state network.”**

## Main result

The unmodified receptor-state network topology itself permits a Braess-like increase in the negative inter-pulse NMDA plateau after complete removal of the B route. Under the broad parameter prior, a rare but numerically stable Base realization also reaches the experimentally observed magnitude range. Route-specific B-path slowing is therefore treated as an **amplifier of an already accessible topological regime**, not as the origin of the possibility.

## Repository layout

- `data/derived/` — public processed electrophysiology data: median inter-pulse plateau values, not raw ABF files.
- `data/model_outputs/` — source tables for the final Base/B-slow/dual analyses.
- `data/model_inputs/step10/` — provenance for the recovered frozen Step10 accepted ensemble used by the final model reruns.
- `code/plateau_observable.py` — frozen plateau observable.
- `code/model/nmdar_braess_model.py` — compact publication implementation of the final kinetic network and B-slow extension.
- `code/reproduce_publication.py` — source-data verification and compact figure regeneration.
- `provenance/` — frozen claim ledger, methods constraints, result hierarchy, and publication summary.
- `REPRODUCIBILITY.md` — exact scope of public reproducibility.

The manuscript is still under author proofreading. The current numerical source data and claim hierarchy are frozen; the final submission TeX/PDF will be synchronized after textual revision is complete.

## Experimental data policy

Raw `.abf` files are not redistributed here. The public experimental dataset consists of the median negative inter-pulse plateau values used in the final analysis.

The primary table is `data/derived/plateau_by_cell_voltage.csv`. It contains one row per cell series and negative holding potential, including the same-day reference plateau, Ro25 plateau, ratio, and QC fields. The frozen -40 mV endpoint is in `data/derived/plateau_minus40_primary.csv`.

## Endpoint

Positive stimulation artifacts are detected only to define pulse timing and are then masked. For each adjacent pulse pair, the median negative inter-pulse current is computed after guard intervals. The first five intervals are discarded, and the plateau is the median of the remaining 19 inter-pulse levels.

The primary magnitude ratio is `r_plateau = |Ro25 plateau| / |same-day reference plateau|`. `r_plateau > 1` is the only binary model direction criterion; magnitude is compared continuously with the five frozen experimental potentiating ratios.

## Recovered Step10 frozen ensemble

The previously missing Step10 archive has been recovered and verified. Its SHA-256 is:

`bffddd17332d16c3048b7e23bcb2ab05a3a4c15213de1aa5eba816b23fd59d99`

The archive contains 50,000 candidate parameter sets per prior and 5,000 accepted sets per prior (`broad` and `reference`). The Step10 run used seed `20260728`. The run summary records a control-only calibration based on early/primary, late/primary, and charge/primary control-shape information; Ro25 outcomes and responder labels were explicitly excluded from calibration.

The final topology-first Steps 06-09 use the frozen 10,000-row accepted ensemble rather than regenerating candidate priors. A column-reduced copy containing every field required by the downstream model code is the publication input. The original Step10 archive also contains legacy outputs that are not used by the final plateau-only analysis.

## Quick reproduction

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r environment/requirements.txt
python code/reproduce_publication.py
```

This verifies the frozen public source tables and regenerates compact diagnostic figures in `reproduced/`.

## License

Analysis software in this repository is released under the MIT License. See `LICENSE`.

## Reproducibility status

See `REPRODUCIBILITY.md`. Raw ABF recordings are intentionally absent from the public repository. The derived experimental endpoint, exact downstream model equations, final source tables, numerical audits, and recovered Step10 ensemble provenance are available for publication-level reproducibility.
