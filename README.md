# NMDAR_Braess

Publication repository for the study **“Route removal reveals a latent Braess-like amplification regime in an NMDA receptor-state network.”**

## Main result

The unmodified receptor-state network topology itself permits a Braess-like increase in the negative inter-pulse NMDA plateau after complete removal of the B route. Under the broad parameter prior, a rare but numerically stable Base realization also reaches the experimentally observed magnitude range. Route-specific B-path slowing is therefore treated as an **amplifier of an already accessible topological regime**, not as the origin of the possibility.

## Repository layout

- `data/derived/` — public processed electrophysiology data: median inter-pulse plateau values, not raw ABF files.
- `data/model_outputs/` — source tables for the final Base/B-slow/dual analyses.
- `code/plateau_observable.py` — frozen plateau observable.
- `code/model/nmdar_braess_model.py` — compact publication implementation of the final kinetic network and B-slow extension.
- `code/reproduce_publication.py` — source-data verification and compact figure regeneration.
- `provenance/` — frozen claim ledger, methods constraints, result hierarchy, and publication summary.
- `REPRODUCIBILITY.md` — exact scope and current limits of public reproducibility.

The manuscript is still under author proofreading and is therefore not yet treated as a frozen repository artifact. Its final TeX/PDF version will be synchronized after textual revision is complete; the numerical source data and claim hierarchy are already frozen here.

## Experimental data policy

Raw `.abf` files are not redistributed here. The public experimental dataset consists of the median negative inter-pulse plateau values used in the final analysis.

The primary table is `data/derived/plateau_by_cell_voltage.csv`. It contains one row per cell series and negative holding potential, including the same-day reference plateau, Ro25 plateau, ratio, and QC fields. The frozen -40 mV endpoint is in `data/derived/plateau_minus40_primary.csv`.

## Endpoint

Positive stimulation artifacts are detected only to define pulse timing and are then masked. For each adjacent pulse pair, the median negative inter-pulse current is computed after guard intervals. The first five intervals are discarded, and the plateau is the median of the remaining 19 inter-pulse levels.

The primary magnitude ratio is `r_plateau = |Ro25 plateau| / |same-day reference plateau|`. `r_plateau > 1` is the only binary model direction criterion; magnitude is compared continuously with the five frozen experimental potentiating ratios.

## Quick reproduction

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r environment/requirements.txt
python code/reproduce_publication.py
```

This verifies the frozen public source tables and regenerates compact diagnostic figures in `reproduced/`.

## Reproducibility status

See `REPRODUCIBILITY.md`. Raw ABF files are intentionally absent. The frozen Step10 accepted model-parameter ensemble is also not present in the current export, so every 5000-set trajectory cannot yet be regenerated de novo from public files alone. The model equations, endpoint code, final source tables, numerical audits, and manuscript-level result provenance are included.
