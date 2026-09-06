# NMDAR_Braess

Publication repository for the study **“Route removal reveals a latent Braess-like amplification regime in an NMDA receptor-state network.”**

## Main result

The unmodified receptor-state network topology itself permits a Braess-like increase in the negative inter-pulse NMDA plateau after complete removal of the B route. Under the broad parameter prior, a rare but numerically stable Base realization also reaches the experimentally observed magnitude range. Route-specific B-path slowing is therefore treated as an **amplifier of an already accessible topological regime**, not as the origin of the possibility.

## Repository layout

- `manuscript/` — current TeX manuscript, Supplementary Information, and bibliography.
- `data/derived/` — public processed electrophysiology data: median inter-pulse plateau values, not raw ABF files.
- `data/model_outputs/` — source tables for the final Base/B-slow/dual analyses.
- `code/step_01` ... `code/step_09` — versioned analysis pipeline.
- `code/model/` — exact kinetic-network source files used by the final reruns.
- `provenance/` — frozen claim ledger, methods constraints, and numerical provenance.
- `REPRODUCIBILITY.md` — what can and cannot currently be regenerated from the public repository alone.

## Experimental data policy

Raw `.abf` files are not redistributed here. The public experimental dataset consists of the median negative inter-pulse plateau values used in the final analysis.

The primary table is `data/derived/plateau_by_cell_voltage.csv`. It contains one row per cell series and negative holding potential, including the same-day reference plateau, Ro25 plateau, ratio, and QC fields.

## Endpoint

Positive stimulation artifacts are detected only to define pulse timing and are then masked. For each adjacent pulse pair, the median negative inter-pulse current is computed after guard intervals. The first five intervals are discarded, and the plateau is the median of the remaining 19 inter-pulse levels.

The primary magnitude ratio is `r_plateau = |Ro25 plateau| / |same-day reference plateau|`. `r_plateau > 1` is the only binary model direction criterion; magnitude is compared continuously with the five frozen experimental potentiating ratios.

## Environment

Python dependencies are listed in `environment/requirements.txt`.

## Reproducibility status

See `REPRODUCIBILITY.md`. Raw ABF files are intentionally absent. The frozen Step10 accepted model-parameter ensemble still needs to be added for a complete de novo regeneration of every model trajectory. All final source tables used for the manuscript are present.
