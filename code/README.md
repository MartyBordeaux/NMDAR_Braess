# Analysis code

The public code is organized around the frozen scientific endpoint and the final topology-first analysis rather than the internal development-step numbering.

- `plateau_observable.py` — exact frozen inter-pulse plateau observable used for the final model comparison.
- `model/nmdar_braess_model.py` — compact publication implementation of the Base network, complete B-route removal, normalized 25-pulse forcing, and the B-slow route-split extension.
- `reproduce_publication.py` — verifies the frozen endpoint/topology-first summary tables and regenerates compact diagnostic figures from the distributed source data.

The original development workflow proceeded through Steps 01–09 (raw reconstruction, QC, plateau extraction, target freeze, source freeze, exact reruns, B-slow boundary analysis, topology-first map, and publication synthesis). Those internal versioned scripts are retained in the project archive, but the public repository exposes the compact publication implementation and the final source tables rather than requiring the historical server directory layout.

The scientific endpoint is the negative inter-pulse plateau during the 25-pulse train. Positive stimulation-artifact peak height is not used as the endpoint.

## Raw ABF extraction

Raw ABF files are not distributed. Consequently, the public repository starts experimental reproduction from `data/derived/plateau_by_cell_voltage.csv`. The raw-data extraction algorithm is documented in the manuscript and is represented by the same interval definition implemented in `plateau_observable.py`.

## Model reruns

A complete rerun of the 5000+5000 model ensemble additionally requires the frozen Step10 accepted parameter table. That table was not present in the export used to initialize this repository; see `../REPRODUCIBILITY.md`.
