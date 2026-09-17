# Validation record

The bundled code was checked before packaging with:

```bash
python -m py_compile step26_strong_connectivity.py step25_engine_snapshot.py
python step26_strong_connectivity.py --self-test
```

The self-test passed.

A real-input preflight was run against the uploaded Step-25 result bundle. It verified:

- Step-25 status `KINETIC_REGION_STRUCTURE_UNRESOLVED`;
- 40 confirmed strong vectors;
- 5000 corrected accepted broad vectors;
- corrected midpoint and strict cutoffs;
- exact corrected target file;
- bundled Step-25 engine replay against representative Step-25 calibration scores and `dt=0.05 ms` response ratios.

A reduced end-to-end smoke run was also executed on the real Step-25 outputs using 39 coarse edges and three path points per edge. It completed successfully through coarse screening, edge refinement, graph-component construction, decision output and figure generation.

The reduced smoke result is not a scientific result; the server run must use the default 21/41-point settings and the intended edge budget.
