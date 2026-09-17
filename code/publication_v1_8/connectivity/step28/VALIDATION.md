# Validation notes

The packaged Step 28 code was checked locally before archiving.

Checks performed:

```text
python -m py_compile step28_residual_connectivity_stress.py
python -m py_compile step25_engine_snapshot.py
python step28_residual_connectivity_stress.py --self-test
```

The built-in self-test verifies:

- exact preservation of kinetic endpoints;
- exact preservation of auxiliary endpoints;
- correct construction of 2/3-waypoint piecewise-log-linear paths;
- both `native_aux` and `flex_aux` path families;
- standardized path-length diagnostics.

The production pipeline also contains runtime hard gates for:

- the Step-25 scientific status;
- the Step-27 scientific status;
- the exact 40/40 confirmed strong cohort;
- the exact Step-27 strict residual component freeze;
- strict validity of all inherited Step-27 network edges;
- numerical replay of corrected calibration scores and response ratios;
- recovery of the original full broad-prior box from Step 10.

A production run was not executed locally because the complete Step-25 corrected accepted broad ensemble and the server-side Step-10 full broad table are required. These are auto-discovered on the project server.
