# Step 31 results — final window robustness and parameter-space geometry

Status: **WINDOW_ROBUSTNESS_AND_GEOMETRY_COMPLETED**

The post-train start was held fixed at 140 ms. The operational boundaries were varied on a prespecified 3x3x3 grid: split = 180/200/220 ms, primary-window end = 230/250/270 ms, late/integral end = 450/500/550 ms. For every scheme, untreated SlowEPSC reference recordings were re-extracted from raw ABF data. The exact five frozen series-level baseline ratios were retained as anchors, while raw ABFs supplied only the within-series changes caused by moving the window boundaries. Experimental targets and floor-bounded robust tolerances were then recomputed, and all 100,000 Step-10 candidates were recalibrated. Model-window sensitivity uses a paired-ratio anchor: exact frozen Step-24 baseline observables are retained, while the prefix integrator supplies only the within-candidate change caused by moving a window boundary.

The baseline scheme 140-200 / 140-250 / 200-500 / 140-500 was required to replay all five frozen series-level ratios after anchoring and the frozen candidate calibration scores before any sensitivity result was accepted. Ro25 and memantine condition files are hard-excluded from calibration mapping.

Headline sensitivity quantities:

- minimum broad-cohort Jaccard versus baseline: **0.651**
- median broad-cohort Jaccard: **0.830**
- minimum retained count among the frozen 40 experiment-scale solutions: **20/40**
- minimum retained count from the dominant 36-solution network: **19/36**
- primary typical representative retained under every scheme: **True**

`Fig31B_parameter_space_geometry` is the manuscript-facing geometry figure. Panel A shows all 5,000 control-compatible broad-prior models in a two-dimensional PCA projection of standardized log-rates, with the 40 experiment-scale solutions overlaid. Panel B adds only verified connectivity edges. Panel C shows calibration score against amplification magnitude for the 40 solutions. Panel D explicitly defines what a verified path means by plotting calibration score and amplification ratio along one path; this is a path through parameter space, not a receptor-state trajectory.

`Fig31C_rate_filtering_marginals` shows sampled broad, control-compatible, and experiment-scale marginal rate distributions. These are filtering distributions from the search/calibration procedure and are **not** interpreted as Bayesian priors or posteriors.
