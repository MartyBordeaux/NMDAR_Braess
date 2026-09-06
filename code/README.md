# Analysis code

The `step_01`-`step_09` directories preserve the versioned analysis pipeline that produced the final manuscript results.

- `step_01`: reconstruct raw ABF inventory and metadata.
- `step_02`: raw-trace QC.
- `step_03`: extract the negative inter-pulse train plateau.
- `step_04`: freeze the experimental plateau target and direction criterion.
- `step_05`: freeze exact legacy model sources and provenance.
- `step_06`: exact plateau rerun of Base, B-slow, and dual extensions.
- `step_07`: B-slow boundary and robustness analysis.
- `step_08`: topology-first full-map analysis and tie-aware post-processing.
- `step_09`: claim-level publication synthesis.

`model/` contains the exact kinetic-network source files used by the final reruns.

The scientific endpoint is the negative inter-pulse plateau during the 25-pulse train. Positive stimulation-artifact peak height is not used as the endpoint.
