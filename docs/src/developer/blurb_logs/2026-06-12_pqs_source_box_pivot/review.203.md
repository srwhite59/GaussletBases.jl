Pass 203 review - accepted

Reviewed the split of `src/CartesianContractedParentMetrics.jl`.

Result:
- Accepted.
- The former 19k-line file is now a 69-line module wrapper.
- The exported metric-packet core is isolated in
  `src/cartesian_contracted_parent_metrics/core.jl`.
- The non-core code is explicitly separated into:
  - `source_box_route_shadow.jl`
  - `current_route_metadata_export.jl`
  - `component_smoke_sidecars.jl`

What changed:
- No public export changes.
- No intended numerical behavior changes.
- No new tests or new route concepts.
- The split exposes that the true contracted-parent metric packet core is much
  smaller than the legacy/private route-shadow, metadata/export, and sidecar
  code around it.

Validation run by manager:
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - passed.
- `julia --project=. -e 'using Test, LinearAlgebra, SparseArrays, GaussletBases; include("test/nested/cartesian_contracted_parent_metric_packet_runtests.jl")'`
  - passed, 29/29.
- `julia --project=. -e 'ENV["TMPDIR"] = "/private/tmp"; include("test/nested/pqs_source_metadata_real_artifact_acceptance_runtests.jl")'`
  - passed the explicit export wrapper, 14/14.
  - artifact opt-in remained skipped/broken as expected without the external artifact directory.
- `julia --project=. test/nested/pqs_component_route_report_adapter_runtests.jl`
  - passed, 605/605.
- `git diff --check`
  - passed.

Line budget:
- The tracked original file loses 18,936 lines and gains 4 wrapper include lines.
- The new include files add 18,931 lines.
- Net source/test/bin line budget: -1.

Deletion/shrinkage assessment:
- deleted: no live helper/test deletion in this pass.
- simplified: `CartesianContractedParentMetrics.jl` is now a small wrapper with
  visible ownership boundaries.
- quarantined: source-box route shadows, current-route metadata/export, and
  component-smoke/CR2 sidecars are no longer mixed with the exported metric
  packet core.
- not deleted because: the split response identified live callers for the
  legacy/private includes; deletion needs a separate stale-caller audit.
- exact remaining caller/blocker: `pqs_projected_q_shell_local_layer_integration_runtests.jl`,
  `pqs_source_metadata_real_artifact_acceptance_support.jl`,
  `pqs_component_route_report_adapter_runtests.jl`,
  `src/cartesian_gto_probes.jl`, `src/pqs_multilayer_shell_source_plan.jl`,
  and route-driver skeleton/helper code.

Manager note:
- This satisfies the first analysis: the coherent metric-packet core is now
  separated from the legacy/private route-shadow and reporting/sidecar code.
- The next simplification should be a caller audit of the legacy includes, not
  another additive seam.

-- repo-manager@macmini
