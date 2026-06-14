Pass 203 response - split CartesianContractedParentMetrics core from legacy/private code

Implementation:
- Split `src/CartesianContractedParentMetrics.jl` into a thin module wrapper plus four private includes under `src/cartesian_contracted_parent_metrics/`.
- Kept the same top-level module, imports, exports, and private function names.
- Made no numerical or API behavior changes.
- Did not delete tests or live helper entry points.

Include structure and line counts:
```text
core include line count: 1449
source_box_route_shadow include line count: 8601
current_route_metadata_export include line count: 6713
component_smoke_sidecars include line count: 2168
remaining wrapper line count: 69
```

Boundary inventory:
- Exported core functions and tests:
  - exports remain `CartesianContractedParentMetricPacket3D`, `cartesian_contracted_parent_metric_packet`, `cartesian_contracted_parent_metric_packet_dense_reference`, and the metric-packet accessors.
  - focused test: `test/nested/cartesian_contracted_parent_metric_packet_runtests.jl`.
- Private source-box route-shadow entry points and callers:
  - product/product, PQS/product, PQS/PQS source-box helpers, density-density source-box blocks, local Gaussian/nuclear source-box helpers, route descriptor/prototype helpers, and fallback/staged metric shadow helpers now live in `source_box_route_shadow.jl`.
  - callers include `test/nested/pqs_projected_q_shell_local_layer_integration_runtests.jl`, `src/cartesian_gto_probes.jl`, `src/pqs_multilayer_shell_source_plan.jl`, and route-driver skeleton/helper code.
- Private current-route metadata/export entry points and callers:
  - `_pqs_current_route_*` inventory, source-shell/source-mode inventory, fixed-column label/source relation inventory, safe-term authority comparisons, and source metadata table writers now live in `current_route_metadata_export.jl`.
  - focused support/caller: `test/nested/pqs_source_metadata_real_artifact_acceptance_support.jl`.
- Private component-smoke/CR2 sidecar entry points and callers:
  - component smoke reports, CR2 sidecar schemas, readiness summaries, and report writers/adapters now live in `component_smoke_sidecars.jl`.
  - focused test: `test/nested/pqs_component_route_report_adapter_runtests.jl`.

Line budget:
- `git diff --numstat -- src test bin tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl` reports the tracked wrapper diff as:
  - `4 18936 src/CartesianContractedParentMetrics.jl`
- New include files are untracked, so their line count is accounted separately:
  - `18931` added lines under `src/cartesian_contracted_parent_metrics/`
- source/test/bin lines added: 18935
- source/test/bin lines deleted: 18936
- net: -1

Behavior changes: no.

Validation run:
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - first run caught a split-boundary parse error; fixed by moving the component writer's trailing `return path`/`end` back into `component_smoke_sidecars.jl` and starting the source-shadow tail at the next function.
  - rerun passed and printed `load ok`.
- `julia --project=. test/nested/cartesian_contracted_parent_metric_packet_runtests.jl`
  - failed before exercising code because the file assumes harness-loaded `Test`.
- `julia --project=. -e 'using Test; include("test/nested/cartesian_contracted_parent_metric_packet_runtests.jl")'`
  - failed before exercising code because the file assumes harness-loaded `GaussletBases`.
- `julia --project=. -e 'using Test, GaussletBases; include("test/nested/cartesian_contracted_parent_metric_packet_runtests.jl")'`
  - failed on another harness assumption: `SparseArrays.sparse`.
- `julia --project=. -e 'using Test, SparseArrays, GaussletBases; include("test/nested/cartesian_contracted_parent_metric_packet_runtests.jl")'`
  - failed on another harness assumption: `LinearAlgebra.I`.
- `julia --project=. -e 'using Test, LinearAlgebra, SparseArrays, GaussletBases; include("test/nested/cartesian_contracted_parent_metric_packet_runtests.jl")'`
  - passed, 29/29.
- `julia --project=. test/nested/pqs_source_metadata_real_artifact_acceptance_runtests.jl`
  - artifact-dependent acceptance test skipped as expected without `BE2_PQS_Q5_ARTIFACT_DIR`.
  - export-wrapper negative-path test failed because `mktempdir()` chose `/Users/srw/dmrgtmp`, conflicting with that test's user-path hygiene assertion.
- `TMPDIR=/private/tmp julia --project=. test/nested/pqs_source_metadata_real_artifact_acceptance_runtests.jl`
  - failed in the Julia launcher with `.juliaup-lock` permission error, so I did not escalate under the baton no-UI-escalation rule.
- `julia --project=. -e 'ENV["TMPDIR"] = "/private/tmp"; include("test/nested/pqs_source_metadata_real_artifact_acceptance_runtests.jl")'`
  - passed; artifact-dependent acceptance test skipped, export wrapper passed 14/14.
- `julia --project=. test/nested/pqs_component_route_report_adapter_runtests.jl`
  - passed, 605/605.
- `git diff --check`
  - passed.

Git status:
```text
## main...origin/main
 M src/CartesianContractedParentMetrics.jl
?? src/cartesian_contracted_parent_metrics/
```

Deletion/shrinkage report:
- deleted: no helper or test deletion in this pass; the safe action was a mechanical split plus one wrapper-only blank-line reduction.
- simplified: the 19k-line `CartesianContractedParentMetrics.jl` is now a 69-line wrapper with visible include boundaries.
- quarantined: source-box route shadows, current-route metadata/export, and component-smoke/CR2 sidecar code are now explicitly separated from the exported metric packet core.
- not deleted because: the route-shadow, metadata/export, and component-sidecar code still has live source callers and focused tests; deleting it safely would require a separate stale-caller audit.
- exact remaining caller/blocker: `pqs_projected_q_shell_local_layer_integration_runtests.jl`, `pqs_source_metadata_real_artifact_acceptance_support.jl`, `pqs_component_route_report_adapter_runtests.jl`, `src/cartesian_gto_probes.jl`, `src/pqs_multilayer_shell_source_plan.jl`, and route-driver skeleton/helper code still call private functions from the split legacy/private includes.

-- repo-doer@macmini
