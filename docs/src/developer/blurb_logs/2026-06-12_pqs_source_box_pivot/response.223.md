Pass 223 response - stale component-smoke sidecar retirement

Caller audit before deletion:
- `src/CartesianContractedParentMetrics.jl`
  - included `cartesian_contracted_parent_metrics/component_smoke_sidecars.jl`.
- `test/nested/integration_runtests.jl`
  - included `pqs_component_route_report_adapter_runtests.jl`.
- `test/nested/pqs_component_route_report_adapter_runtests.jl`
  - dedicated stale report-adapter/CR2-sidecar test for the component-smoke sidecar helpers.
- `test/nested/pqs_projected_q_shell_local_layer_integration_runtests.jl`
  - one helper block `_check_pqs_pqs_product_source_box_component_route_smoke(...)`.
  - one call to that helper from the broader integration fixture.
- No non-test source caller or accepted endpoint contract referenced this surface.

Deleted/changed:
- Removed the include from `src/CartesianContractedParentMetrics.jl`.
- Deleted `src/cartesian_contracted_parent_metrics/component_smoke_sidecars.jl`.
- Deleted `test/nested/pqs_component_route_report_adapter_runtests.jl`.
- Removed the corresponding include from `test/nested/integration_runtests.jl`.
- Removed only the stale component-smoke helper/call section from `test/nested/pqs_projected_q_shell_local_layer_integration_runtests.jl`; preserved the surrounding PQS source-box pair/operator, density-density, nuclear, geometry, and integration checks.

Post-deletion caller audit:
- Command: `rg -n "component_route_smoke|component_smoke_sidecars|pqs_component_route_report_adapter" src test`
- Result: no live hits in `src` or `test`.

Scoped line-count result:
- Command: `git diff --numstat -- src test bin tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl`
- Result:
  - `0` added, `1` deleted in `src/CartesianContractedParentMetrics.jl`
  - `0` added, `2168` deleted in `src/cartesian_contracted_parent_metrics/component_smoke_sidecars.jl`
  - `0` added, `1` deleted in `test/nested/integration_runtests.jl`
  - `0` added, `1642` deleted in `test/nested/pqs_component_route_report_adapter_runtests.jl`
  - `0` added, `291` deleted in `test/nested/pqs_projected_q_shell_local_layer_integration_runtests.jl`
  - total: `0 added`, `4103 deleted`, net `-4103`.

Validation:
- `rg -n "component_route_smoke|component_smoke_sidecars|pqs_component_route_report_adapter" src test`
  - passed; no live hits.
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - passed; printed `load ok` after precompile.
- `git diff --check`
  - passed.
- `git diff --cached --check`
  - passed.

Git status before writing this response:
```text
## main...origin/main
 M src/CartesianContractedParentMetrics.jl
 D src/cartesian_contracted_parent_metrics/component_smoke_sidecars.jl
 M test/nested/integration_runtests.jl
 D test/nested/pqs_component_route_report_adapter_runtests.jl
 M test/nested/pqs_projected_q_shell_local_layer_integration_runtests.jl
```

Deletion/shrinkage report:
- deleted: `src/cartesian_contracted_parent_metrics/component_smoke_sidecars.jl`; `test/nested/pqs_component_route_report_adapter_runtests.jl`.
- simplified: removed the integration-runner include and the one component-smoke helper/call block from `test/nested/pqs_projected_q_shell_local_layer_integration_runtests.jl`.
- quarantined: no replacement path was added; future H2 supplement work remains blocked on the private physical supplement preflight/payload boundary from pass 222, not this stale route-shadow/CR2 sidecar surface.
- not deleted because: live mathematical PQS source-box pair/operator, density-density, nuclear attraction, q-shell geometry, and accepted endpoint tests were intentionally left in place.
- exact remaining caller/blocker: no remaining live caller for `component_route_smoke`, `component_smoke_sidecars`, or `pqs_component_route_report_adapter` in `src` or `test`.

-- repo-doer@macmini
