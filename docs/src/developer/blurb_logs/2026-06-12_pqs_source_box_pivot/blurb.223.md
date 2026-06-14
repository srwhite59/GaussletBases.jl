Pass 223 - audit and retire stale component-smoke sidecar surface

Role:
You are `repo-doer@macmini` doing one bounded cleanup/audit pass for
GaussletBases. Follow `AGENTS.md`, `JuliaStyle.md`, `BlurbStyle.md`, and the
unattended baton rules in `.agent_handoffs/RUN.md`.

Loop/approval rule:
- Unattended baton mode is active.
- Do not ask the user for permission through an escalation prompt.
- If a needed command requires approval, write `.agent_handoffs/ATTENTION.md`
  with the exact command, why it is needed, and the blocking condition, then
  stop.

Current state:
- Head before this pass should be:
  `cfe470db Retire H2 source-box diagnostic scaffold`
- H2 R=4 q=n_s=5 gausslet-only PQS/WL physical endpoint is accepted and
  comparison-ready.
- Pass 222 identified the future supplement seam as a private H2 physical
  supplement preflight/payload boundary, but no supplement implementation should
  be added yet.

Why this pass:
Before adding new supplement code, remove or quarantine old component-smoke/CR2
sidecar machinery if it is no longer live. That surface mixes source-box
diagnostic route-shadow vocabulary with final-residual MWG sidecar concepts and
can pull new supplement work toward the wrong contract.

Primary audit surfaces:

```text
src/cartesian_contracted_parent_metrics/component_smoke_sidecars.jl
test/nested/pqs_component_route_report_adapter_runtests.jl
test/nested/pqs_projected_q_shell_local_layer_integration_runtests.jl
test/nested/integration_runtests.jl
src/CartesianContractedParentMetrics.jl
```

Key helper families to audit:

```text
_pqs_pqs_product_source_box_component_route_smoke
_pqs_pqs_product_component_route_smoke_summary
_pqs_pqs_product_component_route_smoke_report_adapter
_pqs_pqs_product_component_route_smoke_cr2_sidecar_schema
_write_pqs_pqs_product_component_route_smoke_report
_write_pqs_pqs_product_component_route_smoke_cr2_sidecar_schema_report
```

Decision rule:
1. If these helpers are only called by stale route-shadow/component-smoke tests
   and historical docs, delete the source include and its test pressure:
   - remove the include from `src/CartesianContractedParentMetrics.jl`;
   - delete `src/cartesian_contracted_parent_metrics/component_smoke_sidecars.jl`;
   - delete or shrink the corresponding sections in
     `test/nested/pqs_component_route_report_adapter_runtests.jl`;
   - remove the integration-runner include if the whole test file is deleted;
   - delete matching component-smoke sections from
     `test/nested/pqs_projected_q_shell_local_layer_integration_runtests.jl`
     if they become the only remaining callers.
2. If deletion is blocked by a live non-test source caller or an accepted
   endpoint contract, do not add adapters. Report the exact caller/blocker and
   shrink only the stale assertions.
3. Preserve live mathematical tests for PQS source-box pair operators,
   density-density normalization, nuclear attraction, q-shell geometry, and the
   accepted H2 physical endpoint.

Do not:
- implement GTO/MWG supplement support;
- add the H2 supplement preflight payload yet;
- change H2 physical endpoint values;
- touch He/H2 accepted endpoint tests except for reference audits;
- add public API/export/HamV6/CR2/HFDMRG/DMRG behavior;
- move stale component-smoke code to another live path just to preserve it.

Line-count rule:
This pass must be net-negative under source/test/bin. Measure:

```sh
git diff --numstat -- src test bin tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
```

Expected line-budget source is deleting stale component-smoke source/test
surface. Do not satisfy line budget by deleting accepted endpoint tests.

Validation:
Run:

```sh
rg -n "component_route_smoke|component_smoke_sidecars|pqs_component_route_report_adapter" src test
julia --project=. -e 'using GaussletBases; println("load ok")'
git diff --check
git diff --cached --check
```

If you remove an include from `test/nested/integration_runtests.jl`, do not run
the full slow integration runner unless a focused failure points there. The load
check and `rg` caller audit are the required validation for this cleanup pass.

Response file:
Write `.agent_handoffs/response.223.md` and also copy it to:

```text
docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.223.md
```

Report:
- exact callers found before deletion/shrinkage;
- whether `component_smoke_sidecars.jl` was deleted, retained, or partially
  shrunk;
- test files deleted or shrunk;
- source/test/bin scoped added/deleted/net line count;
- validation commands and results;
- deletion/shrinkage report:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

Continue polling after writing the response. Manager will review, commit/push,
and publish the next blurb unless a real blocker appears.

-- repo-manager@macmini
