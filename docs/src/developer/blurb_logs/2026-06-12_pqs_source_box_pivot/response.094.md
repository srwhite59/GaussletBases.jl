Pass 094 response: PQS complete core/shell region/source-plan seam

Summary:
- Added a private route-owned complete core/shell source-plan payload helper.
- The helper consumes the selected terminal shellification/lowering plans from the driver and the carried `parent.parent_axis_bundle_object`.
- `cartesian_assembly(...)` now passes available region/source plans and the standard Coulomb expansion into the existing pass-091 H1/J diagnostic helper.
- No new report fields were added. Blocked H1/J reports still expose only compact summary fields and do not retain dense source-plan payloads.

Files changed:
- `src/pqs_source_box_route_driver_helpers.jl`

Driver-owned region/source-plan status:
- Region plan: available through `pqs_multilayer_shell_region_plan(shellification_plan, lowering_plan; ...)`.
- Source plan: available through `pqs_multilayer_shell_source_plan(parent_axis_bundle_object, region_plan; ...)`.
- Source-plan helper status: `:available_complete_core_shell_source_plan_payload`.
- Source-plan blocker: `nothing`.
- Active source authority remains shellification-backed geometry; the explicit-box bridge was not used.

Focused dry-run smoke:

```
diatomic.pqs_layer.lowdin_cleanup: 0.1085 seconds
h1j_status=blocked_missing_complete_core_shell_h1_j_route_inputs
h1j_blocker=missing_complete_core_shell_h1_j_route_inputs
h1j_missing=(:pqs_multilayer_complete_core_shell_final_basis, :pqs_multilayer_complete_core_shell_h1_payload, :axis_weights, :raw_pair_factor_terms)
shellization_source=terminal_cartesian_shellification_geometry
terminal_selected=true
lowering_plan_status=available_terminal_lowering_plan
driver_route_materialized=false
```

H1/J diagnostic advancement:
- The diagnostic missing-input list no longer includes:
  - `:pqs_multilayer_shell_region_plan`
  - `:pqs_multilayer_shell_source_plan`
  - `:coulomb_expansion`
- It remains correctly blocked because this pass did not build final-basis or H1 payloads and did not feed IDA density inputs.

Remaining missing boundary:
- Final-basis stage: `pqs_multilayer_complete_core_shell_final_basis(...)`.
- H1 payload stage: `pqs_multilayer_complete_core_shell_h1_payload(...)`.
- Density/J diagnostic inputs: final axis weights and raw pair numerator terms for the complete core/shell retained density interaction.

Validation:
- Ran compact dry-run smoke/probe for the source-plan seam: passed with the missing-input list reduced to final-basis/H1/axis-weights/raw-pair-terms.
- Ran `julia --project=. -e 'using GaussletBases; println("load ok")'`: passed.
- Ran `git diff --check`: passed.
- Did not run broad CPBM/report/materialization tests.
- Did not run `test/nested/pqs_direct_retained_final_h1_runtests.jl` because final-basis/H1 construction was not changed.

Deletion/shrinkage report:
- No old code, tests, metadata, or compatibility paths became unnecessary in this pass.
- Nothing was deleted because this pass adds the first route-owned source-plan seam; the existing blocked H1/J diagnostic remains necessary to localize the next missing boundary.
- No new tests were added.
- Remaining stale/duplicate surface to retire later: after driver-owned final-basis/H1 payload construction is connected, fixture-local complete core/shell H1 construction should be kept as oracle-only or shrunk if it duplicates driver coverage.

-- repo-doer@macmini
