Pass 093 response: PQS source-box shell-stage selection

Summary:
- Added a narrow shell-stage seam for one-center `:pqs_source_box` diagnostics.
- The driver now selects the existing terminal shellification payload for one-center PQS source-box routes instead of reporting shellification as `:not_applicable`.
- The terminal shellification payload now uses `CartesianShellification.OneCenterShellification` for one-center parents, matching the complete core/shell H1 fixture policy.
- No H1/final-basis/source-plan realization was added in this pass.

Files changed:
- `src/pqs_source_box_route_driver_helpers.jl`

Focused dry-run result:

```
shellization_source=terminal_cartesian_shellification_geometry
terminal_selected=true
terminal_plan_available=true
terminal_region_count=2
lowering_plan_available=true
lowering_plan_status=available_terminal_lowering_plan
h1j_status=blocked_missing_complete_core_shell_h1_j_route_inputs
h1j_missing=(:pqs_multilayer_shell_region_plan, :pqs_multilayer_shell_source_plan, :pqs_multilayer_complete_core_shell_final_basis, :pqs_multilayer_complete_core_shell_h1_payload, :axis_weights, :raw_pair_factor_terms, :coulomb_expansion)
```

Shell-stage status:
- `:pqs_source_box` one-center shell stage now selects a complete core/shell terminal shellification plan.
- Status advanced from `low_order_shellization_source = :not_applicable` / `terminal_shellification_selected = false` to:
  - `low_order_shellization_source = :terminal_cartesian_shellification_geometry`
  - `low_order_terminal_shellification_selected = true`
  - `low_order_terminal_shellification_region_count = 2`
- The existing downstream lowering seam also now reports:
  - `low_order_terminal_shellification_lowering_plan_available = true`
  - `low_order_terminal_shellification_lowering_plan_status = :available_terminal_lowering_plan`

Remaining missing boundary for H1/J:
- Region/source stage: `pqs_multilayer_shell_region_plan(...)` and `pqs_multilayer_shell_source_plan(...)` are still not produced from driver-owned objects.
- Final-basis stage: `pqs_multilayer_complete_core_shell_final_basis(...)` is still missing.
- H1 payload stage: `pqs_multilayer_complete_core_shell_h1_payload(...)` is still missing.
- Density/Coulomb input seam: axis weights, raw pair numerator terms, and a named Coulomb expansion route input are still not fed into the pass-091 H1/J diagnostic payload helper.

Validation:
- Ran focused dry-run smoke/probe for the shell-stage seam: passed with terminal shellification and terminal lowering available.
- Ran `julia --project=. -e 'using GaussletBases; println("load ok")'`: passed.
- Ran `git diff --check`: passed.
- Did not run slow broad route-driver/report tests.
- Did not run `test/nested/pqs_direct_retained_final_h1_runtests.jl` because this pass did not touch final-basis or H1 construction.

Deletion/shrinkage report:
- No old code, tests, metadata, or compatibility paths became unnecessary yet.
- Nothing was deleted because this pass only turns on the first route-owned shell-stage selection seam; the previous blocked diagnostic path remains needed to localize the next missing stage.
- No new tests were added.
- Remaining stale/duplicate surface to retire later: once the driver owns the complete core/shell region/source/final-basis path, fixture-local reconstruction in H1 probes should be kept as oracle coverage only or shrunk.

-- repo-doer@macmini
