Pass 097 response

Implemented the narrow private density-input adapter for the complete core/shell
H1/J diagnostic path.

Files changed:

- `src/pqs_source_box_route_driver_helpers.jl`
- `.agent_handoffs/response.097.md`
- `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.097.md`

Helper added:

- `_pqs_source_box_route_driver_complete_core_shell_density_inputs(source_plan, coulomb_expansion; metadata = (;))`

Returned fields:

- `status`
- `blocker`
- `axis_weights`
- `raw_pair_factor_terms`
- `missing_inputs`
- `summary`
- `metadata`

Implementation notes:

- Requires `source_plan.object_kind === :pqs_multilayer_shell_source_plan`.
- Requires `source_plan.status === :available_pqs_multilayer_shell_source_plan`.
- Requires structured `source_plan.bundles`.
- Requires `coulomb_expansion.coefficients`.
- Calls `CartesianContractedParentMetrics._pqs_source_box_ida_factor_provenance(bundles; expected_term_count = length(coulomb_expansion.coefficients))`.
- Returns `provenance.axis_weights` and `provenance.raw_axis_pair_factor_terms`.
- On missing/invalid structured inputs, returns blocked density inputs and lets the H1/J diagnostic payload remain blocked.
- Does not fall back to density-normalized pair terms, retained diagnostic weights, report scalar aliases, fixed-block data, or shell/support-row oracle authority.

Wiring:

- `_pqs_source_box_route_driver_complete_core_shell_diagnostic_route_payload(...)` now builds `density_inputs` after the source-plan and H1 payload helpers.
- It passes:
  - `axis_weights = density_inputs.axis_weights`
  - `raw_pair_factor_terms = density_inputs.raw_pair_factor_terms`
  into `_pqs_source_box_route_driver_complete_core_shell_h1_j_diagnostic_payload(...)`.
- Report-facing fields and aliases are unchanged.

H1/J diagnostic status:

- Before this pass, the structured source-plan/H1 fixture remained blocked for H1/J at `:axis_weights` and `:raw_pair_factor_terms`.
- After this pass, the focused direct structured source-plan probe reports:
  - `density_inputs_status=available_complete_core_shell_density_inputs`
  - `h1_j_status=materialized_pqs_multilayer_complete_core_shell_h1_j_payload`
  - `missing_inputs=()`
  - `final_dimension=223`
  - `h1_energy=-0.48047934800387126`
  - `self_coulomb=0.6397851751855723`
  - `density_gauge=pre_final_localized_positive_weight`

Focused route-driver dry-run status:

- A focused `cartesian_assembly` route probe with `route_family=:pqs_source_box` and terminal shellification policy still blocks before source-plan construction:
  - `parent_axis_bundle_object_available=true`
  - `terminal_route_available=false`
  - `source_status=blocked_missing_complete_core_shell_source_plan_inputs`
  - `source_missing=(:terminal_shellification_plan, :terminal_lowering_plan)`
  - `h1_j_status=blocked_missing_complete_core_shell_h1_j_route_inputs`
- This means the route-driver dry run did not reach the new density-input helper with an available route-owned source plan. The direct structured source-plan probe validates the new helper itself.

Route/global flags:

- Focused route-driver probe stayed globally unmaterialized:
  - `driver_route_materialized=false`
  - `rhf_materialized=false`
  - `exports_materialized=false`
  - `artifacts_materialized=false`
  - `ida_data_materialized=false`
  - `density_density_materialized=false`
- The direct structured source-plan probe materialized only the private H1/J diagnostic payload; RHF, exports, and artifacts remained false there as well.

Validation commands/results:

- `julia --project=. test/nested/pqs_direct_retained_final_h1_runtests.jl`
  - passed: `46/46`, about `5.7s` after final edit
- focused direct structured source-plan H1/J diagnostic probe via `julia --project=. -e 'include("test/nested/pqs_direct_retained_final_h1_runtests.jl"); ...'`
  - passed focused H1 test during include
  - `density_inputs_status=available_complete_core_shell_density_inputs`
  - `h1_j_status=materialized_pqs_multilayer_complete_core_shell_h1_j_payload`
  - `final_dimension=223`
  - `h1_energy=-0.48047934800387126`
  - `self_coulomb=0.6397851751855723`
- focused `cartesian_assembly` route probe via `julia --project=. -e 'using GaussletBases; ...'`
  - completed
  - route remained blocked before complete core/shell source-plan construction as noted above
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - passed: `load ok`
- `git diff --check`
  - passed

Intermediate fix:

- An initial direct probe exposed that `_pqs_source_box_ida_factor_provenance` needed module qualification from `pqs_source_box_route_driver_helpers.jl`; fixed to `CartesianContractedParentMetrics._pqs_source_box_ida_factor_provenance(...)`.

Git status before writing this response:

```text
## main...origin/main
 M src/pqs_source_box_route_driver_helpers.jl
```

Deletion/shrinkage report:

- deleted: none
- simplified: H1/J diagnostic wiring now consumes one compact density-input payload instead of leaving axis weights and raw pair terms as adjacent missing scalar inputs.
- quarantined: none
- not deleted because: report-facing H1/J aliases and existing blocked-path summaries remain compatibility/report surfaces.
- exact remaining caller/blocker: current focused `cartesian_assembly` route-driver probe does not select terminal route state for `route_family=:pqs_source_box`; complete core/shell source-plan payload remains blocked on `(:terminal_shellification_plan, :terminal_lowering_plan)` before density inputs can be consumed in that route-driver path.

-- repo-doer@macmini
