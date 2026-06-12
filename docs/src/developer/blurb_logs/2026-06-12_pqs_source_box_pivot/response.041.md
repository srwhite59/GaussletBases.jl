Pass 041 complete.

Summary:
- Created and ran an ignored audit/probe script:
  - `tmp/work/pqs_complete_one_center_h1_probe.jl`
- Wrote the ignored artifact:
  - `tmp/work/pqs_complete_one_center_h1_probe_summary.txt`
- No tracked source, test, or docs files were changed.

Complete-basis strategy:
- Preferred target selected:
  - `5 x 5 x 5` inner core plus one surrounding PQS boundary shell.
- Expected retained/final dimension for that target:
  - core count: `125`
  - boundary shell count: `98`
  - combined complete one-center count: `223`
- This is the correct next physical target; the pass 040 boundary-shell-only `98`-function H1 fixture should not be promoted as a physical acceptance target.

Audit result:
- Existing production surfaces cannot yet express the preferred complete one-center PQS/final basis cleanly.
- Exact blocker:
  - `:missing_combined_direct_core_boundary_shell_final_basis_realization`
- Missing surfaces recorded by the probe:
  - `:combined_direct_core_boundary_shell_final_basis_realization`
  - `:route_owned_complete_one_center_shell_projection_lowdin`
  - `:direct_retained_operator_blocks_for_combined_core_shell_basis`
  - `:all_source_mode_retained_rule_if_single_box_fallback_is_desired`

All-source fallback audit:
- Diagnostic fallback considered:
  - one `5 x 5 x 5` source box with all `125` source modes retained.
- Result:
  - `all_source_fallback_status = :blocked_as_expected`
  - `all_source_fallback_blocker = :missing_all_source_mode_retained_rule_contract`
  - exact rejection: `ArgumentError: PQS final basis requires the boundary source-mode retained rule`
- Reason:
  - `CartesianFinalBasisRealization.pqs_source_shell_realization_final_basis(...)` currently validates specifically for `:boundary_comx_product_mode_selection`.
  - Recent direct retained operators also require the boundary source-mode rule, so a non-boundary all-source fallback is not a valid active path without a deliberate new rule/operator contract.

Fixture/extent facts:
- source-mode dims/count:
  - `(5, 5, 5) / 125`
- current boundary shell retained count:
  - `98`
- diagnostic all-source retained count:
  - `125`
- physical extent is known:
  - axis centers: `(-7.999999999999998, -1.9955875976247728, 0.0, 1.9955875976247728, 7.999999999999998)`
  - axis min/max: about `-8.0 / 8.0`

H1 result:
- H1 was not materialized for the complete one-center target.
- exact H reference retained in artifact: `-0.5`
- exact-reference error: unavailable because the complete final basis is blocked.
- shell-support oracle delta: unavailable for the complete target.

Route/nonclaim flags:
- `_pqs_current_route_safe_term_matrices(...)` active path used: `false`
- generalized-overlap final solve used: `false`
- IDA materialized: `false`
- density-density materialized: `false`
- RHF materialized: `false`
- driver route materialized: `false`
- GTO/export/artifact work: `false`

Validation:
- `julia --project=. tmp/work/pqs_complete_one_center_h1_probe.jl` passed and produced the blocked audit artifact.
- `julia --project=. -e 'using GaussletBases; println("load ok")'` passed.
- `git diff --check` passed.

Deletion/shrinkage report:
- Nothing was deleted or shrunk.
- This was probe/audit-only because the preferred complete one-center route is blocked by missing final-basis representation, not by a small numerical bug.
- Exact next implementation target:
  - add a narrow combined direct-core plus PQS boundary-shell final-basis realization object/path for the one-center `125 + 98 = 223` fixture, with route-owned shell projection/Lowdin data and direct retained overlap/kinetic/by-center nuclear block placement for the combined basis.
- Do not proceed to final IDA weights/RHF until that complete one-center H1 path exists or a reviewed alternative complete-basis contract is chosen.

-- repo-doer@macmini
