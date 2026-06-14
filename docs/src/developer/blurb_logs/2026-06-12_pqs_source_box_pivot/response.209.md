Pass 209 response - physical H2 source-plan construction input audit

No source, test, docs, input, or generated files were edited for the audit.

Git status before response files:

```text
## main...origin/main
```

Files/functions inspected:
- `src/cartesian_shellification_plan.jl`
  - `_cartesian_shellification_plan_atom_growth_complete_rectangular_low_order`
  - `_cartesian_materialize_atom_growth_complete_rectangular_sequence_low_order`
  - `_cartesian_materialize_atom_growth_complete_rectangular_shellification_low_order`
  - `_cartesian_materialize_atom_local_child_shellification_low_order`
  - `_cartesian_materialize_shared_complete_shell_region`
  - `_cartesian_materialize_direct_box_region`
  - `_cartesian_materialize_source_backed_diatomic_shellification`
- `src/pqs_source_box_low_order_materialization.jl`
  - `_pqs_source_box_route_driver_diatomic_atom_growth_materializer_probe`
  - `_pqs_source_box_route_driver_diatomic_atom_growth_basis_adapter`
- `src/bond_aligned_diatomic_geometry.jl`
  - `_bond_aligned_nested_shell_provenance`
  - `_bond_aligned_source_region_points`
- `src/pqs_source_box_route_driver_skeletons.jl`
  - `_pqs_source_box_route_driver_physical_gausslet_core_shell_skeleton`
- `src/pqs_source_box_diatomic_complete_core_shell.jl`
  - `_pqs_source_box_route_driver_diatomic_physical_gausslet_source_plan_payload`
  - existing diagnostic `_PQSDiatomicCompleteCoreShellSourcePlan`
- Context inspected:
  - `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.200.md`
  - `test/nested/bond_aligned_diatomic_high_order_recipe_opt_in_source_construction_integration_runtests.jl`
  - `test/nested/bond_aligned_diatomic_atom_growth_construction_plan_runtests.jl`

Atom-contact-core rows/coefs:
- The exact old 463-column source-backed H2 fixed-source inventory was already
  observed in pass 200 from `bond_aligned_diatomic_nested_fixed_source(...;
  nside = 5)`: one child/core sequence with support count `275` and retained
  count `251`.
- In the source-backed/old object, the rows and coefs are represented by:
  - `source.child_sequences[*].support_indices`
  - `source.child_sequences[*].coefficient_matrix`
  - `source.sequence.core_column_range`
  - `source.child_column_ranges`
  - plus direct midpoint/contact support if present through
    `source.geometry.shared_midpoint_box`
- For the old gausslet-only 463 target, pass 200 says there is no separate
  midpoint retained unit; the 25-row contact plane belongs inside
  `atom_contact_core`.
- The newer atom-growth complete-rectangular materializer can also construct
  analogous pieces:
  - `_cartesian_materialize_atom_growth_complete_rectangular_sequence_low_order`
    builds `left_sequence`, optional `contact_cap_data`, and `right_sequence`;
  - it returns `core_support_blocks`, `core_support_indices`,
    `child_sequences`, `contact_cap_data`, `sequence`, `child_column_ranges`,
    and `contact_cap_column_range`;
  - coefficient blocks exist inside `left_sequence.coefficient_matrix`,
    `right_sequence.coefficient_matrix`, `contact_cap_data.coefficient_matrix`,
    and the merged `sequence.coefficient_matrix`.
- However, this is not currently carried by the H2 physical target source-plan
  seam or skeleton. The physical skeleton still carries only reviewed counts:
  `source_box = nothing`, `source_dimensions = nothing`, support count `275`,
  retained range `1:251`, retained count `251`.

Shared-shell rows/coefs:
- The old/source-backed object carries shared shells as:
  - `source.shared_shell_layers`
  - each layer's `support_indices`, optional `support_states`, and
    `coefficient_matrix`
  - `source.sequence.layer_column_ranges`
  - provenance through `_bond_aligned_nested_shell_provenance`
  - exported points through `_bond_aligned_source_region_points`
- Pass 200 confirmed the old 463 target shared-shell support counts are
  `(578, 362)` and retained counts are `(98, 114)`, with column ranges
  `252:349` and `350:463`.
- The atom-growth complete-rectangular materializer constructs shared layers
  through `_cartesian_materialize_shared_complete_shell_region`; each layer is
  checked against the region support/retained counts and returned in
  `materialization.assembly.shared_shell_layers`, with merged
  `shared_shell_column_ranges = Tuple(sequence.layer_column_ranges)`.
- The route-configured/high-order opt-in tests show related atom-growth/shared
  policies can produce different retained dimensions such as `523` or `589`.
  Those are not the 463 physical gausslet target and should not be silently
  treated as the physical source plan.

Can the existing materialization path be source authority?
- Not yet.
- The old/source-backed fixed-source path is the count/order oracle for 463,
  not the new PQS route authority.
- The atom-growth complete-rectangular materializer is useful construction
  machinery, but it is explicitly marked as low-order/private:
  - `private_development_only = true`
  - `active_source_authority = false`
  - `route_behavior_changed = false`
- The route-configured probe in
  `_pqs_source_box_route_driver_diatomic_atom_growth_materializer_probe` also
  computes `active_source_authority` from the scaffold/materialization flags,
  and the current H2 physical target run does not request or consume this probe
  as the physical source-plan payload.
- So the coherent status is: rows/coefs exist in old/source-backed and optional
  low-order materializer objects, but no route-owned physical H2 PQS source-plan
  producer currently packages the exact `(275, 578, 362)` support split and
  `(251, 98, 114)` retained split as authority.

Recommended pass 210 seam:
- Add a narrow private physical-source-plan producer that consumes an explicit
  source-backed 463 oracle object or an atom-growth materialization only behind
  a checked adapter.
- The first safe implementation should not claim full success unless it can
  verify:
  - support order exactly
    `(:atom_contact_core, :shared_shell_1, :shared_shell_2)`;
  - support counts exactly `(275, 578, 362)`;
  - retained counts exactly `(251, 98, 114)`;
  - retained order/ranges exactly core then shared layers;
  - no residual/GTO supplement;
  - no H2 221 diagnostic object kind or support order;
  - the materializer path is explicitly labeled adapter/oracle until reviewed.
- Smallest useful implementation pass:
  - build an internal adapter summary that locates candidate rows/coefs from
    either `source.child_sequences`/`source.shared_shell_layers` or
    `materialization.assembly`;
  - validate counts/order and expose only compact source-plan readiness fields
    to the existing physical source-plan payload;
  - keep `source_plan` blocked if any count/order/convention mismatch remains.

Blocker labels to use:
- If the next pass does not materialize the old/source-backed or atom-growth
  candidate at all:
  - `:missing_physical_gausslet_source_plan_materializer`
- If atom-contact support rows cannot be located:
  - `:missing_atom_contact_core_support_rows`
- If atom-contact rows exist but coefficients/range cannot be located or
  validated:
  - `:missing_atom_contact_core_coefficients`
- If shared-shell support rows exist but layer coefficients/ranges do not:
  - `:missing_shared_shell_layer_coefficients`
- If the only candidate is the low-order materializer and it is not yet adopted
  as route authority:
  - `:low_order_materializer_not_source_authority`
- If candidate counts do not match the 463 physical target:
  - `:physical_gausslet_source_plan_count_mismatch`

Deletion/shrinkage opportunity noticed:
- The route-configured atom-growth materializer/report path still exposes a
  large number of flat diagnostic aliases in
  `src/pqs_source_box_low_order_materialization.jl`. Once a compact physical
  source-plan candidate/adaptor object exists, those aliases are candidates for
  shrinkage into one summary object. No deletion was performed in this no-edit
  audit.

Validation:
- Read-only `rg`, `sed`, and `git status` commands only.
- No Julia command was run.

-- repo-doer@macmini
