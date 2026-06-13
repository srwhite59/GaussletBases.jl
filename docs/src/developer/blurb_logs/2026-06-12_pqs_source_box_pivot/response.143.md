Pass 143 response - raw-box to source-plan realization mapping audit

Role/task boundary:
- repo-doer@macmini.
- No-edit audit only; no source/test edits, no commits, no broad tests.
- Governing framework kept source-box-first PQS as route framing and left shell/support-row contraction as oracle/debug.

Files/helpers inspected:
- `src/pqs_source_box_route_driver_helpers.jl`
  - `_PQSDiatomicCompleteCoreShellSupportWindowPayload`
  - `_PQSDiatomicRawBoxRoutePayload`
  - `_PQSDiatomicCompleteCoreShellSourcePlanPayload`
  - `_pqs_source_box_route_driver_diatomic_complete_core_shell_support_window_payload`
  - `_pqs_source_box_route_driver_diatomic_raw_box_route_payload`
  - `_pqs_source_box_route_driver_diatomic_complete_core_shell_source_plan_payload`
  - `cartesian_assembly`
- `src/CartesianContractedParentMetrics.jl`
  - `_pqs_pqs_product_raw_box_route_producer`
  - `_pqs_product_doside_identity_slab_unit`
  - `_pqs_parent_coefficient_matrix_from_raw_plan`
  - `_product_doside_parent_coefficient_matrix`
  - `_pqs_pqs_product_route_parent_coefficient_matrix`
- `src/pqs_multilayer_shell_source_plan.jl`
  - `_pqs_multilayer_realize_shell_source_plan`
  - `pqs_multilayer_shell_source_plan`
- `src/pqs_multilayer_complete_core_shell_h1.jl`
  - `pqs_multilayer_complete_core_shell_final_basis`
  - `pqs_multilayer_complete_core_shell_h1_payload`
- `src/pqs_multilayer_support_one_body.jl`
  - support-space kinetic and electron-nuclear helpers
- `src/pqs_multilayer_support_density.jl`
  - support-state and support-weight helpers
- `src/cartesian_final_basis_realization/pqs_complete_core_shell_final_basis.jl`
  - `pqs_complete_core_shell_final_basis`
- `test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl`

Proposed raw-box to source-plan mapping:
- `producer.product_unit.support_states` and `producer.product_unit.support_indices` can honestly serve as the direct core/body support states and indices for a Be2-specific complete core/shell realization if the route contract explicitly defines the product slab as the core/body sector.
- This is not the same semantic `core_box` used by the existing one-center/multilayer shell source plan. It is a product/doside slab sector in a diatomic PQS route.
- Left/right raw PQS plan views provide enough data for shell/source sectors:
  - each raw PQS plan has source intervals, source-mode dimensions, axis-local coefficients, and boundary selector modes/count;
  - `_pqs_parent_coefficient_matrix_from_raw_plan` already expands selected boundary modes to parent rows;
  - a source-realization object can derive support-local row sets by the raw source boxes and retain the corresponding selected-mode coefficient blocks.
- The existing retained route order is `(:pqs_left, :pqs_right, :product)`.
- The complete core/shell support order expected by current final-basis/support helpers is core first, then shell: `(:product, :pqs_left, :pqs_right)`.
- Therefore the retained-order/support-order permutation is real and should be explicit route-owned metadata, not implicit in final-basis or H1 helpers.

Expected shell coefficient shape:
- Product/core support count: 25.
- Left raw PQS support count: 125; retained boundary count: 98.
- Right raw PQS support count: 125; retained boundary count: 98.
- `shell_support_states = vcat(left_support_states, right_support_states)` should have 250 rows.
- `shell_final_coefficients` should be block diagonal over left/right PQS sectors with shape `(250, 196)`:
  - left block `(125, 98)`;
  - right block `(125, 98)`;
  - zero cross blocks.
- The complete pre-final support layout would be product/core rows first, then left/right shell rows. Pre-cleanup retained columns would be product/core 25 plus shell 196 = 221.

Where the permutation should live:
- Put it in the next private diatomic source-realization payload.
- Suggested fields:
  - `retained_order = (:pqs_left, :pqs_right, :product)`
  - `support_order = (:product, :pqs_left, :pqs_right)`
  - `core_unit_key = :product`
  - `shell_unit_keys = (:pqs_left, :pqs_right)`
  - `route_retained_ranges = producer.descriptor.expected_ranges`
  - `source_plan_precleanup_ranges = (product = 1:25, pqs_left = 26:123, pqs_right = 124:221)`
  - `retained_to_support_order_permutation_required = true`
  - a compact permutation/range summary, not a flat cloud on `cartesian_assembly`.

`bundles` mapping:
- For current density/H1/final-basis consumers, `bundles` means the parent axis bundle object and `metrics` means parent-axis PGDG factors derived from that bundle.
- It should not be the raw PQS plan objects and does not require a new diatomic bundle object for the next pass.
- The raw PQS plans are source/retained sector objects; the support-space operators still need parent-axis metrics.

Is returning `:pqs_multilayer_shell_source_plan` honest yet?
- No.
- The existing object is a one-center/multilayer rectangular shell plan:
  - `core_box` and `outer_box` describe nested boxes;
  - `layer_count` and `shell_records` describe shell layers;
  - the realization checks duplicate shell support, core-shell overlap, and coverage of a rectangular `outer_box`.
- The Be2/PQS route is a diatomic source-box route with product slab core/body sector plus separated left/right PQS source sectors. Its natural object vocabulary is product/left/right units, not nested core/outer shell layers.
- Even though the support arrays and coefficient block can likely be built, returning `object_kind = :pqs_multilayer_shell_source_plan` now would overload `core_box`, `outer_box`, `layer_count`, and `shell_records`.

Recommended next coding pass:
- Add a private route-owned source-realization payload, not a source-plan materializer yet.
- Suggested name:
  - `_PQSDiatomicCompleteCoreShellSourceRealizationPayload`
  - helper `_pqs_source_box_route_driver_diatomic_complete_core_shell_source_realization_payload`
- Required inputs:
  - `parent`
  - `route_skeleton`
  - `recipe`
  - `_PQSDiatomicCompleteCoreShellSupportWindowPayload`
  - `_PQSDiatomicRawBoxRoutePayload`
- Suggested fields:
  - `status`
  - `blocker`
  - `route_family`
  - `system_classification`
  - `bond_axis`
  - `support_window_payload_status`
  - `raw_box_route_payload_status`
  - `core_unit_key`
  - `shell_unit_keys`
  - `retained_order`
  - `support_order`
  - `route_retained_ranges`
  - `source_plan_precleanup_ranges`
  - `core_support_count`
  - `shell_support_counts`
  - `shell_support_count`
  - `shell_retained_counts`
  - `shell_retained_count`
  - `precleanup_retained_dimension`
  - `shell_final_coefficients_shape`
  - `shell_coefficient_block_structure`
  - `parent_axis_bundle_object_available`
  - `bundles_role = :parent_axis_bundle`
  - `object_kind_claim = :not_pqs_multilayer_shell_source_plan`
  - `available_objects`
  - `missing_objects`
  - `summary`
  - `metadata`
- Status behavior:
  - default/no parent bundle path stays blocked through raw payload with `:missing_parent_axis_bundle_object`;
  - probe-enabled path can become `:available_diatomic_complete_core_shell_source_realization_payload`;
  - source-plan payload should still block with `:missing_pqs_multilayer_shell_source_plan_adapter_contract` or equivalent until the adapter/new consumer contract is approved.
- Focused validation for that pass:
  - `julia --project=. test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl`
  - `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - `git diff --check`

What it would replace or shrink:
- It would stop `_pqs_source_box_route_driver_diatomic_complete_core_shell_source_plan_payload` from carrying only a broad missing `:diatomic_complete_core_shell_source_realization_contract` without a compact intermediate authority object.
- It would concentrate source-realization ordering, support counts, coefficient block shape, and object-kind nonclaim in one route-owned payload.
- It would let later source-plan/materializer work consume one payload instead of manually coordinating support-window payload, raw-box producer summaries, descriptor ranges, and ad hoc ordering facts.

Blockers:
- Need an explicit adapter/new-object contract before claiming `:pqs_multilayer_shell_source_plan`.
- Need a reviewed decision on whether the existing one-center complete-core/shell consumers should accept a new diatomic source-realization object or whether a narrow adapter should emit the old consumer shape with documented semantic translation.
- Need explicit retained-order to support-order permutation/range mapping before final-basis/H1 materialization.

Validation:
- Read-only inspection only, per blurb.
- Commands used:
  - `sed` over startup/framework docs and audited source/test files.
  - `rg` over audited helper names and source-plan consumer fields.
  - `git status --short --branch`
- No Julia commands or tests were run.

Git status:
```text
## main...origin/main
```

Deletion/shrinkage forecast:
- deleted: none in this no-edit audit; next implementation should not delete behavior.
- simplified: `_PQSDiatomicCompleteCoreShellSourcePlanPayload` missing-object/report logic can shrink by consuming one source-realization payload summary.
- quarantined: `:pqs_multilayer_shell_source_plan` return claim; final-basis/H1/H1-J/Ham materialization; support-row contraction authority; retained diagnostic weights as IDA weights.
- not deleted because: existing complete-core/shell H1/final-basis helpers are still the validated oracle-backed consumer shape and the Be2 route does not yet have an approved adapter/new consumer.
- exact remaining caller/blocker: `pqs_multilayer_complete_core_shell_final_basis`, `pqs_multilayer_complete_core_shell_h1_payload`, support one-body, and support density helpers require an available `object_kind = :pqs_multilayer_shell_source_plan`; Be2/PQS needs a private diatomic source-realization payload plus adapter/new-object contract before those callers can consume it honestly.

-- repo-doer@macmini
