Pass 145 response - audit diatomic source-plan adapter contract

Scope:
- No source edits.
- No tests or Julia commands.
- Read-only audit plus `git status --short --branch`.

Read-note:
- Requested file `docs/src/developer/fixture_role_policy_pqs_h1_h1j_rhf.md` is not present.
- The framework points to `docs/src/developer/pqs_source_box_fixture_policy.md`, so I read that available fixture-policy file instead.

Existing source-plan consumer contract:
- `pqs_multilayer_complete_core_shell_final_basis(plan)`:
  - Requires `object_kind === :pqs_multilayer_shell_source_plan` and `status === :available_pqs_multilayer_shell_source_plan` (`src/pqs_multilayer_complete_core_shell_h1.jl:18-29`).
  - Consumes `plan.metrics`, `plan.core_support_states`, and `plan.shell_support_states` to build core/core, core/shell, and shell/shell overlap blocks (`src/pqs_multilayer_complete_core_shell_h1.jl:31-52`).
  - Passes `plan.core_support_indices`, `plan.shell_support_indices`, and `plan.shell_final_coefficients` to the lower final-basis helper (`src/pqs_multilayer_complete_core_shell_h1.jl:69-80`).
  - Reports `plan.object_kind`, `plan.status`, `plan.layer_count`, and support-index lengths (`src/pqs_multilayer_complete_core_shell_h1.jl:81-99`).
- `pqs_multilayer_complete_core_shell_h1_payload(plan; final_basis, ...)`:
  - Requires the same old source-plan object kind/status plus an available final basis (`src/pqs_multilayer_complete_core_shell_h1.jl:121-127`).
  - Calls `pqs_multilayer_support_kinetic_matrix(plan)` and `pqs_multilayer_support_electron_nuclear_by_center_matrices(plan; ...)` (`src/pqs_multilayer_complete_core_shell_h1.jl:129-137`).
  - Transfers support operators through `pqs_complete_core_shell_final_one_body_matrix(final_basis, ...)` (`src/pqs_multilayer_complete_core_shell_h1.jl:138-171`).
- H1/J density path:
  - Requires old object kind/status, available final basis, and materialized H1 payload (`src/pqs_multilayer_complete_core_shell_h1.jl:242-254`).
  - Calls support-density helpers for weights and raw pair numerator matrix (`src/pqs_multilayer_complete_core_shell_h1.jl:256-275`).
  - Driver density input helper separately requires an available `:pqs_multilayer_shell_source_plan` with `bundles` (`src/pqs_source_box_route_driver_helpers.jl:10841-10924`).
- Support-one-body helpers:
  - Require `object_kind === :pqs_multilayer_shell_source_plan` and available status (`src/pqs_multilayer_support_one_body.jl:196-201`, `src/pqs_multilayer_support_one_body.jl:247-251`).
  - Use `vcat(plan.core_support_states, plan.shell_support_states)` and `plan.metrics` (`src/pqs_multilayer_support_one_body.jl:203-205`, `src/pqs_multilayer_support_one_body.jl:257-278`).
  - Label support ordering as `:core_support_states_then_shell_support_states` (`src/pqs_multilayer_support_one_body.jl:287-288`, `src/pqs_multilayer_support_one_body.jl:331-332`).
- Support-density helpers:
  - Require old object kind/status and use `vcat(plan.core_support_states, plan.shell_support_states)` (`src/pqs_multilayer_support_density.jl:3-10`).
  - Weight and pair-factor helpers assume that same combined support ordering (`src/pqs_multilayer_support_density.jl:17-32`, `src/pqs_multilayer_support_density.jl:66-92`).
- Lower complete-core/shell final-basis helper:
  - Does not require the old plan object, but it does require actual core/shell support indices, overlap blocks, and actual `shell_final_coefficients` (`src/cartesian_final_basis_realization/pqs_complete_core_shell_final_basis.jl:15-24`).
  - Requires disjoint core/shell support and coefficient row count matching shell support count (`src/cartesian_final_basis_realization/pqs_complete_core_shell_final_basis.jl:26-59`).
  - Builds `support_row_order = :core_then_shell`, pre-final coefficients, combined Lowdin cleanup, and final coefficients (`src/cartesian_final_basis_realization/pqs_complete_core_shell_final_basis.jl:116-146`, `src/cartesian_final_basis_realization/pqs_complete_core_shell_final_basis.jl:170-223`).

Old source-plan semantic contract:
- `_pqs_multilayer_realize_shell_source_plan(...)` is a one-center nested shell plan:
  - Inputs are `bundles`, rectangular `core_box`, rectangular `outer_box`, and shell layer specs (`src/pqs_multilayer_shell_source_plan.jl:41-50`).
  - Builds core support from `core_box` and each shell layer from inner/current boxes (`src/pqs_multilayer_shell_source_plan.jl:52-106`).
  - Collapses shell records and checks duplicates, core-shell overlap, and rectangular `outer_box` coverage (`src/pqs_multilayer_shell_source_plan.jl:109-130`).
  - Returns `object_kind = :pqs_multilayer_shell_source_plan` with `core_box`, `outer_box`, `layer_count`, `shell_records`, support arrays, and `shell_final_coefficients` (`src/pqs_multilayer_shell_source_plan.jl:154-180`).
- The focused H1 fixture validates exactly that one-center shape:
  - `inner_box = (2:6, 2:6, 2:6)`, `current_box = (1:7, 1:7, 1:7)`, source plan from shellification-backed region plan (`test/nested/pqs_direct_retained_final_h1_runtests.jl:33-80`).
  - It asserts region/core/outer box agreement, shell-layer count, outer support coverage, support counts, and final dimension 223 (`test/nested/pqs_direct_retained_final_h1_runtests.jl:114-139`, `test/nested/pqs_direct_retained_final_h1_runtests.jl:188-194`).

What the new diatomic source-realization payload already satisfies:
- It records the diatomic route family/system/bond-axis status and upstream payload statuses (`src/pqs_source_box_route_driver_helpers.jl:12338-12367`).
- It records the intended unit mapping:
  - core/body unit `:product`;
  - shell/source units `(:pqs_left, :pqs_right)`;
  - retained order `(:pqs_left, :pqs_right, :product)`;
  - support order `(:product, :pqs_left, :pqs_right)`;
  - explicit retained/support-order permutation requirement (`src/pqs_source_box_route_driver_helpers.jl:12369-12385`).
- On the probe-enabled path it records route retained ranges, precleanup support-order ranges, counts, shell retained count, precleanup retained dimension, and shell coefficient shape (`src/pqs_source_box_route_driver_helpers.jl:12421-12444`).
- It records the nonclaim and convention flags:
  - `object_kind_claim = :not_pqs_multilayer_shell_source_plan`;
  - `bundles_role = :parent_axis_bundle`;
  - `source_plan_materialized = false`;
  - `returns_pqs_multilayer_shell_source_plan = false`;
  - no final-basis/H1/H1-J/Ham/public/export/artifact materialization (`src/pqs_source_box_route_driver_helpers.jl:12469-12530`).
- The focused Be2 fingerprint test asserts the available probe path and the source-plan blocker `:missing_pqs_multilayer_shell_source_plan_adapter_contract` (`test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl:441-519`).

Required fields still missing as real structured data:
- Real `bundles` and `metrics` fields on a source-plan object.
  - The payload records `bundles_role`, but does not carry `bundles` or metrics.
- Real `core_support_indices` and `core_support_states`.
  - Counts exist; arrays from `producer.product_unit` are not surfaced on a source-plan object.
- Real `shell_support_indices` and `shell_support_states`.
  - Counts exist; left/right raw-box support arrays are not materialized into a combined source-plan order.
- Real `shell_final_coefficients`.
  - Shape and block-structure metadata exist; the actual `(250, 196)` block-diagonal coefficient matrix does not.
- Explicit final-basis/pre-final transform semantics.
  - The lower final-basis helper will build combined Lowdin cleanup from the pre-final overlap, but the diatomic object needs to define the pre-final column order and route-retained-to-pre-final mapping as a structured object, not only ranges.
- Product/core/body placement as real source-plan data.
  - The payload says `:product` is core/body; it does not carry a source-plan object with product support arrays and convention labels consumed by one-body/density helpers.
- Object kind and convention labels.
  - The old helpers hard-gate on `:pqs_multilayer_shell_source_plan`; the diatomic payload explicitly refuses that claim.
- Old one-center fields.
  - `core_box`, `outer_box`, `layer_count`, and `shell_records` are not meaningfully represented for product slab + separated left/right PQS sectors.

Decision:
- Returning `object_kind = :pqs_multilayer_shell_source_plan` for Be2/PQS is not semantically correct.
- A narrow tuple could be made to satisfy current guard clauses by filling support arrays and coefficient matrices, but it would falsely overload old meanings:
  - `core_box` would not be the product slab/body sector in the old sense;
  - `outer_box` would be a rectangle containing geometry not owned by the left/product/right source sectors;
  - `layer_count` and `shell_records` would no longer mean nested shellification layers;
  - support coverage checks from the old source-plan constructor would not describe the diatomic route.
- This would violate the framework's source-box-first boundary and the warning that shell/support-row realization is compatibility/oracle, not PQS route authority (`docs/src/developer/pqs_source_box_operator_framework.md:18-36`, `docs/src/developer/pqs_source_box_operator_framework.md:185-189`).

Smallest new private object/consumer seam:
- Add a private Be2/PQS source-plan object instead of adapting to `:pqs_multilayer_shell_source_plan`.
- Suggested object/helper:
  - `_PQSDiatomicCompleteCoreShellSourcePlan`
  - `_pqs_source_box_route_driver_diatomic_complete_core_shell_source_plan`
- It should materialize real structured source-plan data only:
  - `object_kind = :pqs_diatomic_complete_core_shell_source_plan`
  - `status`
  - `blocker`
  - `bundles`
  - `metrics`
  - `core_unit_key = :product`
  - `shell_unit_keys = (:pqs_left, :pqs_right)`
  - `core_support_indices`
  - `core_support_states`
  - `shell_support_indices`
  - `shell_support_states`
  - `shell_final_coefficients`
  - `support_order = (:product, :pqs_left, :pqs_right)`
  - `route_retained_order = (:pqs_left, :pqs_right, :product)`
  - structured retained/pre-final permutation or range map
  - source unit summaries/references
  - convention labels and nonpromotion flags
  - `metadata`
- It should not build final basis, H1, H1-J, Ham, RHF/SCF, exports, or artifacts.

Recommended next implementation pass:
- Implement the private `_PQSDiatomicCompleteCoreShellSourcePlan` materializer from the source-realization/raw-box payload, still with no final-basis/H1/H1-J/Ham materialization.
- Update the focused Be2 fingerprint test to assert:
  - available probe-enabled diatomic source plan;
  - real support array lengths and disjointness;
  - real shell coefficient matrix shape `(250, 196)`;
  - object kind is not `:pqs_multilayer_shell_source_plan`;
  - source-plan payload now blocks on `:missing_diatomic_complete_core_shell_final_basis_consumer` or equivalent, not on the old adapter contract.
- Leave existing one-center complete-core/shell H1 fixture untouched.
- A later pass can decide whether to:
  - add a generic complete-core/shell support-plan consumer shared by one-center and diatomic objects, or
  - add a diatomic-specific final-basis/H1 adapter that calls lower-level helpers without claiming old source-plan object kind.

Validation:
- Read-only inspection only.
- `git status --short --branch`:
```text
## main...origin/main
```

Deletion/shrinkage forecast:
- deleted: none in this no-edit audit.
- simplified: next pass can replace `:missing_pqs_multilayer_shell_source_plan_adapter_contract` with a sharper private diatomic source-plan object and a later consumer blocker.
- quarantined: old `:pqs_multilayer_shell_source_plan` adapter claim for Be2/PQS; final-basis/H1/H1-J/Ham/RHF/public/export/artifact behavior.
- not deleted because: one-center `pqs_multilayer_shell_source_plan` remains the validated compact H1 fixture contract and should not be disturbed by Be2/PQS route-shape work.
- exact remaining caller/blocker: `pqs_multilayer_complete_core_shell_final_basis`, `pqs_multilayer_complete_core_shell_h1_payload`, `pqs_multilayer_support_*`, and complete-core/shell density/H1-J helpers still hard-require `object_kind = :pqs_multilayer_shell_source_plan`; Be2/PQS should first get `_PQSDiatomicCompleteCoreShellSourcePlan`, then a new/generic consumer seam.

-- repo-doer@macmini
