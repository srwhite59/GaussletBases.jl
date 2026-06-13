Pass 134 response - Be2 PQS complete-core/shell seam audit

Role/task boundary:
- repo-doer@macmini, pass 134 only.
- No source, test, or docs edits were made beyond this handoff response.
- No tests or Julia probes were run; read-only inspection was sufficient.

Exact helpers/files inspected:
- `src/pqs_source_box_route_driver_helpers.jl`
  - `_cartesian_one_center_parent_basis_object`
  - `_cartesian_one_center_parent_axis_bundle_object`
  - `_cartesian_parent_object_carry`
  - `_pqs_source_box_route_driver_shell_stage_terminal_shellification`
  - `_pqs_source_box_route_driver_shell_stage_low_order_shellization`
  - `_pqs_source_box_route_driver_unit_stage_low_order_summary`
  - `_pqs_source_box_route_driver_pair_stage_low_order_summary`
  - `_pqs_source_box_route_driver_complete_core_shell_source_plan_payload`
  - `_pqs_source_box_route_driver_complete_core_shell_center_records`
  - `_pqs_source_box_route_driver_complete_core_shell_axis_layers`
  - `_pqs_source_box_route_driver_complete_core_shell_h1_payload`
  - `_pqs_source_box_route_driver_complete_core_shell_density_inputs`
  - `_pqs_source_box_route_driver_complete_core_shell_h1_j_diagnostic_payload`
  - `_pqs_source_box_route_driver_complete_core_shell_ham_payload`
  - `_PQSCompleteCoreShellDiagnosticRoutePayload`
  - `_pqs_source_box_route_driver_complete_core_shell_diagnostic_route_payload`
  - `cartesian_assembly`
- `src/pqs_source_box_route_driver_skeletons.jl`
  - `_pqs_source_box_route_driver_route_skeleton`
- `src/CartesianContractedParentMetrics.jl`
  - `_pqs_pqs_product_source_box_route_skeleton`
- `src/pqs_multilayer_shell_region_plan.jl`
  - `pqs_multilayer_shell_region_plan`
- `src/pqs_multilayer_shell_source_plan.jl`
  - `pqs_multilayer_shell_source_plan`
- `src/pqs_multilayer_complete_core_shell_h1.jl`
  - `pqs_multilayer_complete_core_shell_final_basis`
  - `pqs_multilayer_complete_core_shell_h1_payload`
  - `pqs_multilayer_complete_core_shell_h1_j_payload`
- `src/cartesian_final_basis_realization/pqs_complete_core_shell_final_basis.jl`
  - `pqs_complete_core_shell_final_basis`
  - final one-body, H1 solve, IDA-weight, and pre-final density helpers
- `test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl`
- `test/nested/pqs_source_box_route_driver_report_runtests.jl`
- Governing docs:
  - `docs/src/developer/pqs_source_box_operator_framework.md`
  - `docs/src/developer/pqs_near_term_final_basis_realization_plan.md`
  - `docs/src/developer/pqs_source_box_fixture_policy.md`
  - `docs/src/developer/manager_doer_collaboration_contract_2026-05-26.md`

Current one-center complete-core/shell payload dataflow:
1. `cartesian_parent` can materialize one-center route objects through `_cartesian_one_center_parent_basis_object` and `_cartesian_one_center_parent_axis_bundle_object`; the parent then carries `parent.parent_axis_bundle_object`.
2. For one-center PQS, `_pqs_source_box_route_driver_shell_stage_low_order_shellization` selects terminal shellification. `cartesian_units` builds a `terminal_route_state` carrying the `CartesianShellification.ShellificationPlan` and `CartesianTerminalLowering.TerminalLoweringPlan`.
3. `_pqs_source_box_route_driver_complete_core_shell_source_plan_payload(parent, transforms, recipe)` consumes `transforms.terminal_route_state.shellification_plan`, `transforms.terminal_route_state.lowering_plan`, and `parent.parent_axis_bundle_object`. It creates:
   - `coulomb_gaussian_expansion(doacc=false)`;
   - `pqs_multilayer_shell_region_plan(shellification_plan, lowering_plan)`;
   - `pqs_multilayer_shell_source_plan(parent_axis_bundle_object, region_plan; bond_axis, term_coefficients)`.
4. `_pqs_source_box_route_driver_complete_core_shell_h1_payload` consumes that source-plan payload, parent center records, and axis layers from the parent axis bundle. It creates:
   - `pqs_multilayer_complete_core_shell_final_basis(source_plan)`;
   - `pqs_multilayer_complete_core_shell_h1_payload(source_plan; final_basis, coulomb_expansion, center_records, axis_layers)`.
5. `_pqs_source_box_route_driver_complete_core_shell_density_inputs` derives diagnostic support-density inputs from `source_plan.bundles` and the Coulomb expansion through `_pqs_source_box_ida_factor_provenance`.
6. `_pqs_source_box_route_driver_complete_core_shell_h1_j_diagnostic_payload` consumes source plan, final basis, H1 payload, axis weights, raw pair-factor terms, and Coulomb expansion to create the private H1/J diagnostic payload.
7. `_pqs_source_box_route_driver_complete_core_shell_ham_payload` consumes the above route-owned objects and exposes the private Ham payload or a structured blocker.
8. `_PQSCompleteCoreShellDiagnosticRoutePayload` is the current route-owned carrier assembled by `_pqs_source_box_route_driver_complete_core_shell_diagnostic_route_payload` and consumed by `cartesian_assembly`.

Structurally one-center-only pieces:
- `_cartesian_one_center_parent_basis_object` and `_cartesian_one_center_parent_axis_bundle_object` explicitly return not-applicable unless `system_classification == :one_center`.
- The complete-core/shell source-plan producer is structurally terminal-shellification backed: it requires a terminal `ShellificationPlan`, terminal `TerminalLoweringPlan`, and a parent axis-bundle object. It is not keyed from Be2 route-skeleton units.
- The selected one-center source-plan shape is `:shellification_backed_repeated_one_cell_projected_q_shell_layers`: one direct core plus ordered complete-shell layers collapsed into one shell sector.
- `pqs_complete_core_shell_final_basis` is conceptually a direct-core plus surrounding-shell final-basis helper with `support_row_order = :core_then_shell`. It does not know about the Be2 `(:pqs_left, :product, :pqs_right)` unit route shape.

Current Be2/PQS available structured objects:
- `cartesian_parent` classifies the Be2 fixture as `:bond_aligned_diatomic`, with `bond_axis = :x`, center table, nuclear charges, atom locations, parent axis counts, route-axis counts, and materialization planning metadata.
- The default Be2 PQS report/test path keeps `probe_parent_axis_construction = false`, so it has parent axis metadata but no `parent_axis_bundle_object` handoff.
- `_pqs_source_box_route_driver_route_skeleton` delegates to `_pqs_pqs_product_source_box_route_skeleton`, which carries structured Be2/PQS source-box data:
  - route shape `(:pqs_left, :product, :pqs_right)`;
  - source boxes for `pqs_left`, `product`, `pqs_right`;
  - source dimensions and retained counts/ranges;
  - retained units with retained-rule provenance and weight semantics;
  - pair entries and pair-family counts;
  - helper metadata for source-box density-density pairs.
- The route-driver report tests already assert the Be2 route construction vocabulary and materialization request metadata, including `:bond_aligned_diatomic_shellization`, `_nested_bond_aligned_diatomic_source`, and the pending materializer input blocker.
- What Be2/PQS does not carry today is a `pqs_multilayer_shell_region_plan`, a `pqs_multilayer_shell_source_plan`, a complete-core/shell final basis, H1 payload, density inputs, or H1/J density interaction. Pass 133's fingerprint records exactly that blocker.

True missing seam:
- The missing seam is best described as a Be2/PQS route-owned source/final-basis readiness or producer boundary, not a scalar report-field expansion.
- Generalizing `_pqs_source_box_route_driver_complete_core_shell_source_plan_payload` directly is not the smallest correct step because it assumes terminal shellification/lowering plus one parent axis-bundle, then produces a single complete core/shell source plan.
- Adding a thin adapter from Be2 route skeleton fields into the existing one-center complete-core/shell producer would be unsafe unless it first creates a real route-owned object with the same support/order semantics as `pqs_multilayer_shell_source_plan`.
- The Be2 route skeleton is source-box-first and structured, but its semantic shape is multi-unit diatomic `pqs_left/product/pqs_right` plus pair inventory. That is not yet equivalent to one direct core plus collapsed surrounding shell sector.
- The correct seam should therefore sit before final-basis/H1 materialization as a compact Be2/PQS complete-core/shell readiness/source-plan payload. It should inventory which structured Be2 facts are available and which producer objects are still missing, without promoting shell/support-row or WL adapter paths.

Recommended smallest next implementation pass:
- Add a private internal readiness payload, for example `_PQSBe2CompleteCoreShellHamReadinessPayload` or `_pqs_source_box_route_driver_be2_complete_core_shell_ham_readiness_payload`.
- Wire it only into the private diagnostic/assembly surface, next to the existing `complete_core_shell_diagnostic_route_payload`.
- Fields should be compact and object/fingerprint oriented:
  - `status`
  - `blocker`
  - `route_family`
  - `system_classification`
  - `bond_axis`
  - `center_summary`
  - `parent_axis_bundle_object_available`
  - `route_skeleton_summary`
  - `source_box_summary`
  - `retained_unit_summary`
  - `pair_inventory_summary`
  - `available_objects`
  - `missing_objects`
  - `summary`
  - `metadata`
- It should consume existing parent and route-skeleton objects, not report aliases, and should preserve the current Ham payload blocker. It should not build final basis, H1, H1/J, RHF, WL payloads, artifacts, exports, or checkpointing.
- After that payload exists, a later pass can choose whether to implement a true diatomic/PQS source-plan producer. That later producer should be reviewed as a route-owned source-plan object, not as a compatibility shim over one-center shellification.

Recommended focused validation for that pass:
- Extend or replace the pass-133 focused standalone Be2 fingerprint test with one compact route-owned readiness-payload assertion.
- Validate only:
  - Be2 route reaches `cartesian_assembly`;
  - readiness payload exists and is private;
  - status/blocker indicate missing Be2 complete-core/shell source-plan producer;
  - source-box, retained-unit, pair-family, center, and parent-axis availability fingerprints match the route skeleton;
  - current `complete_core_shell_ham_payload` blocker remains behavior-preserving;
  - no public API/export/artifact/RHF claims are set.
- Avoid broad `pqs_source_box_route_driver_report_runtests.jl` pressure and avoid adding scalar report aliases. The existing report tests can remain compatibility coverage until a later shrink pass.

Validation commands/results:
- Read-only inspection only: `sed`, `rg`, and `git status --short --branch`.
- No Julia command was run; no probe was needed to answer the seam question.
- `git status --short --branch`:
  - `## main...origin/main`

Deletion/shrinkage forecast:
- deleted: none in this no-edit audit.
- simplified: the next readiness payload can shrink the pass-133 test from checking scattered downstream blocker facts to checking one route-owned readiness object plus the preserved Ham blocker.
- quarantined: shell/support-row contraction, WL materializers/adapters, H1/J density interaction, RHF, exports, artifacts, and report scalar aliases remain outside Be2 route authority.
- not deleted because: current report tests still preserve existing route-construction vocabulary and materialization metadata compatibility.
- exact remaining caller/blocker: Be2 PQS has no route-owned `pqs_multilayer_shell_region_plan`/`pqs_multilayer_shell_source_plan` producer and no complete-core/shell final-basis/H1/density-input producer for the `(:pqs_left, :product, :pqs_right)` diatomic route skeleton. The current complete-core/shell Ham payload therefore correctly blocks on missing final basis, H1 payload, final one-body Hamiltonian, density inputs, H1/J diagnostic payload, and pre-final density interaction.

-- repo-doer@macmini
