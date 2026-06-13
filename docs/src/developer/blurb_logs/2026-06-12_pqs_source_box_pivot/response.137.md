Pass 137 response - diatomic PQS source-plan producer contract audit

Scope:
- No source edits.
- No test edits.
- No Julia probes or tests run; read-only inspection was sufficient.
- No commit.

Files/helpers inspected:
- `src/pqs_source_box_route_driver_helpers.jl`
  - `_PQSDiatomicCompleteCoreShellHamReadinessPayload`
  - `_pqs_source_box_route_driver_diatomic_complete_core_shell_ham_readiness_payload`
  - `_pqs_source_box_route_driver_complete_core_shell_source_plan_payload`
  - `_pqs_source_box_route_driver_complete_core_shell_h1_payload`
  - `_pqs_source_box_route_driver_complete_core_shell_density_inputs`
  - `_pqs_source_box_route_driver_complete_core_shell_h1_j_diagnostic_payload`
  - `_pqs_source_box_route_driver_complete_core_shell_ham_payload`
  - `_pqs_source_box_route_driver_complete_core_shell_diagnostic_route_payload`
- `src/pqs_source_box_route_driver_skeletons.jl`
  - `_pqs_source_box_route_driver_route_skeleton`
- `src/CartesianContractedParentMetrics.jl`
  - `_pqs_pqs_product_source_box_route_skeleton`
  - raw product box probe helper surfaces around `_pqs_explicit_core_spacing_route_raw_product_box_plan_probe`
- `src/pqs_multilayer_shell_region_plan.jl`
  - `PQSMultilayerShellRegionPlan`
  - `pqs_multilayer_shell_region_plan`
- `src/pqs_multilayer_shell_source_plan.jl`
  - `pqs_multilayer_shell_source_plan`
  - `_pqs_multilayer_realize_shell_source_plan`
- `src/pqs_multilayer_complete_core_shell_h1.jl`
  - `pqs_multilayer_complete_core_shell_final_basis`
  - `pqs_multilayer_complete_core_shell_h1_payload`
  - `pqs_multilayer_complete_core_shell_h1_j_payload`
- `src/pqs_multilayer_support_one_body.jl`
  - `pqs_multilayer_support_kinetic_matrix`
  - `pqs_multilayer_support_electron_nuclear_by_center_matrices`
- `src/pqs_multilayer_support_density.jl`
  - `pqs_multilayer_support_weights`
  - `pqs_multilayer_support_pair_raw_numerator_matrix`
- `src/cartesian_final_basis_realization/pqs_complete_core_shell_final_basis.jl`
  - `pqs_complete_core_shell_final_basis`
  - `pqs_complete_core_shell_pre_final_density_interaction`
- `test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl`
- Governing docs:
  - `docs/src/developer/pqs_source_box_operator_framework.md`
  - `docs/src/developer/pqs_near_term_final_basis_realization_plan.md`
  - `docs/src/developer/pqs_source_box_fixture_policy.md`

Existing source-plan fields and consumers:

`pqs_multilayer_complete_core_shell_final_basis(plan)` requires:
- `object_kind == :pqs_multilayer_shell_source_plan`
- `status == :available_pqs_multilayer_shell_source_plan`
- `metrics.x/y/z.overlap`
- `core_support_states`
- `shell_support_states`
- `core_support_indices`
- `shell_support_indices`
- `shell_final_coefficients`
- `metadata`
- summary fields copied after materialization: `layer_count`, support counts, status.

`pqs_multilayer_complete_core_shell_h1_payload(plan; final_basis, ...)` requires:
- same `object_kind` and available status;
- `pqs_multilayer_support_kinetic_matrix(plan)`, which consumes:
  - `core_support_states`
  - `shell_support_states`
  - `metrics.x/y/z.overlap`
  - `metrics.x/y/z.kinetic`
- `pqs_multilayer_support_electron_nuclear_by_center_matrices(plan; ...)`, which consumes:
  - `core_support_states`
  - `shell_support_states`
  - `coulomb_expansion.coefficients`
  - center records
  - either axis layers or explicit Gaussian factor terms by center.

Density-input/H1-J helpers require:
- `_pqs_source_box_route_driver_complete_core_shell_density_inputs(source_plan, coulomb_expansion)`:
  - `object_kind == :pqs_multilayer_shell_source_plan`
  - `status == :available_pqs_multilayer_shell_source_plan`
  - `bundles`
  - `coulomb_expansion.coefficients`
  - then derives `axis_weights` and `raw_pair_factor_terms` through `_pqs_source_box_ida_factor_provenance(bundles; expected_term_count)`.
- `pqs_multilayer_support_weights(plan; axis_weights)`:
  - `core_support_states`
  - `shell_support_states`
  - support order is core first, shell second.
- `pqs_multilayer_support_pair_raw_numerator_matrix(plan; raw_pair_factor_terms, coulomb_expansion)`:
  - `core_support_states`
  - `shell_support_states`
  - `coulomb_expansion.coefficients`
  - term-first raw pair-factor arrays.
- `pqs_complete_core_shell_pre_final_density_interaction(final_basis, support_pair_raw, support_weights)`:
  - final basis in `core_then_shell` support order
  - support weights and pair matrix in that same order.

Universal route concepts versus one-center-specific fields:

Universal concepts a diatomic producer can provide:
- source-plan object/status/blocker;
- parent axis bundles and axis metrics;
- Coulomb expansion coefficients/factor provenance;
- center records and center locations/charges;
- support row indices/states in one declared order;
- partition of support rows into direct/core-like rows and retained/source-realized rows, if the final-basis helper will stay complete-core/shell shaped;
- retained coefficients for the source-realized sector;
- support overlap/kinetic/electron-nuclear/density factor data derived from source-box route authority;
- compact summary/fingerprint and non-promotion metadata.

One-center or current shellification-specific concepts:
- `PQSMultilayerShellRegionPlan` itself: terminal `ShellificationPlan` plus `TerminalLoweringPlan`.
- `shell_layers`, `terminal_region_key`, `lowering_kind`, and `source_cpbs` from terminal shellification/lowering.
- `source_kind = :shellification_backed_repeated_one_cell_projected_q_shell_layers`.
- "one direct core plus surrounding complete-shell layers" geometry.
- `shell_final_coefficients` currently built by `_nested_projected_q_shell_layer` followed by `_pqs_shell_realization_plan`, which is shell-realization/Lowdin machinery.
- `layer_count` is only universal if reinterpreted as source-realized sector count; it is not naturally the Be2 `pqs_left/product/pqs_right` unit count.

Existing Be2/PQS route-skeleton facts relevant to a producer:
- `route_family = :pqs_source_box`
- `route_shape = (:pqs_left, :product, :pqs_right)`
- `parent_axis_counts`
- `q`
- source boxes:
  - `pqs_left`
  - `product`
  - `pqs_right`
- source dimensions for each unit.
- retained units with:
  - `unit_key`
  - `unit_role`
  - `retained_unit_kind`
  - `source_family`
  - `source_box`
  - `source_dimensions`
  - `source_dimension`
  - `retained_rule_kind`
  - `retained_rule_derivation`
  - `retained_range`
  - `retained_count`
  - `weight_semantics`
- retained dimension/ranges.
- pair entries and pair-family counts.
- source-box density-density helper metadata.
- diagnostics that keep PQS source-box-first and mark retained PQS weights as non-quadrature.

Existing parent/route facts relevant to a producer:
- with `probe_parent_axis_construction = :auto`, Be2 readiness has `parent_axis_bundle_object_available = true`.
- center metadata: `:bond_aligned_diatomic`, `bond_axis = :x`, two centers, charges/locations.
- a Coulomb expansion can be constructed by the same private route code path, but it is not yet attached to a diatomic source-plan producer.
- raw product box plan probe machinery exists, but the focused Be2 fixture currently has `probe_raw_product_box_plans = false`; it should not become a producer dependency unless a later pass explicitly chooses that route-owned object.

Recommended source-plan producer contract:

Add a private payload/helper, not a public API:
- `_PQSDiatomicCompleteCoreShellSourcePlanPayload`
- `_pqs_source_box_route_driver_diatomic_complete_core_shell_source_plan_payload(parent, route_skeleton, recipe; source_plan_mode = :readiness_only)`

Required input objects:
- `parent`
  - `system_classification`
  - `bond_axis`
  - `center_table`
  - `parent_axis_bundle_object`
  - parent axis counts/readiness metadata
- `route_skeleton`
  - source boxes/dimensions
  - retained units/counts/ranges
  - retained dimension
  - pair entries/family counts
  - route shape
- `recipe`
  - route family/kind
  - retained rules
  - pair factor normalization
- Coulomb expansion source, initially local `coulomb_gaussian_expansion(doacc=false)` or an explicit argument if manager wants injection.

Output fields for the first pass:
- `status`
- `blocker`
- `route_family`
- `system_classification`
- `bond_axis`
- `parent_axis_bundle_object_available`
- `route_skeleton_summary`
- `source_box_summary`
- `retained_unit_summary`
- `pair_inventory_summary`
- `center_summary`
- `coulomb_expansion_summary`
- `source_plan` initially `nothing`
- `source_plan_status`
- `available_objects`
- `missing_objects`
- `summary`
- `metadata`

Explicit blockers for the first pass:
- not PQS source-box route
- not bond-aligned diatomic
- missing parent axis-bundle object
- missing center metadata
- missing route skeleton source boxes/retained units/pair inventory
- `:missing_diatomic_complete_core_shell_source_plan_materializer` or `:missing_diatomic_complete_core_shell_source_realization_contract`

The first implementation should not materialize final basis, H1, H1-J, or Ham. I recommend it also not materialize a fake `pqs_multilayer_shell_source_plan`. It should return a structured blocked source-plan payload that narrows the blocker to the actual materializer/realization contract. That is the smallest honest pass and avoids forcing Be2 into the one-center terminal-shellification producer.

Next object-shape decision:
- Longer term, the producer should probably return an object with `object_kind = :pqs_multilayer_shell_source_plan` only if it truly satisfies the existing duck-typed consumer shape:
  - `bundles`
  - `metrics`
  - `core_support_indices`
  - `core_support_states`
  - `shell_support_indices`
  - `shell_support_states`
  - `shell_final_coefficients`
  - `metadata`
  - available status.
- If Be2 cannot honestly express the `core_then_shell` support order, then do not force it. Add a later final-basis adapter for a diatomic source-plan object instead.
- For the next pass, keep the object blocked/readiness-only so the design does not accidentally promote shell/support-row contraction or WL old-kernel adapter authority.

Smallest next implementation pass:
1. Add `_PQSDiatomicCompleteCoreShellSourcePlanPayload` plus helper in `src/pqs_source_box_route_driver_helpers.jl`.
2. Call it from `_pqs_source_box_route_driver_diatomic_complete_core_shell_ham_readiness_payload` or adjacent assembly wiring.
3. In the probe-enabled Be2 fixture, assert:
   - parent axis-bundle is available;
   - source-plan payload exists;
   - source-plan payload blocks on `:missing_diatomic_complete_core_shell_source_realization_contract`;
   - existing Ham readiness still blocks on `:missing_diatomic_complete_core_shell_source_plan_producer` or is renamed to the narrower source-plan materializer blocker if manager approves;
   - no final-basis/H1/H1-J/Ham materialization flags are set.
4. Do not build a final basis and do not run H1.

Focused validation:
- `julia --project=. test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git diff --check`

Test pressure to avoid:
- Do not broaden `test/nested/pqs_source_box_route_driver_report_runtests.jl`.
- Do not add a physics endpoint.
- Do not add report alias assertions.
- Do not run broad nested suites.

Git status:
```text
## main...origin/main
```

Deletion/shrinkage forecast:
- deleted: none in this no-edit audit.
- simplified: next pass can move Be2 source-plan readiness out of scattered Ham-readiness missing-object lists into a route-owned source-plan payload.
- quarantined: shell/support-row contraction, WL adapter payloads, final-basis/H1/H1-J materialization, RHF/SCF, public API, exports, artifacts, hfdmrg, and CR2 remain out of scope.
- not deleted because: current Be2 fingerprint tests and readiness payload are still the only compact validation that separates parent axis-bundle availability from source-plan producer availability.
- exact remaining caller/blocker: `cartesian_assembly` has a Be2 readiness payload and, with probe enabled, a parent axis-bundle object, but there is no route-owned diatomic complete-core/shell source-plan payload/materializer that can supply the existing final-basis/H1/H1-J consumer shape. The remaining blocker is the source-plan materializer/realization contract, not parent axis-bundle construction.

-- repo-doer@macmini
