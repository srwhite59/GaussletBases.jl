# Pass 230 response - Independent H2 PQS source-box recovery audit

## Audit verdict

The independent PQS path has useful lower-level source-box machinery, but it is
not ready to materialize the physical H2 463 source plan. The current 463 route
still reaches its retained transform through the fake-PQS source-backed WL/QW
adapter. I did not find an independent PQS rule that generates the
`atom_contact_core` retained count `251`, and I did not find an independent PQS
rule that generates `shared_shell_2` retained count `114`.

The smallest honest next step is a target/readiness pass for a separate real H2
PQS route that records support vocabulary and blocks before final basis/H1.

## Source surfaces inspected

- `src/pqs_multilayer_shell_region_plan.jl`
- `src/pqs_multilayer_shell_source_plan.jl`
- `src/pqs_multilayer_support_one_body.jl`
- `src/pqs_multilayer_support_density.jl`
- `src/pqs_source_box_diatomic_complete_core_shell.jl`
- `src/pqs_source_box_route_driver_skeletons.jl`
- `src/cartesian_raw_product_sources/CartesianRawProductSources.jl`
- `src/cartesian_raw_product_sources/records.jl`
- `src/cartesian_raw_product_sources/source_mode_indices.jl`
- `src/cartesian_terminal_lowering/region_contracts.jl`
- `src/cartesian_terminal_lowering/selection.jl`
- live replacements for absent flat module files:
  - `src/cartesian_route_core/CartesianRouteCore.jl`
  - `src/cartesian_shellification/CartesianShellification.jl`
  - `src/cartesian_terminal_lowering/CartesianTerminalLowering.jl`
  - `src/cartesian_retained_units/CartesianRetainedUnits.jl`
  - `src/cartesian_pair_operator_plans/CartesianPairOperatorPlans.jl`
  - `src/cartesian_pair_block_materialization/CartesianPairBlockMaterialization.jl`
- history/guardrail docs:
  - `docs/src/developer/pqs_manager_running_log.md`
  - `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/summary.md`
  - `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.229.md`
  - selected earlier pass responses named by `rg`

## Answers

1. Existing independent PQS construction that is still usable:

   - `CartesianRawProductSources` owns raw product source-box facts:
     `RawProductBoxPlan`, source CPB, source mode dimensions/order, and
     `PQSBoundaryProductModeRetainedRule`.
   - The q=5 boundary retained rule is independent and produces the familiar
     `98` retained boundary product modes for a `5 x 5 x 5` source box.
   - `CartesianTerminalLowering._pqs_complete_shell_contract` builds a
     complete-shell lowering contract from `filled_cpb(raw.outer_box...)` with
     `:pqs_filled_source_cpb` and retained rule
     `:pqs_boundary_comx_product_modes`.
   - `pqs_multilayer_shell_region_plan` can organize one direct core plus
     complete PQS shell layers from shellification/lowering data.
   - `pqs_multilayer_shell_source_plan` can realize those shell layers with
     `_nested_projected_q_shell_layer`, shell projection, and Lowdin cleanup
     without using the fixed-source WL/QW oracle.
   - `CartesianPairBlockMaterialization` has source-space safe PQS one-body
     helpers and bridge/readiness summaries, but its facade says it does not
     build final shell-realized PQS pair blocks or Hamiltonian artifacts.

2. Old 221-dimensional H2 diagnostic route:

   Reusable:

   - source-box-first H2 vocabulary around `:pqs_left`, `:product`,
     `:pqs_right`;
   - raw product box/source-mode setup for `pqs_left` and `pqs_right`;
   - the independent q=5 boundary selection count `98`;
   - the guardrail that fixed-source oracle data is not source-plan authority.

   Diagnostic-only:

   - route shape `(:pqs_left, :product, :pqs_right)`;
   - support order `(:product, :pqs_left, :pqs_right)`;
   - middle/product retained `25` and total dimension `221`;
   - any final-basis/H1 interpretation of that 221 route as physical H2.

   The old tracked diagnostic input/test
   `test/driver_inputs/h2_pqs_q5_source_box_diagnostic_r4.jl` and
   `test/nested/cartesian_ham_builder_h2_pqs_q5_source_box_diagnostic_r4_runtests.jl`
   are no longer present; `response.222.md` records their deletion.

3. `atom_contact_core` construction without WL/QW coefficients:

   The current physical/fake inventory hard-codes support count `275` and
   retained count `251` from reviewed WL/QW data. The source-backed candidate
   path slices `source.sequence.coefficient_matrix` into
   `core_coefficient_matrix` and shared-shell matrices. That is not independent
   PQS construction.

   A real route should first construct the `atom_contact_core` support rows
   from route geometry, then define an independent retained rule for that
   support block. I did not find such a rule. Keeping all `275` rows as direct
   core support is possible as a different route choice, but it does not
   reproduce the current 463 contract and should not be silently normalized as
   the old `251`.

4. `shared_shell_1` and `shared_shell_2` from filled PQS source CPBs:

   - `shared_shell_1`: plausibly yes at the retained-rule level. The standard
     q=5 PQS boundary rule independently generates `98` modes, matching the
     current first shared-shell retained count. The missing part is the H2
     support-to-source-shell materializer for the physical support count `578`.
   - `shared_shell_2`: not with the current standard q=5 boundary rule. The
     fake route expects retained count `114`, while the existing boundary rule
     gives `98` for a cubic `5 x 5 x 5` source box. I found no independent PQS
     rule that explains `114`; it remains WL/QW-derived or at least unproven.

5. Real target route kind and driver input:

   Recommended route kind:

   ```julia
   :bond_aligned_diatomic_independent_pqs_source_box_core_shell
   ```

   Recommended driver input:

   ```text
   test/driver_inputs/h2_pqs_q5_independent_source_box_r4.jl
   ```

   Initial artifact/readiness fields should include:

   ```text
   fake_pqs/enabled = false
   source_backed_fixed_source_oracle_used = false
   retained_transform_authority = :pqs_source_box_construction
   physics/endpoint_ready = false
   physics/endpoint_blocker = :missing_independent_pqs_atom_contact_core_retained_rule
   ```

6. Exact blocker:

   Primary blocker:

   ```julia
   :missing_independent_pqs_atom_contact_core_retained_rule
   ```

   Secondary blocker:

   ```julia
   :missing_independent_pqs_shared_shell_2_retained_rule
   ```

   Broader route-level blocker:

   ```julia
   :missing_independent_pqs_physical_source_plan_materializer
   ```

## Recommended next pass

Create a blocked target/readiness surface only, not final basis/H1:

- add a separate real H2 PQS route kind and small driver input name;
- record common support vocabulary `(:atom_contact_core, :shared_shell_1,
  :shared_shell_2)` and support counts `(275, 578, 362)`;
- explicitly set `fake_pqs/enabled = false` and
  `source_backed_fixed_source_oracle_used = false`;
- block before source-plan materialization with
  `:missing_independent_pqs_atom_contact_core_retained_rule`;
- optionally report `shared_shell_1` as retained-rule-plausible and
  `shared_shell_2` as blocked on an independent retained rule.

Do not build final basis, H1, H1-J, RHF, supplements, CR2, exports, or public
API in that pass.

## Validation and status

Read-only inspection only. No Julia commands or tests were run.

Commands included `sed`, `nl`, `rg`, `rg --files`, and:

```text
git status --short --branch
## main...origin/main
```

No source, test, driver, or tracked design files were edited in this audit.

## Carrying-cost/deletion note

deleted: none.

simplified: none; read-only audit.

quarantined: none in code this pass; fake-PQS remains the existing quarantine.

not deleted because: the fake-PQS 463 route is still a useful golden regression
for WL/QW fixed-source reproduction, and deletion was out of scope.

exact remaining caller/blocker: the independent H2 PQS target cannot materialize
a source plan until `atom_contact_core` retained rows are generated by an
independent PQS rule rather than imported from the WL/QW fixed-source oracle;
`shared_shell_2` also lacks an independent retained rule for count `114`.

-- repo-doer@macmini
