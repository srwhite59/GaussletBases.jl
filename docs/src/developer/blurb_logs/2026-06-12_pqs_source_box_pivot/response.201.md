Pass 201 response - physical H2 target inventory surface.

Implemented only the physical H2 target/inventory readiness layer. No source-plan construction, final-basis construction, H1, H1-J, density interaction, or RHF was added.

Files changed:

- `src/pqs_source_box_route_driver_skeletons.jl`
- `src/pqs_source_box_diatomic_complete_core_shell.jl`
- `src/pqs_source_box_route_driver_helpers.jl`
- `src/pqs_source_box_route_driver_reporting.jl`
- `test/driver_inputs/h2_pqs_q5_physical_gausslet_r4.jl`
- `test/nested/cartesian_ham_builder_h2_pqs_q5_physical_gausslet_r4_target_runtests.jl`
- `test/nested/cartesian_pair_block_one_body_consumer_smoke_runtests.jl`
- deleted `test/nested/cartesian_pair_block_one_body_lw_plan_batch_runtests.jl`
- deleted `test/nested/cartesian_pair_block_one_body_accessors_contract_runtests.jl`

Target inventory source:

- The target inventory is hard-coded from the reviewed pass-200 WL/QW gausslet-only contract.
- It does not construct the old supplemented route and does not derive by building H1 or final-basis objects.

Target inventory numbers:

- route kind: `:bond_aligned_diatomic_physical_gausslet_core_shell_pqs`
- parent axis counts: `(x = 9, y = 9, z = 15)`
- support units/order: `(:atom_contact_core, :shared_shell_1, :shared_shell_2)`
- support counts: `(275, 578, 362)`
- retained units/order: `(:atom_contact_core, :shared_shell_1, :shared_shell_2)`
- retained counts: `(251, 98, 114)`
- expected final dimension: `463`
- retained atom-core interiors: `true`
- source plan role: `:atom_contact_core_plus_pqs_shared_shells`
- supplement policy: `:none`

What changed:

- Added a route-kind-specific metadata skeleton for the physical H2 target.
- Added private payload `_PQSDiatomicPhysicalGaussletCoreShellTargetPayload`.
- Added report fields `physical_gausslet_target_summary/status/blocker`.
- Added artifact group:
  - `target/status`
  - `target/blocker`
  - `target/support_units`
  - `target/support_counts`
  - `target/retained_units`
  - `target/retained_counts`
  - `target/retained_order`
  - `target/expected_final_dimension`
  - `target/retained_atom_core_interiors`
  - `target/source_plan_role`
  - `target/supplement_policy`
- Added driver input `test/driver_inputs/h2_pqs_q5_physical_gausslet_r4.jl`.
- Added focused artifact test `test/nested/cartesian_ham_builder_h2_pqs_q5_physical_gausslet_r4_target_runtests.jl`.
- Kept the existing 221-dimensional diagnostic route untouched.

Downstream status after this pass:

- target inventory: available
- source plan: intentionally not available
- final basis: intentionally not available
- H1: intentionally not materialized
- H1-J: intentionally not materialized
- private RHF: intentionally not materialized
- next blocker: `:missing_physical_gausslet_source_plan`

Validation:

- `julia --project=. test/nested/cartesian_ham_builder_h2_pqs_q5_physical_gausslet_r4_target_runtests.jl`
  - passed, 37 tests.
  - first run: 2m05.9s including about 58s package precompilation.
  - final rerun: 1m06.7s; still compilation-heavy in assembly.
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - passed, printed `load ok`.
- `git diff --check`
  - passed.
- Deleted-test live search:
  - `rg -n "cartesian_pair_block_one_body_lw_plan_batch_runtests|CartesianPairBlockMaterialization mixed White-Lindsey plan batch|CPBMMixedLWPlan|mixed_lw_plan_batch|cartesian_pair_block_one_body_accessors_contract_runtests|CartesianPairBlockMaterialization mixed one-body accessor/collection contract" test src bin`
  - no matches; exit code 1 as expected.
- `git diff --numstat -- src test bin tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl`
  - `171	0	src/pqs_source_box_diatomic_complete_core_shell.jl`
  - `39	0	src/pqs_source_box_route_driver_helpers.jl`
  - `50	0	src/pqs_source_box_route_driver_reporting.jl`
  - `128	0	src/pqs_source_box_route_driver_skeletons.jl`
  - `0	1438	test/nested/cartesian_pair_block_one_body_accessors_contract_runtests.jl`
  - `1	4	test/nested/cartesian_pair_block_one_body_consumer_smoke_runtests.jl`
  - `0	214	test/nested/cartesian_pair_block_one_body_lw_plan_batch_runtests.jl`
- New untracked source/test files:
  - `39 test/driver_inputs/h2_pqs_q5_physical_gausslet_r4.jl`
  - `81 test/nested/cartesian_ham_builder_h2_pqs_q5_physical_gausslet_r4_target_runtests.jl`
- Line budget including untracked new files:
  - added `509`
  - deleted `1656`
  - net `-1147`

Git status:

```text
## main...origin/main
 M src/pqs_source_box_diatomic_complete_core_shell.jl
 M src/pqs_source_box_route_driver_helpers.jl
 M src/pqs_source_box_route_driver_reporting.jl
 M src/pqs_source_box_route_driver_skeletons.jl
 D test/nested/cartesian_pair_block_one_body_accessors_contract_runtests.jl
 M test/nested/cartesian_pair_block_one_body_consumer_smoke_runtests.jl
 D test/nested/cartesian_pair_block_one_body_lw_plan_batch_runtests.jl
?? test/driver_inputs/h2_pqs_q5_physical_gausslet_r4.jl
?? test/nested/cartesian_ham_builder_h2_pqs_q5_physical_gausslet_r4_target_runtests.jl
```

Deletion/shrinkage report:

- deleted:
  - `test/nested/cartesian_pair_block_one_body_lw_plan_batch_runtests.jl`
  - `test/nested/cartesian_pair_block_one_body_accessors_contract_runtests.jl`
- simplified:
  - removed a stale comment in the mixed one-body smoke that pointed at the deleted accessor scaffold.
- quarantined:
  - none.
- not deleted because:
  - the remaining mixed one-body consumer smoke and LW record dispatch test still carry the compact live contract for direct/PQS/LW dispatch.
  - H2 diagnostic and physical target tests are current endpoint/target surfaces.
- exact remaining caller/blocker:
  - physical H2 PQS target now has route-owned inventory, but cannot build source plans, final basis, H1, H1-J, or RHF until a physical gausslet source-plan producer exists for `:atom_contact_core_plus_pqs_shared_shells`.

-- repo-doer@macmini
