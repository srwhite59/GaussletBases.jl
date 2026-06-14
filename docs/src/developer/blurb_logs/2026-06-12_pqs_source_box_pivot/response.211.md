Pass 211 response - private H2 physical source-plan object

Implemented a private route-owned source-plan object:

- object: `_PQSDiatomicPhysicalGaussletCoreShellSourcePlan`
- object kind: `:pqs_diatomic_physical_gausslet_core_shell_source_plan`
- source-plan status: `:available_pqs_diatomic_physical_gausslet_core_shell_source_plan`
- next blocker: `:missing_physical_gausslet_final_basis_builder`
- authority status advanced: yes, to `:private_source_backed_adapter_authority`

The object is built only from the checked source-backed candidate and carries:

- parent basis and axis bundles
- atom-contact core support indices/states
- shared-shell support indices/states
- core coefficient matrix slice
- shared-shell coefficient matrix slices
- support/retained order `(:atom_contact_core, :shared_shell_1, :shared_shell_2)`
- retained ranges and final dimension
- convention labels for private source-backed route authority, no supplement,
  and no final-basis/H1/H1-J/RHF/export/public materialization

Counts:

- support counts: `(275, 578, 362)`
- retained counts: `(251, 98, 114)`
- final dimension: `463`

Artifact/report fields changed:

- `target/source_plan_status` now reports
  `:available_pqs_diatomic_physical_gausslet_core_shell_source_plan`
- `target/source_plan_blocker` now reports
  `:missing_physical_gausslet_final_basis_builder`
- `target/source_plan_authority_status` now reports
  `:private_source_backed_adapter_authority`
- `route/source_plan_status` now reports
  `:available_pqs_diatomic_physical_gausslet_core_shell_source_plan`
- `physics/endpoint_blocker` now reports
  `:missing_physical_gausslet_final_basis_builder`
- `physics/endpoint_ready` remains `false`
- endpoint manifest row updated to the new blocker

Files changed:

- `src/pqs_source_box_diatomic_complete_core_shell.jl`
- `test/driver_inputs/h2_pqs_q5_physical_gausslet_r4.jl`
- `test/nested/cartesian_ham_builder_h2_pqs_q5_physical_gausslet_r4_target_runtests.jl`
- `test/nested/bond_aligned_diatomic_high_order_recipe_opt_in_source_construction_integration_runtests.jl`
- `docs/src/developer/cartesian_driver_endpoint_manifest.md`

Source/test/bin scoped line budget:

```text
128  8    src/pqs_source_box_diatomic_complete_core_shell.jl
1    1    test/driver_inputs/h2_pqs_q5_physical_gausslet_r4.jl
0    127  test/nested/bond_aligned_diatomic_high_order_recipe_opt_in_source_construction_integration_runtests.jl
5    5    test/nested/cartesian_ham_builder_h2_pqs_q5_physical_gausslet_r4_target_runtests.jl
```

Totals for scoped files: 134 added, 141 deleted, net -7.

Validation:

```sh
julia --project=. -e 'using GaussletBases; println("load ok")'
```

Result: passed, printed `load ok`.

```sh
julia --project=. -e 'using Test; t = @elapsed include("test/nested/cartesian_ham_builder_h2_pqs_q5_physical_gausslet_r4_target_runtests.jl"); println("elapsed_s=", t)'
```

Result: passed, 40/40, `elapsed_s=73.989651291`.

The runtime remains around the expected 75 seconds. The fixed-source candidate
phase printed `diatomic.fixed_source.total: 11.18 seconds`; total runtime is
still dominated by compilation/report construction rather than a new repeated
runtime cost.

```sh
git diff --check
git diff --cached --check
```

Result: both passed.

Git status:

```text
## main...origin/main
 M docs/src/developer/cartesian_driver_endpoint_manifest.md
 M src/pqs_source_box_diatomic_complete_core_shell.jl
 M test/driver_inputs/h2_pqs_q5_physical_gausslet_r4.jl
 M test/nested/bond_aligned_diatomic_high_order_recipe_opt_in_source_construction_integration_runtests.jl
 M test/nested/cartesian_ham_builder_h2_pqs_q5_physical_gausslet_r4_target_runtests.jl
```

Deletion/shrinkage report:

- deleted: none
- simplified: removed 127 stale diagnostic-vocabulary assertions from the old
  non-default high-order opt-in source-construction integration test; retained
  the active block-size, finite-output, oracle-comparison, dimension, and
  support-coverage checks
- quarantined: source-backed authority remains explicitly private and labeled
  as a source-backed adapter, with no public/export/final-basis/H1/H1-J/RHF
  materialization claims
- not deleted because: the H2 target artifact test is the active seam guard, and
  the source-backed candidate remains needed as the private adapter input until
  the physical final-basis builder exists
- exact remaining caller/blocker:
  `:missing_physical_gausslet_final_basis_builder`

-- repo-doer@macmini
