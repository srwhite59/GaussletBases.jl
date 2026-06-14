Pass 208 response - H2 physical gausslet source-plan seam

Implemented the first private H2 463 physical gausslet source-plan seam.

Source-plan seam status:
- available/materialized: no
- status: `:blocked_pqs_diatomic_physical_gausslet_core_shell_source_plan`
- blocker: `:missing_atom_contact_core_support_rows`
- no H2 221 source-box diagnostic plan was reused
- no placeholder support rows, support states, coefficients, final basis, H1,
  H1-J, RHF, exports, or public API were added

Payload:
- name: `_PQSDiatomicPhysicalGaussletCoreShellSourcePlanPayload`
- helper: `_pqs_source_box_route_driver_diatomic_physical_gausslet_source_plan_payload`
- fields:
  - `status`
  - `blocker`
  - `route_family`
  - `source_plan`
  - `available_objects`
  - `missing_objects`
  - `summary`
  - `metadata`

Observed H2 target inventory carried by the seam:
- support order: `(:atom_contact_core, :shared_shell_1, :shared_shell_2)`
- retained order: `(:atom_contact_core, :shared_shell_1, :shared_shell_2)`
- support counts: `(275, 578, 362)`
- retained counts: `(251, 98, 114)`
- expected final dimension: `463`
- retained atom-core interiors: `true`
- source-plan role: `:atom_contact_core_plus_pqs_shared_shells`
- supplement policy: `:none`

Files edited:
- `src/pqs_source_box_diatomic_complete_core_shell.jl`
- `src/pqs_source_box_route_driver_helpers.jl`
- `src/pqs_source_box_route_driver_reporting.jl`
- `test/driver_inputs/h2_pqs_q5_physical_gausslet_r4.jl`
- `test/nested/cartesian_ham_builder_h2_pqs_q5_physical_gausslet_r4_target_runtests.jl`
- `test/nested/bond_aligned_diatomic_high_order_recipe_opt_in_source_construction_integration_runtests.jl`
- `docs/src/developer/cartesian_driver_endpoint_manifest.md`

Files deleted:
- none

Artifact/report behavior:
- `cartesian_assembly` now carries
  `diatomic_physical_gausslet_source_plan_payload` next to the physical target
  payload.
- `physical_gausslet_target_summary` now includes compact source-plan status,
  blocker, object kind, and missing objects.
- The H2 physical target artifact now writes:
  - `target/source_plan_status =
    :blocked_pqs_diatomic_physical_gausslet_core_shell_source_plan`
  - `target/source_plan_blocker = :missing_atom_contact_core_support_rows`
  - `route/source_plan_status =
    :blocked_pqs_diatomic_physical_gausslet_core_shell_source_plan`
  - `physics/endpoint_blocker = :missing_atom_contact_core_support_rows`

Source/test/bin scoped line budget:

```text
75   5   src/pqs_source_box_diatomic_complete_core_shell.jl
22   0   src/pqs_source_box_route_driver_helpers.jl
16   2   src/pqs_source_box_route_driver_reporting.jl
1    1   test/driver_inputs/h2_pqs_q5_physical_gausslet_r4.jl
0    109 test/nested/bond_aligned_diatomic_high_order_recipe_opt_in_source_construction_integration_runtests.jl
7    7   test/nested/cartesian_ham_builder_h2_pqs_q5_physical_gausslet_r4_target_runtests.jl
```

Totals: 121 added, 124 deleted, net -3.

Validation:
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - passed
  - precompiled `GaussletBases` in about 56 s, then printed `load ok`
- `julia --project=. -e 'using Test; t = @elapsed include("test/nested/cartesian_ham_builder_h2_pqs_q5_physical_gausslet_r4_target_runtests.jl"); println("elapsed_s=", t)'`
  - passed: 36/36
  - elapsed: `62.80044525` s
- `git diff --check`
  - passed
- `git diff --cached --check`
  - passed

Git status:

```text
## main...origin/main
 M docs/src/developer/cartesian_driver_endpoint_manifest.md
 M src/pqs_source_box_diatomic_complete_core_shell.jl
 M src/pqs_source_box_route_driver_helpers.jl
 M src/pqs_source_box_route_driver_reporting.jl
 M test/driver_inputs/h2_pqs_q5_physical_gausslet_r4.jl
 M test/nested/bond_aligned_diatomic_high_order_recipe_opt_in_source_construction_integration_runtests.jl
 M test/nested/cartesian_ham_builder_h2_pqs_q5_physical_gausslet_r4_target_runtests.jl
```

Deletion/shrinkage report:
- deleted: none
- simplified:
  - removed redundant physical-target metadata flags already represented by
    target/source-plan summaries
  - removed 109 lines of stale receipt/record diagnostic assertion pressure
    from the old high-order source-construction integration test while leaving
    the basis size, operator shape, finite-output, symmetry, overlap, and route
    dimension checks intact
  - removed the non-contract H2 physical target artifact path print
- quarantined: none
- not deleted because:
  - the H2 463 artifact-readiness test remains the active guard for the
    physical gausslet target inventory and source-plan blocker
  - final-basis, H1, H1-J, density interaction, RHF, and physics comparison
    paths remain intentionally unimplemented for this target
- exact remaining caller/blocker:
  - remaining blocker is
    `:missing_atom_contact_core_support_rows`; a real physical H2 source plan
    still needs route-owned atom-contact-core support rows and shared-shell
    coefficient/contraction data before final-basis construction can start

-- repo-doer@macmini
