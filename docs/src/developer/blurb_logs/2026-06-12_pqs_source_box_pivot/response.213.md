Pass 213 response - H2 physical gausslet H1 seam

Result:
- H1 materialized for the H2 physical q5 gausslet-only final basis.
- `route/h1_status = :materialized_pqs_physical_gausslet_h1_solve`
- `route/h1_materialized = true`
- `physics/endpoint_blocker = :missing_physical_gausslet_h1_j_builder`

Implemented:
- Added private `_PQSDiatomicPhysicalGaussletH1Payload`.
- The payload consumes:
  - `_PQSDiatomicPhysicalGaussletCoreShellSourcePlan`
  - `_PQSDiatomicPhysicalGaussletFinalBasisPayload`
- The payload materializes the physical support/final one-electron path in support order:
  - `:atom_contact_core`
  - `:shared_shell_1`
  - `:shared_shell_2`
- Wired `cartesian_assembly` to build and carry the H2 physical H1 payload.
- Extended the compact physical route report/artifact fields for H1 status, materialized flag, lowest H1 energy, finiteness, and symmetry error.
- Turned on `run_h1 = true` for `test/driver_inputs/h2_pqs_q5_physical_gausslet_r4.jl`.
- Updated the H2 463 driver artifact endpoint test to assert the H1 route facts and finite negative H1 sanity check.
- Updated the endpoint manifest blocker from missing H1 to missing H1-J.

Observed artifact facts:
- final dimension: `463`
- H1 lowest energy: `-0.7946609179724462`
- H1 Hamiltonian finite: `true`
- H1 Hamiltonian symmetry error: `7.105427357601002e-15`

Support/final one-body statuses:
- support kinetic: `:materialized_pqs_physical_gausslet_support_kinetic_matrix`
- support electron-nuclear by-center: `:materialized_pqs_physical_gausslet_support_electron_nuclear_by_center_matrix_set`
- final kinetic: `:materialized_pqs_physical_gausslet_final_one_body_matrix`
- final electron-nuclear by-center: `:materialized_pqs_physical_gausslet_final_electron_nuclear_by_center`
- final H1 Hamiltonian: `:materialized_pqs_physical_gausslet_final_one_electron_hamiltonian`
- H1 solve: `:materialized_pqs_physical_gausslet_h1_solve`

Files changed:
- `src/pqs_source_box_diatomic_complete_core_shell.jl`
- `src/pqs_source_box_route_driver_helpers.jl`
- `src/pqs_source_box_route_driver_reporting.jl`
- `test/driver_inputs/h2_pqs_q5_physical_gausslet_r4.jl`
- `test/nested/cartesian_ham_builder_h2_pqs_q5_physical_gausslet_r4_target_runtests.jl`
- `test/nested/integration_runtests.jl`
- `docs/src/developer/cartesian_driver_endpoint_manifest.md`
- deleted `test/nested/bond_aligned_diatomic_endcap_panel_shared_shell_source_policy_runtests.jl`

Artifact fields changed:
- `route/h1_status`
- `route/h1_materialized`
- `physics/h1_lowest`
- `physics/h1_hamiltonian_matrix_finite`
- `physics/h1_hamiltonian_symmetry_error`
- `physics/endpoint_blocker` now reports `:missing_physical_gausslet_h1_j_builder` for this target.

Source/test/bin scoped line budget:
- added: `556`
- deleted: `853`
- net: `-297`

Validation:
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - passed; precompile/load output ended with `load ok`
  - elapsed reported by Julia precompile output: `56439.5 ms`
- `julia --project=. -e 'using Test; t = @elapsed include("test/nested/cartesian_ham_builder_h2_pqs_q5_physical_gausslet_r4_target_runtests.jl"); println("elapsed_s=", t)'`
  - passed, `43/43`
  - `elapsed_s=77.847293917`
  - Previous pass was about 74-76 seconds. This grew only slightly while adding H1 support/final matrix materialization; assembly still reports high compilation share.
- focused artifact scalar read:
  - `h1_lowest=-0.7946609179724462`
  - `h1_finite=true`
  - `h1_symmetry_error=7.105427357601002e-15`
  - `route_h1_status=materialized_pqs_physical_gausslet_h1_solve`
  - `route_h1_materialized=true`
  - `physics_endpoint_blocker=missing_physical_gausslet_h1_j_builder`
  - `final_dimension=463`
- `git diff --check`
  - passed
- `git diff --cached --check`
  - passed
- `git diff --numstat -- src test bin tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl`
  - confirms net `-297` in scoped files.

Deletion/shrinkage report:
- deleted:
  - `test/nested/bond_aligned_diatomic_endcap_panel_shared_shell_source_policy_runtests.jl`
- simplified:
  - `test/nested/integration_runtests.jl` no longer includes that stale policy test.
  - H2 physical endpoint artifact now carries compact H1 facts instead of stopping at final-basis availability.
- quarantined:
  - H1 materialization remains private/internal on the H2 physical diagnostic route.
  - No H1-J, density, RHF, supplement, CR2, export, or public API work was added.
- not deleted because:
  - The H2 driver endpoint test remains the active physics seam guard.
  - The source-plan, final-basis, and H1 payloads are the current route-owned objects needed before the H1-J seam.
- exact remaining caller/blocker:
  - `physics/endpoint_blocker = :missing_physical_gausslet_h1_j_builder`

Git status at response time:
- dirty worktree with the pass-213 source/test/doc edits above.
- no commits made.

-- repo-doer@macmini
