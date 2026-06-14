Pass 197 response - H2 H1-only materialization.

Files changed:

- `test/driver_inputs/h2_pqs_q5_gausslet_only_r4.jl`
- `src/pqs_source_box_route_driver_helpers.jl`
- `src/pqs_source_box_route_driver_reporting.jl`
- `test/nested/cartesian_ham_builder_h2_pqs_q5_gausslet_only_r4_readiness_runtests.jl`
- deleted `test/nested/cartesian_pair_block_one_body_record_dispatch_runtests.jl`

H1 request surface used:

- The H2 driver input now requests:
  - `run_final_basis = true`
  - `run_h1 = true`
  - `run_h1_j = false`
  - `run_private_rhf = false`
- H1 uses the existing diatomic complete core/shell route path:
  - `_pqs_source_box_route_driver_diatomic_complete_core_shell_h1_payload`
- No second H1 path was added.

Artifact/report updates:

- The diatomic readiness artifact now records compact existing H1 scalars:
  - `physics/h1_lowest`
  - `physics/h1_hamiltonian_matrix_finite`
  - `physics/h1_hamiltonian_symmetry_error`
- These are read from the existing H1 payload/final Hamiltonian objects.
- No broad H1 matrix is saved.

H1 status and scalar values:

- `route/h1_status = :materialized_pqs_complete_core_shell_final_h1_solve`
- `route/h1_materialized = true`
- `physics/h1_lowest = 0.14582426982296057`
- `physics/h1_hamiltonian_matrix_finite = true`
- `physics/h1_hamiltonian_symmetry_error = 3.6637359812630166e-15`

Route readiness after this pass:

- `route/source_plan_status = :available_pqs_diatomic_complete_core_shell_source_plan`
- `route/final_basis_status = :available_pqs_complete_core_shell_final_basis`
- `basis/final_dimension = 221`
- `basis/final_overlap_identity_error = 2.510171831978689e-13`
- `route/h1_j_materialized = false`
- `private_rhf/materialized = false`
- `comparison/ready = false`
- Actual next blocker preserved:
  - `route/readiness_blocker = :missing_diatomic_complete_core_shell_hamiltonian_handoff_payload`
- Readback also showed:
  - `route/ham_input_status = :available_diatomic_complete_core_shell_ham_input_payload`
  - `route/hamiltonian_handoff_status = :blocked_diatomic_complete_core_shell_hamiltonian_handoff_payload`

Validation and timing:

- `julia --project=. -e 't = @elapsed include("test/nested/cartesian_ham_builder_h2_pqs_q5_gausslet_only_r4_readiness_runtests.jl"); println("h2_h1_readiness_elapsed_s=", t)'`
  - passed after updating the test to the actual next blocker
  - 47 tests passed
  - `h2_h1_readiness_elapsed_s=93.344068042`
- `julia --project=. test/nested/cartesian_ham_builder_h2_pqs_q5_gausslet_only_r4_readiness_runtests.jl`
  - passed
  - 47 tests passed
  - test summary time `1m33.7s`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - passed: `load ok`
- `git diff --check`
  - passed with no output
- `git diff --numstat -- src test bin tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl`
  - `12	0	src/pqs_source_box_route_driver_helpers.jl`
  - `7	1	src/pqs_source_box_route_driver_reporting.jl`
  - `1	1	test/driver_inputs/h2_pqs_q5_gausslet_only_r4.jl`
  - `9	3	test/nested/cartesian_ham_builder_h2_pqs_q5_gausslet_only_r4_readiness_runtests.jl`
  - `0	389	test/nested/cartesian_pair_block_one_body_record_dispatch_runtests.jl`
- `git status --short --branch`
  - `## main...origin/main`
  - modified:
    - `src/pqs_source_box_route_driver_helpers.jl`
    - `src/pqs_source_box_route_driver_reporting.jl`
    - `test/driver_inputs/h2_pqs_q5_gausslet_only_r4.jl`
    - `test/nested/cartesian_ham_builder_h2_pqs_q5_gausslet_only_r4_readiness_runtests.jl`
  - deleted:
    - `test/nested/cartesian_pair_block_one_body_record_dispatch_runtests.jl`

Timing note:

- The focused H2 readiness command exceeds 60 seconds.
- The first run spent about 57 seconds precompiling `GaussletBases`.
- Fresh-process route phases remained compilation dominated:
  - transform stage about 0.26s, about 99.9% compilation
  - pair terms about 1.0s, about 99.9% compilation
  - assembly about 32.6-33.1s, about 99.2-99.5% compilation
- Based on the printed phase timings, the observed cost is dominated by compilation, not repeated parent/final-basis/H1 runtime.

Scoped line-budget arithmetic:

- Scoped files: `src test bin tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl`
- Insertions: 29
- Deletions: 394
- Net: -365 lines

Deletion/shrinkage report:

- deleted:
  - `test/nested/cartesian_pair_block_one_body_record_dispatch_runtests.jl`
- simplified:
  - H2 readiness artifact/test now checks H1 materialization and scalar H1 diagnostics from the driver artifact.
  - The scalar H1 diagnostic report path reuses the existing H1 payload/final Hamiltonian data.
- quarantined:
  - none
- not deleted because:
  - The H2 readiness test is the active driver-owned artifact gate.
  - The existing diatomic complete core/shell H1 payload remains the route-owned H1 path used by this pass.
  - H1-J, private RHF, supplemented comparison, exports, and public solver behavior remain out of scope.
- exact remaining caller/blocker:
  - `:missing_diatomic_complete_core_shell_hamiltonian_handoff_payload`

-- repo-doer@macmini
