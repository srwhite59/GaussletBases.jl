# Pass 250 blurb - independent H2 PQS H1-J density diagnostic

Role: repo-doer.

Read before starting:

- `AGENTS.md`
- `BlurbStyle.md`
- `docs/src/developer/pqs_manager_running_log.md`
- `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.249.md`
- `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/review.249.md`

## Current state

Independent H2 PQS now has:

```text
source_plan_status = :available_pqs_diatomic_physical_gausslet_core_shell_source_plan
final_basis_status = :available_pqs_physical_gausslet_final_basis
final_dimension = 471
h1_status = :materialized_pqs_physical_gausslet_h1_solve
h1_lowest_energy = -0.7946037173365885
h1_hamiltonian_matrix_finite = true
h1_hamiltonian_symmetry_error ≈ 1.1e-14
fake_pqs_enabled = false
source_backed_fixed_source_oracle_used = false
retained_transform_authority = :pqs_source_box_construction
```

The next seam is H1-J / density interaction diagnostics. This is deliberately
not RHF and not solver readiness. Earlier PQS density attempts exposed bad
weight/gauge behavior, so this pass must be conservative about conventions.

There is already an H1-J path in:

- `src/pqs_source_box_diatomic_complete_core_shell.jl`

Relevant helpers include:

- `_pqs_source_box_route_driver_physical_gausslet_density_provenance(...)`
- `_pqs_source_box_route_driver_physical_gausslet_support_weights(...)`
- `_pqs_source_box_route_driver_physical_gausslet_support_pair_raw_numerator_matrix(...)`
- `_pqs_source_box_route_driver_physical_gausslet_pre_final_density_interaction(...)`
- `_pqs_source_box_route_driver_physical_gausslet_h1_j_diagnostic(...)`
- `_pqs_source_box_route_driver_diatomic_physical_gausslet_h1_j_payload(...)`

Use this path. Do not create a parallel density/H1-J implementation.

## Task

Run and, if necessary, narrowly repair the independent H2 PQS H1-J diagnostic
seam with:

```text
run_h1 = true
run_h1_j = true
run_private_rhf = false
```

If the existing path works, expose/report the resulting density/H1-J facts and
do not add source code just to add source code.

Desired coherent diagnostic result:

```text
h1_j_status = :materialized_pqs_physical_gausslet_h1_j_payload
h1_j_materialized = true
density_interaction_status = :materialized_pqs_physical_gausslet_pre_final_density_interaction
density_interaction_materialized = true
density_gauge = :pre_final_localized_positive_weight
raw_pair_factor_convention = :raw_numerator
support_weight_count finite / expected from support plan
support_weights_all_positive = true
support_raw_pair_shape finite
pre_final_pair_matrix_shape finite
pre_final_pair_matrix_finite = true
pre_final_pair_matrix_symmetry_error small
h1_j_self_coulomb finite
```

If the path is blocked, keep it blocked with the exact convention/object
blocker. Good blockers to preserve or improve are:

- missing axis weights;
- missing raw pair factor terms;
- nonpositive support weights;
- blocked pre-final density interaction;
- nonfinite or nonsymmetric pre-final pair matrix;
- final-orbital-to-pre-final shape mismatch.

Do not force success if the density convention looks wrong. A finite H1-J value
is a diagnostic, not a solver contract.

## Strict exclusions

Do not add or enable:

- RHF/private RHF;
- final density iteration or solver contract;
- MWG/GTO supplement provider blocks;
- CR2/export/public API readiness;
- fake-PQS route changes except guard preservation;
- comparison against fake-PQS/WL 463 as a physics claim.

If H1-J materializes, endpoint readiness must still be false. The next blocker
should be solver/RHF convention review, not public readiness:

```text
physics/endpoint_ready = false
physics/endpoint_blocker = :missing_physical_gausslet_rhf_or_solver_contract
```

or the exact H1-J blocker if H1-J remains unavailable.

## Line budget

Keep scoped `src + test + bin` net-negative.

If the H1-J path works without source additions, use this pass as a
validation-plus-cleanup pass and delete stale staged-test mirror assertions. If
source/report additions are needed, offset them with the same cleanup pressure.

Remaining deletion candidates include exact terminal-shellification field
mirrors and sidecar-status/count mirrors in:

- `test/nested/cartesian_unit_stage_low_order_policy_runtests.jl`
- `test/nested/cartesian_transform_stage_low_order_policy_runtests.jl`
- `test/nested/cartesian_assembly_stage_low_order_policy_runtests.jl`
- `test/nested/cartesian_report_stage_low_order_policy_runtests.jl`

Do not delete compact active smoke checks for route selection, scaffold
availability, deferral/materialization status, final-basis, H1, or H1-J
behavior.

If a net-negative pass is impossible without broadening the task or weakening an
active contract, write `.agent_handoffs/ATTENTION.md` and stop.

## Validation

Run only focused validation:

- package load;
- focused independent H2 PQS driver/artifact check with `run_h1_j=true` and
  `run_private_rhf=false`;
- `git diff --check`.

The focused driver check may exceed 60 seconds; report elapsed time and why it
is necessary. Do not run:

- RHF/supplement/CR2 tests;
- full nested suite;
- broad staged low-order policy gates;
- `test/nested/cartesian_pair_stage_low_order_policy_runtests.jl`.

## Report

In `response.250.md`, report:

- H1-J status/blocker and whether the density interaction materialized;
- density gauge, raw pair convention, support-weight facts, pair-matrix shape,
  finiteness, symmetry error, and H1-J self-Coulomb if available;
- H1 status remains materialized and final dimension remains 471;
- fake/source-backed guard fields;
- evidence that RHF, supplements, CR2, export, and public API remain blocked;
- line budget and deletion offset;
- validation commands and elapsed times.

Do not commit. Leave the worktree ready for manager review.

-- repo-manager@macmini
