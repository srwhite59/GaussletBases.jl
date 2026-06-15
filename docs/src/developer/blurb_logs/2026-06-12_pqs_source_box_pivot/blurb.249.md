# Pass 249 blurb - independent H2 PQS H1 one-body seam

Role: repo-doer.

Read before starting:

- `AGENTS.md`
- `BlurbStyle.md`
- `docs/src/developer/pqs_manager_running_log.md`
- `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.248.md`
- `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/review.248.md`

## Current state

The independent H2 PQS route now has:

```text
source_plan_status = :available_pqs_diatomic_physical_gausslet_core_shell_source_plan
final_basis_status = :available_pqs_physical_gausslet_final_basis
final_dimension = 471
retained_counts = (275, 98, 98)
final_overlap_identity_error ≈ 1.3e-13
fake_pqs_enabled = false
source_backed_fixed_source_oracle_used = false
retained_transform_authority = :pqs_source_box_construction
```

The next seam is one-body H1 materialization for this independent final basis.
There is already a physical-gausslet H1 path in:

- `src/pqs_source_box_diatomic_complete_core_shell.jl`

Relevant helpers include:

- `_pqs_source_box_route_driver_physical_gausslet_support_kinetic(...)`
- `_pqs_source_box_route_driver_physical_gausslet_support_electron_nuclear(...)`
- `_pqs_source_box_route_driver_physical_gausslet_final_one_body(...)`
- `_pqs_source_box_route_driver_physical_gausslet_h1_hamiltonian(...)`
- `_pqs_source_box_route_driver_physical_gausslet_h1_solve(...)`
- `_pqs_source_box_route_driver_diatomic_physical_gausslet_h1_payload(...)`

Use this path. Do not create a parallel H1 implementation.

## Task

Run and, if necessary, narrowly repair the independent H2 PQS H1 one-body seam
with `run_h1=true`.

The desired result is:

```text
h1_status = :materialized_pqs_physical_gausslet_h1_solve
h1_materialized = true
final_dimension = 471
h1_lowest_energy finite and negative
h1_hamiltonian_matrix_finite = true
h1_hamiltonian_symmetry_error small
support_kinetic_status materialized
support_electron_nuclear_status materialized
final_kinetic_status materialized
final_electron_nuclear_status materialized
```

If the H1 path is already working, do not add source code just to add source
code. Report the result and spend any needed line budget on stale staged-test
cleanup only.

If it is blocked, fix only the narrow H1 seam. Good blockers to resolve in this
pass are:

- missing route/report artifact fields for existing H1 payload facts;
- stale assumptions about the independent final-basis dimension;
- missing propagation of H1 support/final status fields.

If the H1 matrix is nonfinite, nonsymmetric, nonnegative at the lowest
eigenvalue, or depends on fake/source-backed data, keep it blocked and report
the exact blocker. Do not force success.

## Strict exclusions

Do not add or enable:

- H1-J or density-density interaction;
- RHF/private RHF;
- MWG/GTO supplement provider blocks;
- CR2/export/public API readiness;
- fake-PQS route changes except guard preservation;
- comparison against fake-PQS/WL 463 as a physics claim.

The endpoint should remain not solver/export/public-ready. If H1 materializes,
the next blocker should be H1-J/density convention, not RHF or export:

```text
physics/endpoint_ready = false
physics/endpoint_blocker = :missing_physical_gausslet_h1_j_builder
```

or an exact H1 blocker if H1 remains unavailable.

## Line budget

Keep scoped `src + test + bin` net-negative.

If the H1 path works without source additions, use this pass as a small
validation-plus-cleanup pass and delete stale terminal-shellification assertion
mirrors in the staged tests. If source additions are needed, offset them with
the same cleanup pressure.

Remaining stale deletion candidates include exact terminal-shellification
field mirrors and sidecar-count mirrors in:

- `test/nested/cartesian_unit_stage_low_order_policy_runtests.jl`
- `test/nested/cartesian_transform_stage_low_order_policy_runtests.jl`
- `test/nested/cartesian_assembly_stage_low_order_policy_runtests.jl`
- `test/nested/cartesian_report_stage_low_order_policy_runtests.jl`

Do not delete compact active smoke checks for route selection, scaffold
availability, deferral/materialization status, or final-basis/H1 behavior.

If a net-negative pass is impossible without broadening the task or weakening an
active contract, write `.agent_handoffs/ATTENTION.md` and stop.

## Validation

Run only focused validation:

- package load;
- focused independent H2 PQS driver/artifact check with `run_h1=true` and
  `run_h1_j=false`, `run_private_rhf=false`;
- `git diff --check`.

The focused driver check may exceed 60 seconds; report elapsed time and why it
is necessary. Do not run:

- H1-J/RHF/supplement/CR2 tests;
- full nested suite;
- broad staged low-order policy gates;
- `test/nested/cartesian_pair_stage_low_order_policy_runtests.jl`.

## Report

In `response.249.md`, report:

- H1 status, blocker, lowest energy, final dimension, finiteness, and symmetry
  error;
- support/final kinetic and electron-nuclear statuses;
- fake/source-backed guard fields;
- evidence that H1-J, RHF, supplements, CR2, export, and public API remain
  blocked;
- line budget and deletion offset;
- validation commands and elapsed times.

Do not commit. Leave the worktree ready for manager review.

-- repo-manager@macmini
