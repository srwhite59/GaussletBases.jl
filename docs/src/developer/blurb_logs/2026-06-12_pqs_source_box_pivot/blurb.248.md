# Pass 248 blurb - independent H2 PQS final-basis seam

Role: repo-doer.

Read before starting:

- `AGENTS.md`
- `BlurbStyle.md`
- `docs/src/developer/pqs_manager_running_log.md`
- `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.247.md`
- `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/review.247.md`

## Current state

Pass 247 assembled an available route-owned source-plan payload for the
independent H2 PQS target:

```text
source_plan_status = :available_pqs_diatomic_physical_gausslet_core_shell_source_plan
support_counts = (275, 578, 362)
retained_counts = (275, 98, 98)
final_dimension = 471
fake_pqs_enabled = false
source_backed_fixed_source_oracle_used = false
retained_transform_authority = :pqs_source_box_construction
```

The default input still has `run_final_basis = false`, and that should remain
the default in this pass:

- `test/driver_inputs/h2_pqs_q5_independent_source_box_r4.jl`

The next blocker is the final-basis seam. The existing helper:

- `_pqs_source_box_route_driver_physical_gausslet_final_basis(...)`

still contains old physical-gausslet assumptions, including:

```text
retained_counts == (251, 98, 114)
```

and its success metadata still has stale source-backed wording. That must not be
allowed to turn the new source plan into another fake/source-backed path.

## Task

Make the independent H2 PQS final-basis seam either:

1. materialize a final basis from the new route-owned `(275, 98, 98)` source
   plan under an explicit independent-PQS authority label; or
2. remain blocked with a narrower, correct blocker explaining the exact missing
   final-basis object.

Prefer option 1 if it is a small local correction.

If final basis materializes, it must report:

```text
final_basis_status = :available_pqs_physical_gausslet_final_basis
final_dimension = 471
retained_counts = (275, 98, 98)
source_backed_fixed_source_oracle_used = false
fake_pqs_enabled = false
retained_transform_authority = :pqs_source_box_construction
final_basis_materialized = true
h1_materialized = false
h1_j_materialized = false
rhf_materialized = false
exports_materialized = false
public_api = false
```

Also record overlap diagnostics:

```text
pre_final_overlap_identity_error
final_overlap_identity_error
rank / full-rank status
```

The final overlap should remain diagnostic/debug data only. Do not introduce a
downstream generalized-overlap workflow.

## Important implementation guidance

- Do not change the default driver input to request final basis. Use a focused
  validation override such as `run_final_basis=true`.
- Do not reuse the fake-PQS/WL source-backed candidate.
- Do not compare this to the fake-PQS/WL 463 endpoint as a physics result.
- Do not start H1, H1-J, RHF, supplement, CR2, export, or public API work.
- If the existing dense support-overlap/final-cleanup code is used, label it as
  first final-basis diagnostic materialization, not a production performance
  claim.
- If the final-basis construction exposes a rank/conditioning problem, keep the
  route blocked and report the numerical diagnostic. Do not force success.

## Line budget

Keep scoped `src + test + bin` net-negative.

If source/test assertions are needed, offset them by deleting stale exact
terminal-shellification assertion pressure in the staged tests. Remaining
candidates include exact lowering-contract and selected-contract mirrors in:

- `test/nested/cartesian_unit_stage_low_order_policy_runtests.jl`
- `test/nested/cartesian_transform_stage_low_order_policy_runtests.jl`
- `test/nested/cartesian_assembly_stage_low_order_policy_runtests.jl`
- `test/nested/cartesian_report_stage_low_order_policy_runtests.jl`

Do not delete active terminal route/scaffold/materialization smoke. Do not run
the broad staged tests as validation gates.

If a net-negative pass is impossible without broadening the task or weakening an
active contract, write `.agent_handoffs/ATTENTION.md` and stop.

## Validation

Run only focused validation:

- package load;
- focused independent H2 PQS driver/artifact check with
  `run_final_basis=true`, asserting final-basis status or the exact blocker;
- `git diff --check`.

The focused driver check may exceed 60 seconds; report elapsed time and why it
is necessary. Do not run:

- H1/HF/RHF/H1-J/supplement/CR2 tests;
- full nested suite;
- broad staged low-order policy gates;
- `test/nested/cartesian_pair_stage_low_order_policy_runtests.jl`.

## Report

In `response.248.md`, report:

- whether final basis materialized or the exact blocker;
- final dimension, retained counts, and overlap/rank diagnostics if available;
- fake/source-backed guard fields;
- evidence that H1, H1-J, RHF, supplements, CR2, export, and public API remain
  blocked;
- line budget and deletion offset;
- validation commands and elapsed times.

Do not commit. Leave the worktree ready for manager review.

-- repo-manager@macmini
