# Pass 247 blurb - assemble independent H2 PQS source-plan payload

Role: repo-doer.

Read before starting:

- `AGENTS.md`
- `BlurbStyle.md`
- `docs/src/developer/pqs_manager_running_log.md`
- `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.241.md`
- `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/review.246.md`

## Current state

The fake-PQS H2 463 route is quarantined and must remain only a
source-backed WL/QW reproduction. Do not use it as evidence for independent
PQS.

The independent H2 PQS route currently has:

- route kind: `:bond_aligned_diatomic_independent_pqs_source_box_core_shell`
- input: `test/driver_inputs/h2_pqs_q5_independent_source_box_r4.jl`
- target support counts: `(275, 578, 362)`
- target retained counts from independent PQS pieces: `(275, 98, 98)`
- expected source-plan final dimension: `471`
- descriptor status:
  `:available_independent_pqs_physical_source_plan_descriptor`
- shared-shell realization status:
  `:available_independent_pqs_shared_shell_realization_payload`
- shared-shell counts: `(98, 98)`
- shared-shell identity errors around `1e-14` and `5e-14`
- current source-plan blocker:
  `:missing_independent_pqs_complete_core_shell_source_plan_assembly`

The relevant implementation seam is in:

- `src/pqs_source_box_diatomic_complete_core_shell.jl`

Look first near these helpers:

- `_pqs_source_box_route_driver_independent_h2_physical_source_plan_descriptor(...)`
- `_pqs_source_box_route_driver_independent_h2_shared_shell_realization_payload(...)`
- `_pqs_source_box_route_driver_diatomic_physical_gausslet_source_plan_payload(...)`
- `_pqs_source_box_route_driver_physical_gausslet_final_basis(...)`

The final-basis helper already expects a complete source-plan object with fields
such as `atom_contact_core_support_indices`, `shared_shell_support_indices`,
`core_coefficient_matrix`, `shared_shell_coefficient_matrices`,
`retained_ranges`, and `axis_bundles`. This pass should assemble that source-plan
payload; it should not invoke or validate final-basis construction.

## Task

Assemble the complete independent H2 PQS core/shell source-plan payload from the
existing independent descriptor, atom-contact support information, and the
shared-shell realization payload.

The source-plan payload should become available if the assembly is complete:

```text
source_plan_status = :available_pqs_diatomic_physical_gausslet_core_shell_source_plan
source_plan_blocker = nothing
```

If you find that this cannot be done safely, keep the source plan blocked but
replace the broad blocker with the exact missing object.

The source plan must record, directly or through a compact summary:

```text
support_order = (:atom_contact_core, :shared_shell_1, :shared_shell_2)
retained_order = (:atom_contact_core, :shared_shell_1, :shared_shell_2)
support_counts = (275, 578, 362)
retained_counts = (275, 98, 98)
final_dimension = 471
source_backed_fixed_source_oracle_used = false
fake_pqs_enabled = false
retained_transform_authority = :pqs_source_box_construction
source_plan_authority_status = :independent_pqs_route_owned_source_plan
```

For `:atom_contact_core`, use the route-owned atom-contact support/core data
already identified by the independent target. Treat it as the direct retained
core in this pass, but avoid materializing a large dense identity if a compact
identity-like representation is enough. If a dense matrix is unavoidable, label
it as an implementation detail of the current source-plan assembly, not a
production dense-parent/operator path.

For the two shared shells, consume the existing shared-shell realization payload
from pass 241. Preserve the shared-shell count and identity-error diagnostics.

## Strict exclusions

Do not add or enable:

- final-basis materialization;
- H1, H1-J, density interaction, or RHF;
- supplement provider blocks or MWG/GTO supplement materialization;
- CR2/export/public API readiness;
- fake-PQS route changes, except if needed to preserve guard assertions;
- source-backed WL/QW fixed-source coefficient or retained-transform use;
- broad report-key clouds or duplicate scalar aliases.

The endpoint must remain not physics-ready:

```text
physics/endpoint_ready = false
physics/endpoint_blocker = :missing_physical_gausslet_final_basis_request
```

or another exact final-basis-not-requested/missing-final-basis blocker. The
important point is that H1/H1-J/RHF/supplement/export remain blocked.

## Line budget

This pass must be net-negative in scoped `src + test + bin`.

If assembling the source-plan object adds source lines, offset it by deleting a
small stale test block. Good deletion-offset candidates are the old exact
terminal-shellification field assertions in:

- `test/nested/cartesian_unit_stage_low_order_policy_runtests.jl`
- `test/nested/cartesian_transform_stage_low_order_policy_runtests.jl`
- `test/nested/cartesian_assembly_stage_low_order_policy_runtests.jl`
- `test/nested/cartesian_report_stage_low_order_policy_runtests.jl`

These have repeatedly shown up as stale exact-field pressure in broad staged
tests. Do not run those broad staged tests as validation gates for this pass.

Target net result: at least `-20` lines in scoped `src + test + bin`.

If you cannot keep the pass net-negative without broadening the task or deleting
an active contract, write `.agent_handoffs/ATTENTION.md` and stop.

## Validation

Run only focused validation:

- package load;
- the focused independent H2 PQS driver/readiness artifact check used in pass
  241, updated to assert source-plan availability and the new source-plan
  fields;
- `git diff --check`.

Do not run:

- stale broad staged unit/transform/assembly/report gates;
- `test/nested/cartesian_pair_stage_low_order_policy_runtests.jl`;
- full nested suite;
- CR2 tests;
- H1/HF/RHF/supplement tests.

If the focused driver check takes more than 60 seconds, report elapsed time and
why it was necessary. That is acceptable here because it is the only focused
artifact route check for the source-plan seam.

## Report

In `response.247.md`, report:

- source-plan status, blocker, support counts, retained counts, and final
  dimension;
- evidence that fake-PQS/source-backed fixed-source paths remain unused;
- evidence that final basis, H1, H1-J, RHF, supplements, CR2, export, and public
  API remain blocked;
- line budget for scoped `src + test + bin`;
- any deletion offset used and why it was stale;
- validation commands and elapsed times.

Do not commit. Leave the worktree ready for manager review.

-- repo-manager@macmini
