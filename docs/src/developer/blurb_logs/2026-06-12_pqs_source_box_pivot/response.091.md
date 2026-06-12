Pass 091 complete: added a driver-facing complete core/shell PQS H1/J diagnostic route payload slot.

Files changed:
- `src/pqs_source_box_route_driver_helpers.jl`
- `docs/src/developer/pqs_near_term_final_basis_realization_plan.md`

What was added:
- Internal payload helper:
  - `_pqs_source_box_route_driver_complete_core_shell_h1_j_diagnostic_payload(...)`
- Internal support helpers:
  - `_pqs_source_box_route_driver_complete_core_shell_h1_j_missing_inputs(...)`
  - `_pqs_source_box_route_driver_blocked_complete_core_shell_h1_j_payload(...)`
  - `_pqs_source_box_route_driver_complete_core_shell_h1_j_summary(...)`
  - `_pqs_source_box_route_driver_complete_core_shell_h1_j_report_fields(...)`

Whether a real driver-facing diagnostic route payload was added:
- Yes, but the active driver state is blocked until the complete core/shell route inputs are present.
- The helper can materialize by calling `pqs_multilayer_complete_core_shell_h1_j_payload(...)` when it is given:
  - `pqs_multilayer_shell_region_plan`
  - `pqs_multilayer_shell_source_plan`
  - complete core/shell final basis
  - complete core/shell H1 payload
  - axis weights
  - raw pair factor terms
  - Coulomb expansion
- In the current driver path, those inputs are not yet carried, so `cartesian_assembly(...)` records:
  - `status = :blocked_missing_complete_core_shell_h1_j_route_inputs`
  - `blocker = :missing_complete_core_shell_h1_j_route_inputs`
  - `driver_route_materialized = false`

Where it enters the driver spine:
- `cartesian_assembly(...)` now carries `complete_core_shell_h1_j_diagnostic_payload` and the matching compact status/blocker summary fields.
- `cartesian_report(...)` exposes compact report fields derived from that assembly payload.
- `cartesian_materialization(...)` is unchanged.

Compact report fields exposed:
- `complete_core_shell_h1_j_diagnostic_summary`
- `complete_core_shell_h1_j_diagnostic_status`
- `complete_core_shell_h1_j_diagnostic_blocker`
- `complete_core_shell_h1_j_final_dimension`
- `complete_core_shell_h1_j_h1_energy`
- `complete_core_shell_h1_j_self_coulomb`
- `complete_core_shell_h1_j_density_gauge`
- `complete_core_shell_h1_j_driver_route_materialized`

What remains private/diagnostic:
- The payload is still a diagnostic route seam, not RHF/SCF.
- It does not build Fock matrices, density iteration, GTO, exports, artifacts, fixture policy, q ladders, or side-13 acceptance.
- Rejected density conventions remain outside the driver path.
- The current active report state is blocked because the driver still needs to feed real shellification/lowering-backed complete core/shell route data into the payload.

Validation:
- Compact driver dry-run smoke passed:
  - constructed a PQS source-box dry-run report;
  - asserted `complete_core_shell_h1_j_diagnostic_status == :blocked_missing_complete_core_shell_h1_j_route_inputs`;
  - asserted `complete_core_shell_h1_j_driver_route_materialized == false`.
- Focused H1 gate passed:
  - `julia --project=. test/nested/pqs_direct_retained_final_h1_runtests.jl`
  - 46 tests passed.
- Load check passed:
  - `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git diff --check` passed.
- Trailing whitespace scan on edited source/docs found no matches.

Validation note:
- I started the broad `test/nested/pqs_source_box_route_driver_report_runtests.jl` file while the test assertions were temporarily in that file. It reached the route-driver standard-unit report testset and then spent several silent minutes in unrelated materialization checks; I interrupted it and removed those permanent assertions from the slow file. It is not counted as a completed validation gate for this pass.

Deletion/shrinkage report:
- No old code path was deleted. The driver did not have a prior complete core/shell H1/J route payload to replace.
- The docs were updated so they no longer say the driver spine lacks any payload slot; they now state that the slot exists but is blocked on missing complete core/shell route inputs.
- No permanent test was added. The compact dry-run smoke is the live validation for this pass; adding assertions to the broad report integration file would have increased carrying cost and tied the new seam to a slow non-focused gate.
- Remaining stale/duplicate surface to retire next: tmp/work H1/J probe glue can shrink only after the driver feeds this payload from real shellification/lowering-backed source-plan, final-basis, H1, and support-density inputs.

-- repo-doer@macmini
