Pass 090 complete: audited the PQS H1/J driver-spine integration point and recorded the design boundary.

Files changed:
- `docs/src/developer/pqs_near_term_final_basis_realization_plan.md`

Where H1/J should enter the driver spine:
- The compact H1/J seam belongs after `cartesian_assembly` has a real complete core/shell diagnostic route payload.
- That payload should carry the shellification/lowering-backed region plan, the PQS multilayer source plan, the complete core/shell final basis, the H1 payload, and the density inputs required by `pqs_multilayer_complete_core_shell_h1_j_payload(...)`.
- `cartesian_report` should expose only a compact diagnostic summary: status, blocker, H1 energy, self-Coulomb value, density gauge, and nonclaim flags.
- `cartesian_materialization` should remain a durable-summary/sanitization boundary unless it creates a real driver-owned route object.

Implementation decision:
- I did not add an internal optional assembly/report payload in this pass.
- The insertion point is clear, but the required route object is not currently present in the driver spine, so a placeholder report hook would add metadata without making the driver consume the H1/J seam.
- `driver_route_materialized` should remain false until a real driver-owned route object exists.

Precise missing route object:
- A small driver-facing complete core/shell diagnostic route payload, conceptually:
  - shellification/lowering-backed region plan;
  - `pqs_multilayer_shell_source_plan`;
  - `pqs_multilayer_complete_core_shell_final_basis`;
  - `pqs_multilayer_complete_core_shell_h1_payload`;
  - support density inputs from `pqs_multilayer_support_weights(...)` and `pqs_multilayer_support_pair_raw_numerator_matrix(...)`;
  - compact H1/J payload from `pqs_multilayer_complete_core_shell_h1_j_payload(...)`.

Audit findings:
- Existing `cartesian_assembly` carries staged skeleton/low-order pair readiness and does not yet carry the complete core/shell final-basis/H1/J route payload.
- Existing `cartesian_report` mirrors the low-order assembly/report summary and can carry compact status once assembly owns the route payload.
- Existing `cartesian_materialization` has a transient route materializer payload, but that payload only carries parent basis and parent axis-bundle handoff objects. It does not carry source plans, final basis data, H1 payloads, support weights, or raw support pair numerators.
- The current H1/J helper is still private/diagnostic. It is route-owned code, but not yet driver-owned route flow.

What remains private/diagnostic:
- `pqs_multilayer_complete_core_shell_h1_j_payload(...)` remains a compact diagnostic seam.
- `tmp/work` probes remain the exploratory consumers.
- Rejected signed-final-weight and raw-no-division comparisons remain probe-local.

Next design boundary:
- The next boundary is the small driver-facing complete core/shell diagnostic route object.
- RHF remains a later design boundary. The driver should not grow SCF/Fock/density iteration until H1/J is actually consumed through the driver route object.

Validation:
- `git diff --check` passed.

Deletion/shrinkage report:
- No source or tests were deleted. This pass was an audit/docs pass because no existing driver object was made obsolete by the H1/J helper.
- The docs were simplified away from the stale claim that lowest H1 orbital coefficients remain a route-owned gap; pass 089 moved that into the H1/J helper.
- No new test was added. A report-stage placeholder test would preserve metadata vocabulary instead of protecting a live driver contract.
- Remaining stale/duplicate surface to retire next: tmp/work H1/J probe glue should shrink only after the driver-facing complete core/shell diagnostic route payload exists and consumes the compact H1/J helper.

-- repo-doer@macmini
