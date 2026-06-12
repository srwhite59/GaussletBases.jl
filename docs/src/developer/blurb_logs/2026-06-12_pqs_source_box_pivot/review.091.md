Accepted with a validation caveat.

The implementation adds a real internal driver-facing H1/J diagnostic payload
slot rather than only a free-floating report flag. The active driver path still
blocks with `:missing_complete_core_shell_h1_j_route_inputs`, and
`driver_route_materialized` remains false, which is the right status until the
driver actually carries the shellification/lowering-backed region plan, source
plan, complete core/shell final basis, H1 payload, axis weights, raw pair
numerator terms, and Coulomb expansion.

The report fields are compact and appropriate for this stage: status, blocker,
final dimension, H1 energy, self-Coulomb, density gauge, and route
materialization status. `cartesian_materialization(...)` was left unchanged,
and the pass did not add RHF, SCF, GTO, fixture policy, exports, artifacts, or
WL/fixed-block authority.

The slow-test issue is real. The broad report/materialization integration file
is not a good pass-091 gate, and the assembly/report low-order policy tests are
also too heavy for this kind of narrow driver-slot change. During review, the
assembly/report policy tests ran silently for more than 90 seconds; the report
test also emitted failures in old low-order route-core assertions before being
stopped. Do not use those broad files as routine baton validation for this
seam.

Validation accepted for this pass:

- compact dry-run smoke of the PQS source-box report H1/J blocked status;
- `julia --project=. test/nested/pqs_direct_retained_final_h1_runtests.jl`;
- `julia --project=. -e 'using GaussletBases; println("load ok")'`;
- `git diff --check`.

Deletion/shrinkage review:

- No old code was made obsolete yet; the driver had no prior H1/J diagnostic
  payload slot.
- No permanent test was added, which is correct. Adding the assertions to the
  broad report integration file would have increased carrying cost.
- The docs were updated from "no driver payload slot exists" to "slot exists
  but is blocked on missing route inputs."
- Next stale surface to shrink remains the `tmp/work` H1/J probe glue, but only
  after the driver feeds this payload from real complete core/shell route data.

Next pass should not add more report fields. It should either feed the slot
from real shellification/lowering-backed complete core/shell route inputs, or
stop with the missing input boundary if that is not small.

-- repo-manager@macmini
