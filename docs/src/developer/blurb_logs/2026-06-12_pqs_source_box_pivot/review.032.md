Review 032

Accepted.

The doc update fixes the CRPS ownership inconsistency and records the PQS
final-basis H1 status at the right level: an oracle-backed validated seam, not
a fully production-owned PQS route. It does not change source code, tests, or
route behavior.

Validation:

- `git diff --check` passed for the edited docs.
- ASCII scan on edited docs found no matches.
- No Julia tests were needed for this documentation-only pass.

Deletion/shrinkage:

- Stale framework wording that implied CRPS owned no retained-rule-adjacent
  source-mode selector was corrected.
- The new H1 status note earns its cost by preventing two wrong readings:
  "PQS still cannot reach H1" and "PQS H1 is production route-owned."
- Remaining stale surfaces: final-basis realization/operator-transfer helpers
  still live inside CPBM, and the raw-source dense block plus selector path is
  still the retained-boundary overlap/kinetic implementation.

Next:

- Audit the mechanical extraction boundary for `CartesianFinalBasisRealization`
  before moving code.

-- repo-manager@macmini
