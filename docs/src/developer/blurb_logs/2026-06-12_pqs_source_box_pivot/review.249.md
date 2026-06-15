# Pass 249 manager review

Decision: accepted.

Commit reviewed:

- pending commit: expose independent H2 PQS H1 artifact facts

Scope reviewed:

- `src/pqs_source_box_route_driver_helpers.jl`
- `src/pqs_source_box_route_driver_reporting.jl`
- `test/nested/cartesian_unit_stage_low_order_policy_runtests.jl`

Findings:

- No blocking findings.
- The existing physical-gausslet H1 path already works for the independent H2
  PQS final basis when requested. This pass exposes the H1 status facts in the
  artifact and preserves the stronger `:missing_physical_gausslet_h1_j_builder`
  endpoint blocker when H1-J is not requested.
- The reported H1 result is finite, symmetric, dimension `471`, and has lowest
  one-electron energy `-0.7946037173365885`.
- H1-J, RHF, supplements, CR2, export, and public API remain blocked/off.
- The deletion offset removes stale selected-terminal-sidecar count/status
  mirror assertions while preserving compact private-development smoke.

Validation accepted:

- Doer ran package load; it passed in about 58 seconds.
- Doer ran the focused independent H2 PQS H1 driver/artifact check; it passed
  in about 80 seconds.
- Doer ran `git diff --check`; it passed.
- Manager reviewed the focused diff and accepted doer validation without
  duplicating the slow route run.

Line budget:

- Scoped `src + test + bin`: `17` added / `26` deleted, net `-9`.

Remaining blocker / next:

- Independent H2 PQS now has source plan, final basis, and H1. The next seam is
  H1-J/density convention review; do not jump straight to RHF or supplements.

-- repo-manager@macmini
