Pass 134 manager review

Accepted.

Plain-language state:

- The one-center complete-core/shell path is terminal-shellification backed:
  parent axis bundle plus shellification/lowering plans produce the region plan,
  source plan, final basis, H1 payload, density inputs, H1/J diagnostic, and
  private Ham payload.
- The Be2/PQS route is already structured as a source-box-first diatomic route
  with `(:pqs_left, :product, :pqs_right)` units, source boxes, retained units,
  pair entries, center metadata, and route-axis counts.
- That Be2/PQS route structure is not the same object as the one-center
  collapsed complete-core/shell source plan, and a thin adapter into the
  one-center producer would blur the route semantics.

Decision:

- Add a compact private diatomic/PQS complete-core/shell Ham readiness payload.
- The payload should fingerprint what the Be2/PQS route already owns and what is
  missing before final-basis/H1/Ham construction.
- It should preserve the current blocked Ham payload behavior.
- It should not build final basis, H1, H1/J, RHF, WL payloads, exports,
  artifacts, or CR2/hfdmrg handoff data.

Deletion/shrinkage assessment:

- deleted: none.
- simplified: next work can replace scattered downstream missing-input checks
  in the Be2 fingerprint with one route-owned readiness object plus the
  preserved Ham blocker.
- quarantined: shell/support-row contraction, WL materializers/adapters, H1/J
  density interaction, RHF, exports, artifacts, CR2 execution, hfdmrg, and
  scalar report aliases remain out of scope.
- not deleted because: existing report tests still protect current route
  vocabulary and materialization metadata until the private payload seam
  matures.
- exact remaining caller/blocker: Be2/PQS lacks a route-owned diatomic
  complete-core/shell source-plan producer and therefore lacks final basis, H1,
  density-input, H1/J, and Ham materialization data for the private payload.

-- repo-manager@macmini
