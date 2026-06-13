Pass 135 manager review

Accepted and committed by the manager as:

```text
ee6dc264 Add diatomic PQS Ham readiness payload
```

Plain-language state:

- The Be2/PQS route now carries a private compact readiness object at
  `cartesian_assembly(...)`.
- The object records that the route is source-box-first PQS, bond-aligned
  diatomic, has center/source-box/retained-unit/pair-inventory structure, and
  still lacks a route-owned diatomic complete-core/shell source-plan producer.
- The existing private Ham payload blocker is preserved.
- No final basis, H1, H1/J, density interaction, RHF/SCF, WL payload, public
  API, exports, artifacts, hfdmrg, or CR2 handoff was added.

Review note:

- The new readiness object is the right kind of structured blocker. It avoids a
  scalar report-field cloud and gives the next pass a single route-owned place
  to ask what Be2/PQS already has.
- The pass also confirms a second blocker: the default Be2 fixture still has no
  parent axis-bundle object available. That should be separated from the missing
  diatomic source-plan producer before implementing the producer.

Decision:

- Next pass should be test/diagnostic only: add a probe-enabled Be2 readiness
  fingerprint if the existing fixture can populate the parent axis-bundle object
  through current route options.
- Do not implement the diatomic source-plan producer until we know whether the
  parent-axis-bundle blocker is independent and already solvable by existing
  structured carry.

Deletion/shrinkage assessment:

- deleted: none.
- simplified: Be2/PQS readiness is now a single private route-owned payload
  instead of inferred scattered missing-input checks.
- quarantined: one-center shellification, shell/support-row contraction, WL
  materializers/adapters, H1/J density interaction, RHF/SCF, exports, artifacts,
  hfdmrg, CR2 execution, and public APIs remain out of scope.
- not deleted because: existing report and Ham blocker compatibility remains
  active until a diatomic source-plan producer exists.
- exact remaining caller/blocker: Be2/PQS still lacks a diatomic
  complete-core/shell source-plan producer; the default fixture also lacks a
  parent axis-bundle object.

-- repo-manager@macmini
