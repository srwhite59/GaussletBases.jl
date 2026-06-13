Pass 132 manager review

Accepted as the right readiness audit for the medium-term CR2 comparison goal.

Plain-language state:

- One-center PQS now has a private object-carrying Ham payload through
  `cartesian_assembly(...)`.
- Be2/diatomic PQS has route skeletons, inventories, and report metadata, but
  still does not have a route-owned final-basis/H1/Ham payload seam.
- WL/low-order Be2 has more materializer and Ham-adapter machinery, but it does
  not yet have a compact WL Ham payload boundary parallel to the PQS payload.
- CR2 cannot yet consume a comparable WL/PQS Be2 handoff without parsing broad
  report/materialization surfaces.

Decision:

- The immediate blocker is the missing Be2/diatomic PQS final-basis/H1/Ham
  payload seam.
- Next pass should be a focused Be2 PQS Ham readiness fingerprint: build the
  Be2 staged PQS route to `cartesian_assembly(...)` and record the exact
  private Ham payload status/blocker/missing inputs.
- Do not broaden to WL payload implementation yet.
- Do not implement export/public API or run CR2/hfdmrg.

Medium-term target:

```text
GaussletBases driver
-> comparable WL and PQS Be2 payloads
-> explicit Hamiltonian/final-basis conventions
-> CR2 agent can run an external comparison/probe
```

Deletion/shrinkage assessment:

- deleted: none.
- simplified: next work is narrowed to a PQS Be2 blocker fingerprint instead
  of broad driver usability.
- quarantined: RHF, public API, exports/artifacts, CR2 execution, hfdmrg, and
  physics endpoint claims remain out of scope.
- not deleted because: current report/materialization/export surfaces remain
  compatibility paths for existing driver tooling.
- exact remaining caller/blocker: Be2 PQS lacks a route-owned final-basis/H1/Ham
  payload seam; the next pass should make that blocker concrete.

-- repo-manager@macmini
