Pass 133 manager review

Accepted.

Plain-language state:

- The Be2 PQS route now has a focused tracked blocker fingerprint at the Ham
  payload seam.
- The staged Be2 route reaches `cartesian_assembly(...)`, and the private
  diagnostic route payload and Ham payload slots exist.
- The route does not yet own the complete-core/shell source plan, final basis,
  H1 payload, density inputs, H1/J diagnostic, or pre-final density interaction
  needed for the private Ham payload to materialize.
- The new fingerprint is standalone and not part of the default nested runner,
  which is the right test pressure for this stage.

Decision:

- Do not broaden to WL or CR2 yet.
- Do not add export/public API, RHF/SCF work, hfdmrg calls, or fixture
  promotion.
- Next pass should be a no-edit seam audit: identify what part of the one-center
  complete-core/shell source-plan/final-basis/H1 path is genuinely one-center,
  what can be generalized to Be2/diatomic PQS, and the smallest implementation
  seam that would move the Be2 Ham payload toward materialization.

Deletion/shrinkage assessment:

- deleted: none.
- simplified: the Be2 readiness gap is now a compact Ham blocker fingerprint,
  not an inferred broad-report state.
- quarantined: WL payload implementation, CR2 execution, hfdmrg, public/export
  surfaces, RHF/SCF, and physics endpoint claims remain outside the active pass.
- not deleted because: the existing report/materialization tests still protect
  compatibility paths while the new private payload seam is being built.
- exact remaining caller/blocker: Be2 PQS lacks the route-owned complete
  core/shell source plan, final basis, H1 payload/final H1 matrix, density
  inputs, H1/J diagnostic payload, and pre-final density interaction needed by
  the private Ham payload.

-- repo-manager@macmini
