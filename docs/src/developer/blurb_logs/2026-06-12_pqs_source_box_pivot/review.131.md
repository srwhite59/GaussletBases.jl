Pass 131 manager review

Accepted and pushed.

Commit:

- `08c851ba Add private PQS complete core shell Ham payload`

The implementation adds the intended private Hamiltonian-constructor seam:

- `_pqs_source_box_route_driver_complete_core_shell_ham_payload(...)`
- `complete_core_shell_ham_payload` on
  `_PQSCompleteCoreShellDiagnosticRoutePayload`
- one focused route-smoke test:
  `test/nested/pqs_source_box_route_driver_complete_core_shell_ham_payload_runtests.jl`

The payload carries route-owned objects instead of report aliases: source
payload, final basis, H1 payload/final Hamiltonian, density inputs, H1/J
payload, and the pre-final density interaction. It keeps convention labels and
nonclaim flags explicit, including no public API, no export/artifact writing,
no RHF product surface, no serious-HF claim, and no IDA/MWG promotion.

Validation reported by doer:

- `julia --project=. test/nested/pqs_source_box_route_driver_complete_core_shell_ham_payload_runtests.jl`
  passed: 21/21, reported test time `1m52.3s`, including about 56.5s package
  precompilation.
- `julia --project=. -e 'using GaussletBases; println("load ok")'` passed.
- `git diff --check` passed.

Manager checks:

- `git show --check HEAD` passed.
- Commit touched only the route helper and the focused test.
- Commit was pushed to `origin/main`.

Decision:

- The driver now has a private PQS Ham payload seam for the compact one-center
  complete core/shell route-smoke.
- The medium-term goal is now CR2 trying a Be2 WL/PQS comparison against usable
  driver payloads.
- Next pass should be no-edit readiness audit for that goal: what WL and PQS
  payloads exist today, what is one-center-only, what is Be2/diatomic-ready,
  and what exact missing seam blocks a CR2 handoff.

Deletion/shrinkage assessment:

- deleted: none.
- simplified: downstream Ham handoff now has one private object boundary rather
  than H1/J report aliases plus RHF validator fields.
- quarantined: RHF remains validator-only; report/export/artifact paths remain
  outside this private payload; pre-final density interaction remains labeled
  and not promoted to IDA/MWG or full Coulomb tensor semantics.
- not deleted because: existing H1/J report aliases and materialization/export
  fields are compatibility surfaces.
- exact remaining caller/blocker: the payload is available through
  `assembly.complete_core_shell_diagnostic_route_payload.complete_core_shell_ham_payload`;
  no CR2/HFDMRG/export caller or Be2 WL/PQS comparison path exists yet.

-- repo-manager@macmini
