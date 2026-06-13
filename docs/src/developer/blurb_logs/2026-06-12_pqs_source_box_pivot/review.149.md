Pass 149 review - accepted and committed by manager

Implementation commit:
- `cd7d3fb0 Add diatomic H1 payload`

Response reviewed:
- `.agent_handoffs/response.149.md`

Summary:
- Adds `_PQSDiatomicCompleteCoreShellH1Payload`.
- Adds `_pqs_source_box_route_driver_diatomic_complete_core_shell_h1_payload(...)`.
- Adds private diatomic support-one-body helpers for:
  - support kinetic
  - support electron-nuclear by-center matrices
- Wires the payload into `cartesian_assembly(...)` as
  `diatomic_complete_core_shell_h1_payload`.

Review result:
- Accepted.
- The implementation keeps the old one-center support/H1 object-kind guards
  unchanged.
- The Be2/PQS route now has a private H1 diagnostic payload while still not
  claiming a full Hamiltonian payload.

Observed probe-enabled facts:
- H1 payload status:
  `:available_diatomic_complete_core_shell_h1_payload`
- Final dimension: `221`
- Lowest H1 energy: `-0.27746109235228694`
- Center count: `2`
- Readiness blocker advanced to:
  `:missing_diatomic_complete_core_shell_h1_j_consumer`

Validation checked by doer and repeated by manager:
- `julia --project=. test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl`
  - manager repeat passed:
    - `Be2 PQS Ham payload readiness fingerprint | Pass 182 Total 182 Time 12.9s`
    - `Be2 PQS probe-enabled Ham readiness fingerprint | Pass 259 Total 259 Time 39.4s`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - passed: `load ok`
- `git diff --check`
  - passed

Boundary status:
- Good. No H1-J, density-density, electron-electron Ham payload, full Ham
  payload, RHF, WL payload, public API, export, artifact, hfdmrg, or CR2
  behavior was added.
- This remains private diagnostic route behavior.

Manager note:
- The next blocker label names H1-J, but the medium-term target is a Hamiltonian
  constructor usable by downstream Be2 WL/PQS comparisons. H1-J is a useful
  diagnostic validator, not necessarily the next production-facing Hamiltonian
  seam. The next pass should audit this choice before coding.

-- repo-manager@macmini
