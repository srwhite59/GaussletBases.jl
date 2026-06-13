Pass 153 review - accepted

Commit reviewed:
- `6cd55566 Add diatomic Hamiltonian handoff payload`

Summary:
- Adds `_PQSDiatomicCompleteCoreShellHamiltonianHandoffPayload`.
- Adds `_pqs_source_box_route_driver_diatomic_complete_core_shell_hamiltonian_handoff_payload(...)`.
- Wires it into `cartesian_assembly(...)` as
  `diatomic_complete_core_shell_hamiltonian_handoff_payload`.
- Carries references to the source plan, final basis, H1 payload, Ham-input
  payload, one-body Hamiltonian, pre-final density interaction, pair matrix,
  final-to-pre-final coefficients, weights, raw-pair data, and center metadata.

Review result:
- Accepted.
- This is the first private inspect-only Be2/PQS Hamiltonian handoff object.
- It correctly does not claim dense `Vee`, public export, HamV6 compatibility,
  CR2/HFDMRG readiness, RHF, WL comparison, or H1-J materialization.

Observed probe-enabled facts:
- Final dimension: `221`
- Pre-final dimension: `221`
- Density gauge: `:pre_final_localized_positive_weight`
- Raw-pair convention: `:raw_numerator`
- Support weight count: `275`
- Pre-final pair matrix shape: `(221, 221)`
- Final-to-pre-final coefficient shape: `(221, 221)`
- Nuclear charges: `(4.0, 4.0)`
- Nuclear coordinates: `((-2.0, 0.0, 0.0), (2.0, 0.0, 0.0))`
- Nuclear repulsion: `4.0`
- Electron count: `8`
- Spin sector: `:closed_shell_singlet`
- Readiness blocker:
  `:missing_diatomic_hamiltonian_consumer_contract`

Validation checked by doer and repeated by manager:
- `julia --project=. test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl`
  - manager repeat passed:
    - `Be2 PQS Ham payload readiness fingerprint | Pass 225 Total 225 Time 13.7s`
    - `Be2 PQS probe-enabled Ham readiness fingerprint | Pass 350 Total 350 Time 42.6s`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - passed: `load ok`
- `git diff --check HEAD~1..HEAD`
  - passed

Boundary status:
- Good. This is private inspect-only, not a downstream-consumer contract yet.

Important follow-up:
- Stop adding feature payloads to `src/pqs_source_box_route_driver_helpers.jl`
  for now. The conceptual route objects are sound, but the file has become a
  private route-driver subsystem and needs a behavior-preserving split before
  more H1-J/Ham/export/CR2 work.

-- repo-manager@macmini
