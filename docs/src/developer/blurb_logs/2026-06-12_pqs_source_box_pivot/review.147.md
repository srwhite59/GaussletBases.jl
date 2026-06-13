Pass 147 review - accepted

Commit reviewed:
- `8339dccd Add diatomic final basis payload`

Summary:
- Adds `_PQSDiatomicCompleteCoreShellFinalBasisPayload`.
- Adds `_pqs_source_box_route_driver_diatomic_complete_core_shell_final_basis_payload(...)`.
- Wires the payload into `cartesian_assembly(...)` as
  `diatomic_complete_core_shell_final_basis_payload`.
- Builds the Be2/PQS private final basis by calling the lower
  `CartesianFinalBasisRealization.pqs_complete_core_shell_final_basis(...)`
  helper directly, without pretending the source plan is the old one-center
  `:pqs_multilayer_shell_source_plan`.

Review result:
- Accepted.
- The probe-enabled Be2/PQS route now has:
  - private diatomic source plan available
  - private diatomic final basis available
  - final retained count `221`
  - next blocker `:missing_diatomic_complete_core_shell_h1_consumer`
- The default Be2/PQS path remains blocked upstream, as expected.

Validation checked by doer and repeated by manager:
- `julia --project=. test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl`
  - manager repeat passed:
    - `Be2 PQS Ham payload readiness fingerprint | Pass 165 Total 165 Time 13.0s`
    - `Be2 PQS probe-enabled Ham readiness fingerprint | Pass 221 Total 221 Time 38.7s`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - passed: `load ok`
- `git diff --check HEAD~1..HEAD`
  - passed

Boundary status:
- Good. No H1, H1-J, Ham, support-density, RHF, WL, public API, export,
  artifact, hfdmrg, or CR2 behavior was added.
- The only support operators materialized were overlap blocks needed by the
  final-basis realization helper.

Next step:
- Audit the diatomic H1 consumer seam before coding. H1 needs kinetic and
  electron-nuclear support operators, center/nuclear convention, and final-basis
  transfer. That should be checked against the existing one-center support
  one-body helpers before broadening any guards.

-- repo-manager@macmini
