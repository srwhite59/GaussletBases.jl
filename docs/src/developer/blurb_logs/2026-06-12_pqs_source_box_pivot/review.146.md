Pass 146 review - accepted

Commit reviewed:
- `7dc7d654 Add diatomic source plan object`

Summary:
- Adds `_PQSDiatomicCompleteCoreShellSourcePlan`.
- Adds `_pqs_source_box_route_driver_diatomic_complete_core_shell_source_plan(...)`.
- Carries a private Be2/PQS source-plan object with:
  - `object_kind = :pqs_diatomic_complete_core_shell_source_plan`
  - real product/core support arrays
  - real left/right PQS shell support arrays
  - real block-diagonal shell coefficient matrix
  - retained/pre-final map summary
  - convention labels and nonpromotion flags
- Wires this object into `_PQSDiatomicCompleteCoreShellSourcePlanPayload` as
  `source_plan` on the probe-enabled path.

Review result:
- Accepted.
- The implementation follows the pass-145 audit decision: it does not fake or
  overload the one-center `:pqs_multilayer_shell_source_plan` object kind.
- The source-plan payload now blocks on
  `:missing_diatomic_complete_core_shell_final_basis_consumer`, which is the
  right next seam.

Validation checked by doer and repeated by manager:
- `julia --project=. test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl`
  - manager repeat passed:
    - `Be2 PQS Ham payload readiness fingerprint | Pass 151 Total 151 Time 12.6s`
    - `Be2 PQS probe-enabled Ham readiness fingerprint | Pass 190 Total 190 Time 38.0s`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - passed: `load ok`
- `git diff --check HEAD~1..HEAD`
  - passed

Boundary status:
- Good. No final-basis, H1, H1-J, Ham, support-one-body, support-density, RHF,
  WL, public API, export, artifact, hfdmrg, or CR2 behavior was added.
- The old one-center source-plan contract remains untouched.

Next step:
- Add a private diatomic final-basis payload that consumes
  `:pqs_diatomic_complete_core_shell_source_plan` and calls the lower
  complete-core/shell final-basis realization seam directly.
- Stop before H1/H1-J/Ham materialization.

-- repo-manager@macmini
