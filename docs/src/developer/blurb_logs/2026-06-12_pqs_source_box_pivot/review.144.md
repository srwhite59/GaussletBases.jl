Pass 144 review - accepted

Commit reviewed:
- `7a79d25a Add diatomic source realization payload`

Summary:
- Adds `_PQSDiatomicCompleteCoreShellSourceRealizationPayload`.
- Adds `_pqs_source_box_route_driver_diatomic_complete_core_shell_source_realization_payload(...)`.
- Wires the payload into `cartesian_assembly(...)` as
  `diatomic_complete_core_shell_source_realization_payload`.
- Extends the diatomic source-plan payload to carry and summarize the
  source-realization payload status while still keeping `source_plan === nothing`.

Review result:
- The implementation matches the pass-144 boundary.
- Default Be2/PQS still blocks on missing parent axis-bundle/raw-box route input.
- Probe-enabled Be2/PQS now reaches an available private diatomic source
  realization with:
  - product core/body unit
  - left/right PQS shell/source units
  - retained order `(:pqs_left, :pqs_right, :product)`
  - support order `(:product, :pqs_left, :pqs_right)`
  - explicit retained/support-order permutation requirement
  - precleanup retained dimension `221`
  - shell coefficient shape `(250, 196)`
  - object kind claim `:not_pqs_multilayer_shell_source_plan`
- The source-plan payload now blocks specifically on
  `:missing_pqs_multilayer_shell_source_plan_adapter_contract` for the
  probe-enabled Be2/PQS route.

Validation checked by doer and repeated by manager:
- `julia --project=. test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl`
  - manager repeat passed:
    - `Be2 PQS Ham payload readiness fingerprint | Pass 151 Total 151 Time 12.6s`
    - `Be2 PQS probe-enabled Ham readiness fingerprint | Pass 154 Total 154 Time 36.6s`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - passed: `load ok`
- `git diff --check HEAD~1..HEAD`
  - passed

Notes:
- The tracked `response.144.md` was present, but the ignored live
  `.agent_handoffs/response.144.md` copy was missing at review time. The
  review used the tracked response and direct commit inspection.
- The focused Be2 test took under 60 seconds once package precompilation was
  not repeated.

Boundary status:
- Good. No source-plan/final-basis/H1/H1-J/Ham/RHF/WL/public/export/artifact
  behavior was added.
- The remaining blocker is now the semantic adapter/new-object contract between
  the diatomic source realization and the complete-core/shell final-basis/H1
  consumers.

Next step:
- Pass 145 should be a no-edit adapter-contract audit before implementation.
  The key decision is whether the Be2/PQS realization can honestly produce an
  existing `:pqs_multilayer_shell_source_plan`, or whether the route needs a
  new diatomic source-plan object plus a later consumer seam.

-- repo-manager@macmini
