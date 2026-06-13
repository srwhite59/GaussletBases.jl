Pass 151 review - accepted

Commit reviewed:
- `8c2186ef Add diatomic Ham input payload`

Summary:
- Adds `_PQSDiatomicCompleteCoreShellHamInputPayload`.
- Adds `_pqs_source_box_route_driver_diatomic_complete_core_shell_ham_input_payload(...)`.
- Wires the payload into `cartesian_assembly(...)` as
  `diatomic_complete_core_shell_ham_input_payload`.
- Builds private Be2/PQS support-density inputs:
  - density provenance
  - support weights
  - raw pair factor terms
  - support raw pair numerator
  - pre-final density interaction

Review result:
- Accepted.
- The route now has private Be2/PQS H1 plus reusable electron-electron/Ham-input
  data.
- The implementation does not call the old one-center H1-J payload, does not
  broaden old density-helper guards, and does not claim a full Ham payload or
  export.

Observed probe-enabled facts:
- Ham-input status:
  `:available_diatomic_complete_core_shell_ham_input_payload`
- Density gauge:
  `:pre_final_localized_positive_weight`
- Raw-pair convention:
  `:raw_numerator`
- Final dimension:
  `221`
- Pre-final dimension:
  `221`
- Support weight count:
  `275`
- Support raw-pair shape:
  `(275, 275)`
- Pre-final pair matrix shape:
  `(221, 221)`
- Readiness blocker advanced to:
  `:missing_diatomic_ham_consumer_contract`

Validation checked by doer and repeated by manager:
- `julia --project=. test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl`
  - manager repeat passed:
    - `Be2 PQS Ham payload readiness fingerprint | Pass 202 Total 202 Time 13.3s`
    - `Be2 PQS probe-enabled Ham readiness fingerprint | Pass 296 Total 296 Time 40.3s`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - passed: `load ok`
- `git diff --check HEAD~1..HEAD`
  - passed

Boundary status:
- Good. No H1-J scalar, full Ham payload, public API, export, artifact, RHF,
  WL payload, hfdmrg execution, or CR2 execution was added.

Next step:
- Audit the Ham consumer/handoff contract before implementation. The next
  object should be shaped by what a downstream Be2 WL/PQS comparison needs, not
  by the historical H1-J scalar diagnostic.

-- repo-manager@macmini
