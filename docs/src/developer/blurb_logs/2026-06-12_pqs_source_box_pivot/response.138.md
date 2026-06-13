Pass 138 response - add blocked diatomic PQS source-plan payload

Implemented and committed:

- commit: `288e0231 Add diatomic PQS source plan payload`

Files changed:

- `src/pqs_source_box_route_driver_helpers.jl`
- `test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl`

Payload/helper names:

- private payload type: `_PQSDiatomicCompleteCoreShellSourcePlanPayload`
- source-plan helper:
  `_pqs_source_box_route_driver_diatomic_complete_core_shell_source_plan_payload`
- small shared summary helper:
  `_pqs_source_box_route_driver_diatomic_route_skeleton_summary`
- Coulomb expansion summary helper:
  `_pqs_source_box_route_driver_diatomic_coulomb_expansion_summary`
- assembly field:
  `diatomic_complete_core_shell_source_plan_payload`

Wiring:

- `cartesian_assembly(...)` now builds
  `diatomic_complete_core_shell_source_plan_payload`.
- The existing private Ham readiness helper now receives that payload and
  reuses its route/source/retained/pair summary state instead of rebuilding the
  same adjacent summary group.
- The report-facing Ham readiness blocker remains unchanged:
  `:missing_diatomic_complete_core_shell_source_plan_producer`.

Default Be2/PQS source-plan payload:

- `status = :blocked_diatomic_complete_core_shell_source_plan`
- `blocker = :missing_parent_axis_bundle_object`
- `source_plan === nothing`
- `source_plan_status = :not_materialized_diatomic_complete_core_shell_source_plan`
- missing objects include:
  - `:parent_axis_bundle_object`
  - `:diatomic_complete_core_shell_source_realization_contract`

Probe-enabled Be2/PQS source-plan payload:

- `status = :blocked_diatomic_complete_core_shell_source_plan`
- `blocker = :missing_diatomic_complete_core_shell_source_realization_contract`
- `parent_axis_bundle_object_available = true`
- `source_plan === nothing`
- `source_plan_status = :not_materialized_diatomic_complete_core_shell_source_plan`
- available objects include `:parent_axis_bundle_object`
- missing objects include
  `:diatomic_complete_core_shell_source_realization_contract`
- missing objects no longer include `:parent_axis_bundle_object`

No materialization claims:

- `source_plan_materialized = false`
- `final_basis_materialized = false`
- `h1_materialized = false`
- `h1_j_materialized = false`
- `ham_payload_materialized = false`
- `route_driver_public_surface = false`
- `exports_materialized = false`
- `artifacts_materialized = false`
- metadata keeps
  `returns_pqs_multilayer_shell_source_plan = false`

Validation:

- `julia --project=. test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl`
  passed.
  - precompile dominated the wall time at about 57 s.
  - focused test summaries:
    - `Be2 PQS Ham payload readiness fingerprint | Pass 64 Total 64 Time 12.2s`
    - `Be2 PQS probe-enabled Ham readiness fingerprint | Pass 42 Total 42 Time 17.7s`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  passed with `load ok`.
- `git diff --check` passed.

Git status:

```text
## main...origin/main [ahead 1]
```

Deletion/shrinkage report:

- deleted: none in this boundary-creation pass.
- simplified: the diatomic complete core/shell source-plan blocker is now a
  route-owned payload, and Ham readiness can summarize that payload instead of
  coordinating a separate ad hoc field group.
- quarantined: final-basis, H1, H1/J, Ham materialization, RHF/SCF/Fock, WL
  payloads, public APIs, exports, artifacts, hfdmrg, and CR2 remain outside
  this private source-plan boundary.
- not deleted because: the existing Ham readiness payload/report fields remain
  compatibility surfaces for the focused Be2 fingerprint and for manager review.
- exact remaining caller/blocker: no private diatomic complete core/shell source
  realization contract/materializer currently produces an honest
  `:pqs_multilayer_shell_source_plan` consumer shape.

-- repo-doer@macmini
