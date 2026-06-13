Pass 146 response - add private diatomic source-plan object

Commit:
- `7dc7d654 Add diatomic source plan object`

Files changed:
- `src/pqs_source_box_route_driver_helpers.jl`
- `test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl`

Exact source-plan object/helper names:
- `_PQSDiatomicCompleteCoreShellSourcePlan`
- `_pqs_source_box_route_driver_diatomic_complete_core_shell_source_plan`
- Supporting private helpers:
  - `_pqs_source_box_route_driver_raw_plan_support_rows`
  - `_pqs_source_box_route_driver_raw_plan_support_coefficients`

Default Be2/PQS source-plan status:
- `source_plan === nothing`
- `source_plan_status = :not_materialized_diatomic_complete_core_shell_source_plan`
- payload status remains `:blocked_diatomic_complete_core_shell_source_plan`
- blocker remains `:missing_parent_axis_bundle_object`

Probe-enabled Be2/PQS source-plan status:
- source-plan object is carried by `source_plan_payload.source_plan`.
- `object_kind = :pqs_diatomic_complete_core_shell_source_plan`
- `status = :available_pqs_diatomic_complete_core_shell_source_plan`
- `blocker = nothing`
- It explicitly does not claim `:pqs_multilayer_shell_source_plan`.
- Source-plan payload now blocks on:
  - `:missing_diatomic_complete_core_shell_final_basis_consumer`
- It no longer blocks on:
  - `:missing_pqs_multilayer_shell_source_plan_adapter_contract`

Real structured data now carried:
- `bundles`
- `metrics`
- `core_unit_key = :product`
- `shell_unit_keys = (:pqs_left, :pqs_right)`
- `core_support_indices`
- `core_support_states`
- `shell_support_indices`
- `shell_support_states`
- `shell_final_coefficients`
- `support_order = (:product, :pqs_left, :pqs_right)`
- `route_retained_order = (:pqs_left, :pqs_right, :product)`
- `retained_pre_final_map`
- `source_unit_summaries`
- `convention_labels`
- compact `summary`
- `metadata`

Probe-enabled array/matrix facts:
- `length(core_support_indices) == 25`
- `length(core_support_states) == 25`
- `length(shell_support_indices) == 250`
- `length(shell_support_states) == 250`
- `length(unique(shell_support_indices)) == 250`
- core/shell support intersection is empty.
- `size(shell_final_coefficients) == (250, 196)`
- all shell coefficient entries are finite.
- shell coefficient block structure is `:block_diagonal_left_right_pqs`.
- precleanup retained dimension is `221`.

Behavior boundary confirmation:
- No `:pqs_multilayer_shell_source_plan` object is returned or faked.
- No final-basis, H1, H1-J, Ham, support-one-body consumer, support-density consumer, RHF/SCF, public API, export, artifact, hfdmrg, or CR2 behavior was added.
- Raw product-box probes, shell/support-row contraction, and old WL adapter paths were not promoted to route authority.
- Retained diagnostic/self-integral weights were not reinterpreted as IDA/quadrature weights.

Validation:
- Initial run of the focused test caught a local parse error in a range expression; fixed before final validation.
- `julia --project=. test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl`
  - passed after fix.
  - precompilation: `56686.3 ms  ✓ GaussletBases`
  - `Be2 PQS Ham payload readiness fingerprint | Pass 151 Total 151 Time 13.0s`
  - `Be2 PQS probe-enabled Ham readiness fingerprint | Pass 190 Total 190 Time 38.3s`
  - total wall time exceeded 60s because package precompilation took about 57s; test body time was about 51.3s.
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - passed: `load ok`
- `git diff --check`
  - passed.

Git status after commit:
```text
## main...origin/main [ahead 1]
```

Deletion/shrinkage report:
- deleted: none; this pass added the private source-plan object needed before shrinking the old adapter blocker.
- simplified: the probe-enabled source-plan payload now has a real private diatomic source plan and a sharper next blocker, replacing the old missing multilayer adapter contract for this route.
- quarantined: final-basis/H1/H1-J/Ham materialization; support-one-body and support-density consumers; old `:pqs_multilayer_shell_source_plan` adapter claim; RHF/SCF/DIIS; WL payloads; public API/exports/artifacts/hfdmrg/CR2.
- not deleted because: one-center `pqs_multilayer_shell_source_plan` remains the validated compact H1 fixture contract, and no diatomic final-basis consumer has been approved yet.
- exact remaining caller/blocker: Be2/PQS now has `_PQSDiatomicCompleteCoreShellSourcePlan`, but no approved consumer seam calls lower-level final-basis/H1 helpers for `:pqs_diatomic_complete_core_shell_source_plan`; source-plan payload blocks on `:missing_diatomic_complete_core_shell_final_basis_consumer`.

-- repo-doer@macmini
