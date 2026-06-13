Pass 144 response - add diatomic PQS source-realization payload

Commit:
- `7a79d25a Add diatomic source realization payload`

Files changed:
- `src/pqs_source_box_route_driver_helpers.jl`
- `test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl`

Payload/helper names:
- `_PQSDiatomicCompleteCoreShellSourceRealizationPayload`
- `_pqs_source_box_route_driver_diatomic_complete_core_shell_source_realization_payload`
- New `cartesian_assembly` field:
  - `diatomic_complete_core_shell_source_realization_payload`

Default Be2/PQS realization status:
- `status = :blocked_diatomic_complete_core_shell_source_realization`
- `blocker = :missing_parent_axis_bundle_object`
- `support_window_payload_status = :available_diatomic_complete_core_shell_support_windows`
- `raw_box_route_payload_status = :blocked_diatomic_raw_box_route_payload`

Probe-enabled Be2/PQS realization status:
- `status = :available_diatomic_complete_core_shell_source_realization`
- `blocker = nothing`
- `object_kind_claim = :not_pqs_multilayer_shell_source_plan`
- `bundles_role = :parent_axis_bundle`

Realization mapping facts:
- `core_unit_key = :product`
- `shell_unit_keys = (:pqs_left, :pqs_right)`
- `retained_order = (:pqs_left, :pqs_right, :product)`
- `support_order = (:product, :pqs_left, :pqs_right)`
- `retained_to_support_order_permutation_required = true`
- route retained ranges are copied from `producer.descriptor.expected_ranges` when the raw-box payload is available.
- source-plan precleanup ranges on the probe-enabled path:
  - `product = 1:25`
  - `pqs_left = 26:123`
  - `pqs_right = 124:221`

Coefficient/count facts:
- `core_support_count = 25`
- `shell_support_counts = (pqs_left = 125, pqs_right = 125)`
- `shell_support_count = 250`
- `shell_retained_counts = (pqs_left = 98, pqs_right = 98)` on the available path.
- `shell_retained_count = 196`
- `precleanup_retained_dimension = 221`
- `shell_final_coefficients_shape = (250, 196)`
- `shell_coefficient_block_structure = :block_diagonal_left_right_pqs`

Source-plan payload wiring:
- `_PQSDiatomicCompleteCoreShellSourcePlanPayload` now carries `source_realization_payload`.
- Its summary includes `source_realization_payload_status`.
- It still keeps `source_plan === nothing`.
- Probe-enabled source-plan status still blocks, now on:
  - `:missing_pqs_multilayer_shell_source_plan_adapter_contract`

Behavior boundary confirmation:
- No `:pqs_multilayer_shell_source_plan` object is returned or faked.
- No source-plan, final-basis, H1, H1-J, Ham, RHF/SCF, public API, export, or artifact behavior was added.
- No support matrices or shell coefficient matrices were materialized; the new payload records shape/block-structure metadata only.
- Raw product-box probes and shell/support-row contraction remain non-authoritative.
- Retained diagnostic/self-integral weights were not reinterpreted as IDA/quadrature weights.

Validation:
- `julia --project=. test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl`
  - passed.
  - precompilation: `57138.7 ms  ✓ GaussletBases`
  - `Be2 PQS Ham payload readiness fingerprint | Pass 151 Total 151 Time 12.9s`
  - `Be2 PQS probe-enabled Ham readiness fingerprint | Pass 154 Total 154 Time 36.8s`
  - total wall time was roughly 99s; the run exceeded 60s because this edit triggered package precompilation, which dominated the startup portion.
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - passed: `load ok`
- `git diff --check`
  - passed.

Git status after commit:
```text
## main...origin/main [ahead 1]
```

Deletion/shrinkage report:
- deleted: none; this pass introduced the private boundary before deleting or shrinking old surface.
- simplified: `diatomic_complete_core_shell_source_plan_payload` now consumes one source-realization status/payload rather than only carrying a broad missing realization contract.
- quarantined: `:pqs_multilayer_shell_source_plan` return claim; final-basis/H1/H1-J/Ham materialization; support-row contraction authority; RHF/SCF; WL payloads; public API/exports/artifacts.
- not deleted because: existing complete-core/shell H1/final-basis helpers remain the validated one-center consumer shape, and Be2/PQS still needs an approved adapter/new-object contract.
- exact remaining caller/blocker: `pqs_multilayer_complete_core_shell_final_basis`, `pqs_multilayer_complete_core_shell_h1_payload`, support one-body, and support density helpers still require an available `object_kind = :pqs_multilayer_shell_source_plan`; the new Be2 payload now blocks specifically on `:missing_pqs_multilayer_shell_source_plan_adapter_contract`.

-- repo-doer@macmini
