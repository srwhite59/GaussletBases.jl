Pass 147 response - add private diatomic final-basis payload

Commit:
- `8339dccd Add diatomic final basis payload`

Files changed:
- `src/pqs_source_box_route_driver_helpers.jl`
- `test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl`

Exact payload/helper names:
- `_PQSDiatomicCompleteCoreShellFinalBasisPayload`
- `_pqs_source_box_route_driver_diatomic_complete_core_shell_final_basis_payload`
- New `cartesian_assembly` field:
  - `diatomic_complete_core_shell_final_basis_payload`

Default Be2/PQS final-basis status:
- `status = :blocked_diatomic_complete_core_shell_final_basis_payload`
- `blocker = :missing_pqs_diatomic_complete_core_shell_source_plan`
- `source_plan === nothing`
- `final_basis === nothing`
- `final_basis_status = :not_materialized_diatomic_complete_core_shell_final_basis`

Probe-enabled Be2/PQS final-basis status:
- `status = :available_diatomic_complete_core_shell_final_basis_payload`
- `blocker = nothing`
- `source_plan_status = :available_pqs_diatomic_complete_core_shell_source_plan`
- `final_basis_status = :available_pqs_complete_core_shell_final_basis`
- Final basis is built by calling `CartesianFinalBasisRealization.pqs_complete_core_shell_final_basis(...)` directly.
- The pass does not call `pqs_multilayer_complete_core_shell_final_basis(...)`.

Observed final-basis facts:
- `final_retained_count = 221`
- `precleanup_retained_dimension = 221`
- `core_support_count = 25`
- `shell_support_count = 250`
- `shell_final_retained_count = 196`
- `support_row_order = :core_then_shell`
- `old_source_plan_object_kind = false`
- source-plan support order remains `(:product, :pqs_left, :pqs_right)`
- shell coefficient block structure remains `:block_diagonal_left_right_pqs`

Readiness/Ham fingerprint:
- Probe-enabled readiness now has:
  - `:diatomic_complete_core_shell_final_basis` in `available_objects`
  - `blocker = :missing_diatomic_complete_core_shell_h1_consumer`
  - `:diatomic_complete_core_shell_h1_consumer` in `missing_objects`
- It no longer reports `:diatomic_complete_core_shell_final_basis_consumer` as missing once the final basis is available.

Behavior boundary confirmation:
- No `:pqs_multilayer_shell_source_plan` object is returned or faked.
- No H1, H1-J, Ham, support-one-body consumer, support-density consumer, RHF/SCF, public API, export, artifact, hfdmrg, or CR2 behavior was added.
- Only overlap blocks needed by the lower final-basis helper are materialized.
- Raw product-box probes, shell/support-row contraction, and old WL adapter paths were not promoted to route authority.
- Retained diagnostic/self-integral weights were not reinterpreted as IDA/quadrature weights.

Validation:
- `julia --project=. test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl`
  - passed.
  - precompilation: `56745.7 ms  ✓ GaussletBases`
  - `Be2 PQS Ham payload readiness fingerprint | Pass 165 Total 165 Time 13.2s`
  - `Be2 PQS probe-enabled Ham readiness fingerprint | Pass 221 Total 221 Time 39.1s`
  - total wall time exceeded 60s because package precompilation took about 57s; test body time was about 52.3s.
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - passed: `load ok`
- `git diff --check`
  - passed.

Git status after commit:
```text
## main...origin/main [ahead 1]
```

Deletion/shrinkage report:
- deleted: none; this pass added the private final-basis payload seam.
- simplified: the probe-enabled readiness path now advances past the missing final-basis consumer and reports a sharper missing H1 consumer.
- quarantined: H1/H1-J/Ham materialization; support-one-body and support-density consumers beyond final-basis overlap blocks; old `:pqs_multilayer_shell_source_plan` adapter claim; RHF/SCF/DIIS; WL payloads; public API/exports/artifacts/hfdmrg/CR2.
- not deleted because: the one-center final-basis/H1 stack remains the validated compact H1 fixture contract, and a diatomic H1 consumer has not been approved yet.
- exact remaining caller/blocker: Be2/PQS now has a private source plan and private final basis; the next missing seam is `:missing_diatomic_complete_core_shell_h1_consumer`.

-- repo-doer@macmini
