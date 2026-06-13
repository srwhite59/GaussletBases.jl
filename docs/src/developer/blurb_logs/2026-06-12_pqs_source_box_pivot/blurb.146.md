Pass 146 - add private diatomic source-plan object

Purpose:

Implement the private Be2/PQS source-plan object identified by pass 145. This
should turn the probe-enabled Be2/PQS path from "source realization available,
old source-plan adapter missing" into "private diatomic source plan available,
final-basis/H1 consumer still missing."

Physics target:

- Medium-term target: make the PQS driver usable enough for a downstream CR2
  agent to compare Be2 WL and PQS Hamiltonian construction.
- This pass only advances the structured source-plan seam. It is not a physics
  endpoint and does not run CR2 or hfdmrg.

Implementation target:

- Work in `src/pqs_source_box_route_driver_helpers.jl`.
- Add one private source-plan object and materializer. Suggested names:
  - `_PQSDiatomicCompleteCoreShellSourcePlan`
  - `_pqs_source_box_route_driver_diatomic_complete_core_shell_source_plan`
- Build it from the existing support-window, raw-box route, and
  source-realization payloads.
- Wire it into `_PQSDiatomicCompleteCoreShellSourcePlanPayload` as its
  `source_plan` when available.

Object contract:

The new source-plan object should claim:

```text
object_kind = :pqs_diatomic_complete_core_shell_source_plan
status = :available_pqs_diatomic_complete_core_shell_source_plan
```

It must not claim:

```text
object_kind = :pqs_multilayer_shell_source_plan
```

Carry real structured data, not only counts, where available:

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
- a compact retained/pre-final range or permutation map
- source unit references or compact summaries
- convention labels and nonpromotion flags
- compact `summary`
- `metadata`

Expected probe-enabled facts:

- real `core_support_indices` length `25`
- real `core_support_states` length `25`
- real `shell_support_indices` length `250`
- real `shell_support_states` length `250`
- `shell_final_coefficients` size `(250, 196)`
- shell coefficient block structure `:block_diagonal_left_right_pqs`
- precleanup retained dimension `221`
- support order `(:product, :pqs_left, :pqs_right)`
- route retained order `(:pqs_left, :pqs_right, :product)`
- old source-plan object kind is false
- final-basis/H1/H1-J/Ham materialization flags remain false

Source-plan payload behavior:

- Default Be2/PQS should remain blocked upstream on missing parent axis-bundle
  / raw-box route payload.
- Probe-enabled Be2/PQS should now carry the private diatomic source-plan object.
- The source-plan payload should no longer block on
  `:missing_pqs_multilayer_shell_source_plan_adapter_contract`.
- Use a sharper next blocker such as:

```text
:missing_diatomic_complete_core_shell_final_basis_consumer
```

or an equally precise label if the existing code vocabulary suggests a better
name.

Decision rules:

- If the raw-box route payload does not expose real support arrays or local
  coefficient matrices, stop and report the exact missing field instead of
  inventing data from counts.
- If `bundles` or `metrics` provenance is ambiguous, keep the source plan
  blocked with a precise blocker.
- Do not fill old one-center fields such as `core_box`, `outer_box`,
  `layer_count`, or `shell_records` just to satisfy old consumers.
- Do not call `pqs_multilayer_complete_core_shell_final_basis(...)` in this
  pass.

Test target:

- Update only:
  `test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl`
- Assert default blocked behavior remains behavior-preserving.
- Assert probe-enabled private source-plan availability and real array/matrix
  shapes.
- Assert the object kind is `:pqs_diatomic_complete_core_shell_source_plan`.
- Assert it is not `:pqs_multilayer_shell_source_plan`.
- Assert source-plan/final-basis/H1/H1-J/Ham/global/public/export/artifact flags
  remain false.

Trust boundary:

- Private/internal only.
- No final-basis/H1/H1-J/Ham materialization.
- No support-one-body or support-density consumers yet.
- No RHF/SCF/DIIS work.
- No WL payload, public API, exports, artifacts, hfdmrg, or CR2 execution.
- No scalar report-field clouds.
- Do not promote shell/support-row contraction, raw product-box probes, or old
  WL adapter paths to route authority.
- Do not reinterpret retained diagnostic/self-integral weights as
  IDA/quadrature weights.
- Do not request interactive command approval during unattended baton work. If
  approval would be required, write `.agent_handoffs/ATTENTION.md` and stop.

Validation:

- `julia --project=. test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git diff --check`
- If the focused test takes more than 60 seconds, report elapsed time and
  whether precompilation dominated.

Report back:

- Commit SHA if committed.
- Exact source-plan object/helper names.
- Default/probe-enabled source-plan status/blocker.
- Real support array lengths and shell coefficient matrix shape.
- Confirmation that no final-basis/H1/H1-J/Ham/public/export/artifact behavior
  was added.
- Validation commands/results.
- Git status.
- Deletion/shrinkage report:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

-- repo-manager@macmini
