Pass 144 - add diatomic PQS source-realization payload

Purpose:

Add a private diatomic complete-core/shell source-realization payload for Be2/PQS.
This payload should record the mapping from raw-box route data to the proposed
core/body plus left/right source-sector realization, while explicitly not
claiming to be a `:pqs_multilayer_shell_source_plan`.

Why now:

Pass 143 found a plausible mapping:

```text
core/body sector: product unit, 25 support rows
shell/source sector: pqs_left, 125 support rows -> 98 retained columns
shell/source sector: pqs_right, 125 support rows -> 98 retained columns
precleanup retained dimension: 25 + 98 + 98 = 221
```

But this is a diatomic source-box realization, not the one-center nested
shell-plan vocabulary. The route needs an explicit private realization payload
before any adapter/materializer is attempted.

Implementation target:

- Work in `src/pqs_source_box_route_driver_helpers.jl`.
- Add one private payload type and helper. Suggested names:
  - `_PQSDiatomicCompleteCoreShellSourceRealizationPayload`
  - `_pqs_source_box_route_driver_diatomic_complete_core_shell_source_realization_payload`
- Inputs:
  - `parent`
  - `route_skeleton`
  - `recipe`
  - support-window payload
  - raw-box route payload
- Wire it into `cartesian_assembly(...)` as one private object field, for
  example `diatomic_complete_core_shell_source_realization_payload`.
- Let the pass-138 source-plan payload summarize this realization payload status
  and keep `source_plan === nothing`.

Expected statuses:

Default Be2/PQS:

```text
status = :blocked_diatomic_complete_core_shell_source_realization
blocker = :missing_diatomic_raw_box_route_payload
```

or, if the implementation chooses to surface the upstream blocker directly:

```text
blocker = :missing_parent_axis_bundle_object
```

Probe-enabled Be2/PQS:

```text
status = :available_diatomic_complete_core_shell_source_realization
blocker = nothing
object_kind_claim = :not_pqs_multilayer_shell_source_plan
```

Payload facts for the available path:

- `core_unit_key = :product`
- `shell_unit_keys = (:pqs_left, :pqs_right)`
- `retained_order = (:pqs_left, :pqs_right, :product)`
- `support_order = (:product, :pqs_left, :pqs_right)`
- `retained_to_support_order_permutation_required = true`
- route retained ranges from the descriptor, if available.
- source-plan precleanup ranges:
  - `product = 1:25`
  - `pqs_left = 26:123`
  - `pqs_right = 124:221`
- `core_support_count = 25`
- `shell_support_counts = (pqs_left = 125, pqs_right = 125)`
- `shell_support_count = 250`
- `shell_retained_counts = (pqs_left = 98, pqs_right = 98)`
- `shell_retained_count = 196`
- `precleanup_retained_dimension = 221`
- `shell_final_coefficients_shape = (250, 196)`
- `shell_coefficient_block_structure = :block_diagonal_left_right_pqs`
- `bundles_role = :parent_axis_bundle`
- `source_plan_materialized = false`
- `returns_pqs_multilayer_shell_source_plan = false`

Test target:

- Update only:
  `test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl`
- Assert default blocked status and probe-enabled available status.
- Assert the source-plan payload now summarizes the realization payload status
  but still blocks on the adapter/materializer contract.
- Keep existing raw-box/support-window/Ham blockers behavior-preserving.

Decision rules:

- Do not return or fake `object_kind = :pqs_multilayer_shell_source_plan`.
- Do not build support matrices, final basis, H1, H1/J, or Ham payloads.
- Keep coefficient data as shape/block-structure metadata only unless already
  carried by the raw-box payload.
- If range/permutation facts are ambiguous, stop and report instead of guessing.

Trust boundary:

- Private/internal only.
- No source-plan/final-basis/H1/H1-J/Ham materialization.
- No RHF/SCF work.
- No WL payload, public API, exports, artifacts, hfdmrg, or CR2 execution.
- No scalar report-field clouds.
- Do not promote raw product-box probes, shell/support-row contraction, or WL
  adapter paths to route authority.
- Do not reinterpret retained diagnostic/self-integral weights as
  IDA/quadrature weights.
- Do not request interactive command approval during unattended baton work. If
  approval would be required, write `.agent_handoffs/ATTENTION.md` and stop.

Validation:

- `julia --project=. test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git diff --check`
- If validation passes, commit with a clear message such as:
  `Add diatomic source realization payload`
- If the focused test takes more than 60 seconds, report elapsed time and
  whether precompilation dominated.

Report back:

- Commit SHA if committed.
- Exact payload/helper names.
- Default/probe-enabled realization status/blocker.
- Realization mapping facts and coefficient shape facts.
- Confirmation that no source-plan/final-basis/H1/H1-J/Ham/public/export/artifact
  behavior was added.
- Validation commands/results.
- Git status.
- Deletion/shrinkage report:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

-- repo-manager@macmini
