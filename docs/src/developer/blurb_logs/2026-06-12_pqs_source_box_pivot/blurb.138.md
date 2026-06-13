Pass 138 - add blocked diatomic PQS source-plan payload

Purpose:

Add a private blocked source-plan payload for the Be2/diatomic PQS route. This
should narrow the remaining blocker to the missing diatomic source-realization
contract without building final basis, H1, H1/J, Ham data, WL payloads, RHF, or
exports.

Why now:

Pass 137 found that the existing final-basis/H1 consumers require a genuine
available `:pqs_multilayer_shell_source_plan` shape. Be2/PQS has route skeleton,
retained units, pair inventory, centers, and parent axis bundle in the
probe-enabled path, but it does not yet have an honest source realization that
satisfies that consumer contract.

Implementation target:

- Work in `src/pqs_source_box_route_driver_helpers.jl`.
- Add one private payload type and helper. Suggested names:
  - `_PQSDiatomicCompleteCoreShellSourcePlanPayload`
  - `_pqs_source_box_route_driver_diatomic_complete_core_shell_source_plan_payload`
- Wire it into `cartesian_assembly(...)` as one private object field, for
  example `diatomic_complete_core_shell_source_plan_payload`.
- Reuse existing diatomic summary helpers from pass 135. Do not copy the same
  field group into another ad hoc summary if a helper already provides it.
- Update the diatomic Ham readiness payload to reference or summarize the new
  source-plan payload if that keeps the state clearer.

Expected current statuses:

For default Be2/PQS:

```text
source_plan_payload.status = :blocked_diatomic_complete_core_shell_source_plan
source_plan_payload.blocker = :missing_parent_axis_bundle_object
```

For probe-enabled Be2/PQS:

```text
source_plan_payload.status = :blocked_diatomic_complete_core_shell_source_plan
source_plan_payload.blocker =
    :missing_diatomic_complete_core_shell_source_realization_contract
```

The exact label can be `:missing_diatomic_complete_core_shell_source_plan_materializer`
if that better fits the implementation, but it must mean "we have structured
route inputs, but no materializer/realization contract yet."

Payload content:

Keep the object compact:

- `status`
- `blocker`
- `route_family`
- `system_classification`
- `bond_axis`
- `parent_axis_bundle_object_available`
- `route_skeleton_summary`
- `source_box_summary`
- `retained_unit_summary`
- `pair_inventory_summary`
- `center_summary`
- `coulomb_expansion_summary`
- `source_plan`
- `source_plan_status`
- `available_objects`
- `missing_objects`
- `summary`
- `metadata`

Nonclaim flags:

- `source_plan_materialized = false`
- `final_basis_materialized = false`
- `h1_materialized = false`
- `h1_j_materialized = false`
- `ham_payload_materialized = false`
- `route_driver_public_surface = false`
- `exports_materialized = false`
- `artifacts_materialized = false`

Test target:

- Update only the focused Be2 fingerprint test:
  `test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl`
- Assert for the default fixture:
  - source-plan payload exists;
  - parent axis bundle is unavailable;
  - blocker is `:missing_parent_axis_bundle_object`;
  - no materialization flags are true.
- Assert for the probe-enabled fixture:
  - parent axis bundle is available;
  - source-plan payload exists;
  - blocker is the missing diatomic source-realization/materializer contract;
  - existing private Ham payload remains blocked;
  - no final-basis/H1/H1-J/Ham/public/export/artifact claims are made.

Decision rules:

- Do not return `object_kind = :pqs_multilayer_shell_source_plan` unless the
  object actually satisfies the existing consumer shape. In this pass it should
  not.
- Do not build final basis or H1.
- Do not make probe-enabled parent construction the default.
- If this requires more than one compact payload and one focused test update,
  stop and report the smaller blocker instead of expanding the implementation.

Trust boundary:

- Private/internal only.
- No final-basis/H1/H1-J/Ham materialization.
- No RHF/SCF work.
- No WL payload, public API, exports, artifacts, hfdmrg, or CR2 execution.
- No scalar report-field clouds.
- Do not force one-center shellification/support-row semantics onto the
  diatomic source-box route.
- Do not reinterpret retained diagnostic/self-integral weights as
  IDA/quadrature weights.
- Do not request interactive command approval during unattended baton work. If
  approval would be required, write `.agent_handoffs/ATTENTION.md` and stop.

Validation:

- `julia --project=. test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git diff --check`
- If validation passes, commit with a clear message such as:
  `Add diatomic PQS source plan payload`
- If the focused test takes more than 60 seconds, report elapsed time and
  whether precompilation dominated.

Report back:

- Commit SHA if committed.
- Exact payload/helper names.
- Default and probe-enabled source-plan payload status/blocker/missing objects.
- Confirmation that no source plan/final basis/H1/Ham materialization occurred.
- Validation commands/results.
- Git status.
- Deletion/shrinkage report:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

-- repo-manager@macmini
