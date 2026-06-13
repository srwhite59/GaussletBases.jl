Pass 135 - add diatomic PQS Ham readiness payload

Purpose:

Add a compact private readiness payload for the Be2/diatomic PQS complete-core
shell Ham seam. This payload should make the current blocker structured and
route-owned without pretending the Be2 route can already build final basis, H1,
H1/J, or a Ham payload.

Why now:

Pass 134 found that the one-center complete-core/shell producer is structurally
terminal-shellification backed, while Be2/PQS is a source-box-first diatomic
route with `(:pqs_left, :product, :pqs_right)` units. The correct next step is a
private readiness object that inventories Be2/PQS route facts and the missing
diatomic source-plan producer.

Implementation target:

- Work in `src/pqs_source_box_route_driver_helpers.jl`.
- Add one private compact payload type and helper. Suggested names:
  - `_PQSDiatomicCompleteCoreShellHamReadinessPayload`
  - `_pqs_source_box_route_driver_diatomic_complete_core_shell_ham_readiness_payload`
- Wire it into the private assembly/diagnostic surface next to
  `complete_core_shell_diagnostic_route_payload`.
- Prefer a single object field, not scalar report fields.
- If the route is not PQS diatomic, return a compact not-applicable payload
  rather than throwing.

Payload content:

Keep fields compact and object/fingerprint oriented. Include at most:

- `status`
- `blocker`
- `route_family`
- `system_classification`
- `bond_axis`
- `center_summary`
- `parent_axis_bundle_object_available`
- `route_skeleton_summary`
- `source_box_summary`
- `retained_unit_summary`
- `pair_inventory_summary`
- `available_objects`
- `missing_objects`
- `summary`
- `metadata`

Expected current Be2/PQS status:

```text
status = :blocked_diatomic_complete_core_shell_ham_readiness
blocker = :missing_diatomic_complete_core_shell_source_plan_producer
```

The exact names can vary if the codebase has a better local label, but keep the
meaning precise: Be2/PQS has structured route skeleton/source-box facts, but no
route-owned diatomic complete-core/shell source-plan producer.

Test target:

- Extend the focused standalone Be2 fingerprint test:
  `test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl`
- Assert:
  - Be2 route reaches `cartesian_assembly(...)`;
  - the new readiness payload exists;
  - readiness status/blocker match the current missing producer state;
  - route family is `:pqs_source_box`;
  - system classification is `:bond_aligned_diatomic`;
  - bond axis is `:x`;
  - center/source-box/retained-unit/pair-inventory summaries are present and
    compact;
  - parent axis-bundle availability reflects the current Be2 fixture path;
  - existing `complete_core_shell_ham_payload` blocker remains behavior
    preserving.

Trust boundary:

- Private/internal only.
- No public API, exports, artifacts, hfdmrg, or CR2 execution.
- No WL payload implementation.
- No final-basis, H1, H1/J, density-interaction, RHF, or SCF work.
- No route-global materialization or physics endpoint promotion.
- No scalar report-field clouds.
- Do not make shell/support-row contraction route authority.
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
- Exact payload/helper names.
- Readiness status/blocker and compact summary facts observed in the Be2 test.
- Confirmation that the existing Ham payload blocker is unchanged.
- Validation commands/results.
- Git status.
- Deletion/shrinkage report:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

-- repo-manager@macmini
