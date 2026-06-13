Pass 140 - add diatomic PQS support-window payload

Purpose:

Add a compact private support-window/order payload for Be2/diatomic PQS. This
should record the source-box windows, support counts, retained order, candidate
support order, and required retained-to-support permutation without building raw
product-box plans, coefficients, final basis, H1, H1/J, Ham data, RHF, or
exports.

Why now:

Pass 139 found that Be2 retained units are metadata-only, while source boxes
carry enough x/y/z ranges to define support windows against the parent. It also
found an important ordering issue: current retained order is
`(:pqs_left, :pqs_right, :product)`, while the plausible complete-core/shell
support order is `(:product, :pqs_left, :pqs_right)`. That needs to become an
explicit route-owned fact before any source-plan materializer is attempted.

Implementation target:

- Work in `src/pqs_source_box_route_driver_helpers.jl`.
- Add one private payload type and helper. Suggested names:
  - `_PQSDiatomicCompleteCoreShellSupportWindowPayload`
  - `_pqs_source_box_route_driver_diatomic_complete_core_shell_support_window_payload`
- Wire it into `cartesian_assembly(...)` as one private object field, for
  example `diatomic_complete_core_shell_support_window_payload`.
- Also make the pass-138 source-plan payload reference or summarize the support
  window payload, so the blocker chain is:

```text
support windows/order available
-> source realization contract still missing
-> source plan still not materialized
```

Payload facts:

For probe-enabled Be2/PQS, assert and report:

- `status = :available_diatomic_complete_core_shell_support_windows`
- `parent_dims` or parent axis counts are available.
- `source_box_windows` include `:pqs_left`, `:product`, and `:pqs_right` in
  x/y/z range form.
- `source_mode_dims` are recorded from the route skeleton.
- `retained_order = (:pqs_left, :pqs_right, :product)`
- `candidate_core_then_shell_support_order = (:product, :pqs_left, :pqs_right)`
- `retained_to_support_order_permutation_required === true`
- `support_counts` are compact counts only, not large support-state arrays.
- `support_states_materialized = false`
- `raw_product_box_plans_materialized = false`
- `source_coefficients_materialized = false`

For default Be2/PQS, accept either:

- available metadata-only windows if parent dimensions/source boxes are enough;
  or
- blocked on `:missing_parent_axis_bundle_object` if the helper chooses to
  require the parent axis bundle for parent dimensions.

Keep the choice explicit in the test.

Missing objects should include:

- `:raw_product_box_plan_objects`
- `:pqs_axis_local_coefficients`
- `:diatomic_complete_core_shell_source_plan_materializer`

Test target:

- Update only:
  `test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl`
- Extend both default and probe-enabled fixture assertions for the support-window
  payload.
- Keep the existing Ham and source-plan blockers behavior-preserving.

Decision rules:

- Do not materialize raw product-box plans.
- Do not use raw product-box probe output as route authority.
- Do not create large support-state arrays in summaries.
- Do not return or fake a `:pqs_multilayer_shell_source_plan`.
- If deriving source-box windows requires more than compact range/count
  extraction, stop and report.

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
  `Add diatomic PQS support window payload`
- If the focused test takes more than 60 seconds, report elapsed time and
  whether precompilation dominated.

Report back:

- Commit SHA if committed.
- Exact payload/helper names.
- Default and probe-enabled support-window status/blocker.
- Source-box windows/order/permutation facts.
- Confirmation that raw product-box plans, support states, coefficients, source
  plan, final basis, H1, H1/J, and Ham data were not materialized.
- Validation commands/results.
- Git status.
- Deletion/shrinkage report:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

-- repo-manager@macmini
