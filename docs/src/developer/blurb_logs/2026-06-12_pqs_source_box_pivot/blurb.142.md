Pass 142 - add private diatomic raw-box route payload

Purpose:

Wrap the existing private raw product-box route producer in a compact private
driver payload for Be2/diatomic PQS. This should make raw product-box plans,
PQS axis-local coefficient objects, product/doside unit, and pair inventory
available as private candidate data without materializing a source plan or
promoting the producer to public route authority.

Why now:

Pass 141 proved the existing private producer can build the raw objects needed
for the next source-realization step:

```text
raw product-box plans
raw PQS/source coefficient objects
product/doside unit
pair inventory
```

The next useful boundary is a private payload that carries those objects and a
compact summary, so later passes do not keep calling the producer ad hoc from
tests.

Implementation target:

- Work in `src/pqs_source_box_route_driver_helpers.jl`.
- Add one private payload type and helper. Suggested names:
  - `_PQSDiatomicRawBoxRoutePayload`
  - `_pqs_source_box_route_driver_diatomic_raw_box_route_payload`
- Use the underlying producer:
  `CartesianContractedParentMetrics._pqs_pqs_product_raw_box_route_producer`
- Do not use the metadata-only raw product-box probe result as route authority.
- Derive axis metrics from the parent axis bundle with a small private helper,
  using the same PGDG fields used in the pass-141 test fingerprint.
- Wire the payload into `cartesian_assembly(...)` as one private object field,
  for example `diatomic_raw_box_route_payload`.
- Let the pass-138 source-plan payload summarize this raw-box payload status if
  that keeps the blocker chain clearer.

Expected statuses:

Default Be2/PQS:

```text
status = :blocked_diatomic_raw_box_route_payload
blocker = :missing_parent_axis_bundle_object
```

Probe-enabled Be2/PQS:

```text
status = :available_diatomic_raw_box_route_payload
blocker = nothing
producer_status = :private_shadow_only
```

Payload content:

Keep object fields and summaries compact:

- `status`
- `blocker`
- `producer`
- `producer_status`
- `descriptor_summary`
- `raw_product_box_plan_summary`
- `raw_pqs_plan_summary`
- `product_unit_summary`
- `pair_inventory_summary`
- `support_window_payload_status`
- `available_objects`
- `missing_objects`
- `summary`
- `metadata`

For the available probe-enabled path, summarize:

- left/right raw product-box plan object kinds and source-mode dims;
- left/right axis-local coefficient matrix shapes;
- left/right raw PQS plan representation and boundary selected counts;
- product unit kind, support count, coefficient matrix shape;
- pair inventory count/family facts if available;
- private/shadow status and nonclaim flags.

Nonclaims:

- `private_candidate_only = true`
- `raw_product_box_probe_authority = false`
- `source_plan_materialized = false`
- `final_basis_materialized = false`
- `h1_materialized = false`
- `h1_j_materialized = false`
- `ham_payload_materialized = false`
- `route_driver_public_surface = false`
- `exports_materialized = false`
- `artifacts_materialized = false`

Test target:

- Update only:
  `test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl`
- Replace or reduce direct producer assertions from pass 141 if the new payload
  makes them redundant. The test should prefer the route-owned private payload
  over ad hoc producer calls.
- Assert default blocked status and probe-enabled available status.
- Assert the existing support-window, source-plan, and Ham blockers remain
  behavior-preserving.

Decision rules:

- Do not materialize or fake `:pqs_multilayer_shell_source_plan`.
- Do not build final basis or H1.
- Do not make probe-enabled parent construction the default.
- If wrapping the producer requires broad route rewiring, stop and report.
- If the producer returns unexpected shape/status, preserve the blocker rather
  than patching with scalar aliases.

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
  `Add diatomic raw box route payload`
- If the focused test takes more than 60 seconds, report elapsed time and
  whether precompilation dominated.

Report back:

- Commit SHA if committed.
- Exact payload/helper names.
- Default/probe-enabled payload status/blocker.
- Raw-product-box/PQS/product-unit/pair-inventory compact facts.
- Confirmation that source plan/final basis/H1/H1-J/Ham data were not
  materialized and no public/export/artifact behavior was added.
- Validation commands/results.
- Git status.
- Deletion/shrinkage report:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

-- repo-manager@macmini
