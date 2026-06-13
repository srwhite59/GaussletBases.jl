Pass 139 response - audit diatomic source-realization inputs

No source/test edits. No commits.

Files/helpers inspected:

- `src/CartesianContractedParentMetrics.jl`
  - `_pqs_pqs_product_source_box_route_skeleton`
  - `_pqs_explicit_core_spacing_route_raw_product_box_plan_probe`
  - `_pqs_raw_product_box_plan*`
  - `_pqs_pqs_product_safe_term_route_descriptor`
  - `_pqs_pqs_product_raw_box_route_producer`
  - `_pqs_product_doside_identity_slab_unit`
  - route-shaped density-density and nuclear-attraction consumers
- `src/pqs_source_box_route_driver_helpers.jl`
  - pass-138 `_PQSDiatomicCompleteCoreShellSourcePlanPayload`
  - raw product-box probe plumbing
  - parent-axis probe/object handoff
- `src/pqs_multilayer_shell_source_plan.jl`
  - `pqs_multilayer_shell_source_plan`
  - `_pqs_multilayer_realize_shell_source_plan`
- `src/pqs_multilayer_support_one_body.jl`
- `src/pqs_multilayer_support_density.jl`
- `src/cartesian_shellization_route.jl`
- `test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl`
- route-driver/report tests mentioning Be2 raw product-box and materializer
  blockers.

Retained-unit records:

- The current Be2/PQS route skeleton retained-unit records are metadata only.
- They carry `unit_key`, role/kind, source box, source dimensions, retained
  rule derivation, retained range/count, provenance, and weight semantics.
- They do not carry retained/source coefficient matrices.
- They do not carry support states or support indices.
- Actual coefficient data exists only after constructing raw product-box plans:
  `axis_local_coefficients`, `source_mode_indices`, and the boundary selector
  live in the raw product-box plan path.
- Product/doside identity slab support states and support indices can be built
  by `_pqs_product_doside_identity_slab_unit`, but that object is not carried by
  the current route-driver skeleton.

Source-box windows/states:

- The source boxes are sufficient to define windows:
  `pqs_left`, `pqs_right`, and `product` each carry x/y/z ranges plus source
  dimensions.
- With parent dimensions, those ranges are enough to derive parent support
  states and flat support indices for a window.
- Source boxes alone are not enough to define retained/source coefficients for
  PQS units. Those require the parent axis bundle and a real raw product-box
  plan construction.
- Therefore the current probe-enabled Be2 path has enough geometry to begin a
  support-window contract, but not enough data in the skeleton to honestly emit
  `:pqs_multilayer_shell_source_plan`.

Raw product-box probe machinery:

- `_pqs_explicit_core_spacing_route_raw_product_box_plan_probe` can construct
  raw product-box plans internally when requested, but it returns only
  `unit_plan_metadata`.
- It deliberately drops the raw product-box plan objects and reports itself as
  private development/probe-only, not production route authority.
- It is useful for checking whether raw-box plan construction is possible, but
  should remain probe-only.
- The stronger candidate for a future route-owned structured input is not the
  probe result; it is the underlying raw-box route construction contract around
  `_pqs_pqs_product_raw_box_route_producer`, which returns raw product-box
  plans, raw PQS plans, a product unit, a descriptor, and pair inventory. That
  producer is currently `:private_shadow_only`, so it still needs an explicit
  route-driver boundary before it can be treated as route-owned.

Density-density helper records:

- The source-box density-density helper records are pair/operator diagnostics
  and consumers.
- They are relevant after source realization exists because they consume raw
  PQS plans/product units and pair factor data.
- They are not source-realization authority and should not define support
  windows, support ordering, or source-plan shape.

Existing source-plan consumer shape:

- `pqs_multilayer_shell_source_plan` exposes:
  - `core_support_indices`
  - `core_support_states`
  - `shell_support_indices`
  - `shell_support_states`
  - `shell_final_coefficients`
  - `core_box`, `outer_box`, `bond_axis`, metrics, and layer records.
- Support one-body/density helpers use:
  `vcat(plan.core_support_states, plan.shell_support_states)`.
- The existing support ordering is therefore
  `:core_support_states_then_shell_support_states`.

Be2 support-order implication:

- If the Be2/PQS source-box route were forced into the current
  `core_then_shell` source-plan consumer shape, the only honest geometric
  analogue would be:
  - core/body first: `product` slab support states;
  - shell second: `pqs_left` support states followed by `pqs_right` support
    states.
- That conflicts with the current retained route order and ranges:
  `(:pqs_left, :pqs_right, :product)`.
- A real materializer would therefore need an explicit ordering/permutation
  contract. It should not simply return `object_kind =
  :pqs_multilayer_shell_source_plan` until the support order and coefficient
  block shape satisfy the consumer exactly.

Recommended next coding pass:

- Add a small private source-box support-window extractor, not a real
  materializer yet.
- Suggested names:
  - `_PQSDiatomicCompleteCoreShellSupportWindowPayload`
  - `_pqs_source_box_route_driver_diatomic_complete_core_shell_support_window_payload`
- Inputs:
  - `parent` fields: system classification, bond axis, parent axis counts or
    axis-bundle dimensions, parent axis bundle availability;
  - `route_skeleton.source_boxes`;
  - `route_skeleton.source_dimensions`;
  - `route_skeleton.retained_units`;
  - `route_skeleton.ranges` / retained dimension.
- Output fields/blockers:
  - `status`, `blocker`;
  - `parent_dims`;
  - `source_box_windows` as 3-tuples of ranges;
  - `source_mode_dims`;
  - `retained_order = (:pqs_left, :pqs_right, :product)`;
  - `candidate_core_then_shell_support_order = (:product, :pqs_left, :pqs_right)`;
  - `retained_to_support_order_permutation_required = true`;
  - `support_counts`;
  - optional metadata-only support index/state counts or hashes, if cheap and
    derived only from ranges/dims;
  - `missing_objects = (:raw_product_box_plan_objects,
    :pqs_axis_local_coefficients,
    :diatomic_complete_core_shell_source_plan_materializer)`.
- Wire this payload into the pass-138 source-plan payload as a structured
  readiness subobject/summary. Keep `source_plan === nothing`.

Why this pass, not the alternatives:

- Not a raw product-box probe fingerprint: the probe is intentionally
  metadata-only and should remain probe/debug.
- Not a real source-plan materializer yet: the support ordering/permutation
  issue and coefficient source boundary are still unresolved.
- Not another generic readiness payload alone: pass 138 already introduced the
  source-plan blocker. The next useful route-owned fact is the concrete
  source-window/order boundary.

Focused validation for that next coding pass:

- `julia --project=. test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl`

Git status:

```text
## main...origin/main
```

Deletion/shrinkage forecast:

- deleted: none recommended for the next pass.
- simplified: source-window/order facts would move out of ad hoc skeleton
  interpretation and into one compact route-owned payload.
- quarantined: raw product-box probe results, shell/support-row contraction,
  final-basis/H1/H1-J/Ham materialization, density-density materialization,
  RHF/SCF/Fock, WL payloads, public APIs, exports, artifacts, hfdmrg, and CR2.
- not deleted because: existing skeleton summaries, Ham readiness, and probe
  fields are still the active compatibility/report surfaces.
- exact remaining caller/blocker: no private diatomic complete core/shell
  source materializer currently converts route source boxes plus parent axis
  bundles into raw product-box plans, coefficient blocks, support-order
  permutation/provenance, and an honest source-plan consumer shape.

-- repo-doer@macmini
