Pass 140 response - add diatomic PQS support-window payload

Implemented and committed:

- commit: `0cebcbf8 Add diatomic PQS support window payload`

Files changed:

- `src/pqs_source_box_route_driver_helpers.jl`
- `test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl`

Payload/helper names:

- private payload type:
  `_PQSDiatomicCompleteCoreShellSupportWindowPayload`
- support-window helper:
  `_pqs_source_box_route_driver_diatomic_complete_core_shell_support_window_payload`
- small local helpers:
  - `_pqs_source_box_route_driver_axis_counts_tuple`
  - `_pqs_source_box_route_driver_source_box_window`
  - `_pqs_source_box_route_driver_window_support_count`
- assembly field:
  `diatomic_complete_core_shell_support_window_payload`

Wiring:

- `cartesian_assembly(...)` now builds
  `diatomic_complete_core_shell_support_window_payload`.
- The pass-138 source-plan payload now carries this object as
  `support_window_payload` and records `support_window_payload_status` in its
  summary/metadata.
- Existing source-plan and Ham blockers remain behavior-preserving.

Default Be2/PQS support-window payload:

- `status = :available_diatomic_complete_core_shell_support_windows`
- `blocker = nothing`
- `parent_dims = (9, 7, 9)`
- `parent_axis_bundle_object_available = false`
- `source_box_windows`:
  - `pqs_left = (1:5, 1:5, 1:5)`
  - `product = (1:5, 1:5, 5:5)`
  - `pqs_right = (1:5, 1:5, 5:9)`
- `source_mode_dims`:
  - `pqs_left = (5, 5, 5)`
  - `product = (5, 5, 1)`
  - `pqs_right = (5, 5, 5)`
- `support_counts = (pqs_left = 125, product = 25, pqs_right = 125)`

Probe-enabled Be2/PQS support-window payload:

- `status = :available_diatomic_complete_core_shell_support_windows`
- `blocker = nothing`
- `parent_axis_bundle_object_available = true`
- `parent_dims` is available from the constructed parent-axis path; by the
  current Be2 parent-axis convention this is `(31, 17, 17)`.
- The derived x/y/z windows remain compact range tuples:
  - `pqs_left = (1:5, 1:5, 1:5)`
  - `product = (1:5, 1:5, 9:9)`
  - `pqs_right = (1:5, 1:5, 13:17)`
- `support_counts = (pqs_left = 125, product = 25, pqs_right = 125)`

Order/permutation facts:

- `retained_order = (:pqs_left, :pqs_right, :product)`
- `candidate_core_then_shell_support_order = (:product, :pqs_left, :pqs_right)`
- `retained_to_support_order_permutation_required = true`

Missing objects still recorded:

- `:raw_product_box_plan_objects`
- `:pqs_axis_local_coefficients`
- `:diatomic_complete_core_shell_source_plan_materializer`

No materialization claims:

- `support_states_materialized = false`
- `raw_product_box_plans_materialized = false`
- `source_coefficients_materialized = false`
- `source_plan_materialized = false`
- `final_basis_materialized = false`
- `h1_materialized = false`
- `h1_j_materialized = false`
- `ham_payload_materialized = false`
- `route_driver_public_surface = false`
- `exports_materialized = false`
- `artifacts_materialized = false`
- metadata keeps `raw_product_box_probe_authority = false`

Validation:

- `julia --project=. test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl`
  passed.
  - precompile dominated with about 56 s.
  - focused test summaries:
    - `Be2 PQS Ham payload readiness fingerprint | Pass 94 Total 94 Time 12.3s`
    - `Be2 PQS probe-enabled Ham readiness fingerprint | Pass 73 Total 73 Time 18.6s`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  passed with `load ok`.
- `git diff --check` passed.

Git status:

```text
## main...origin/main [ahead 1]
```

Deletion/shrinkage report:

- deleted: none in this boundary-creation pass.
- simplified: source-window/order facts are now one private payload instead of
  implicit interpretation of route skeleton fields.
- quarantined: raw product-box probe outputs, raw product-box plan objects,
  support states, source coefficients, source-plan/final-basis/H1/H1-J/Ham
  materialization, RHF/SCF/Fock, WL payloads, public APIs, exports, artifacts,
  hfdmrg, and CR2 remain outside this payload.
- not deleted because: the existing source-plan and Ham readiness payloads are
  still the active compatibility blockers and report-facing private surfaces.
- exact remaining caller/blocker: no private diatomic complete core/shell source
  materializer currently converts these windows plus parent axis bundles into
  raw product-box plans, coefficient blocks, support-order permutation
  provenance, and an honest source-plan consumer shape.

-- repo-doer@macmini
