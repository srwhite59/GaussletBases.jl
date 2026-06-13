Pass 142 response - add private diatomic raw-box route payload

Implemented and committed:

- commit: `8a96e308 Add diatomic raw box route payload`

Files changed:

- `src/pqs_source_box_route_driver_helpers.jl`
- `test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl`

Payload/helper names:

- private payload type: `_PQSDiatomicRawBoxRoutePayload`
- payload helper:
  `_pqs_source_box_route_driver_diatomic_raw_box_route_payload`
- producer used:
  `CartesianContractedParentMetrics._pqs_pqs_product_raw_box_route_producer`
- axis-metric helper:
  `_pqs_source_box_route_driver_diatomic_axis_metrics`
- compact summary helpers:
  - `_pqs_source_box_route_driver_raw_product_box_plan_summary`
  - `_pqs_source_box_route_driver_raw_pqs_plan_summary`
  - `_pqs_source_box_route_driver_product_unit_summary`
  - `_pqs_source_box_route_driver_raw_box_pair_inventory_summary`
  - `_pqs_source_box_route_driver_raw_box_descriptor_summary`
- assembly field: `diatomic_raw_box_route_payload`

Default Be2/PQS payload:

- `status = :blocked_diatomic_raw_box_route_payload`
- `blocker = :missing_parent_axis_bundle_object`
- `producer = nothing`
- `producer_status = :not_available`
- support-window status remains
  `:available_diatomic_complete_core_shell_support_windows`
- missing objects include:
  - `:parent_axis_bundle_object`
  - `:raw_product_box_plan_objects`
  - `:pqs_axis_local_coefficients`
  - `:product_doside_unit`
  - `:pair_inventory`

Probe-enabled Be2/PQS payload:

- `status = :available_diatomic_raw_box_route_payload`
- `blocker = nothing`
- `producer_status = :private_shadow_only`
- descriptor:
  - object kind `:pqs_pqs_product_safe_term_route_descriptor`
  - roles `(:pqs_left, :pqs_right, :product)`
  - retained dimension `221`
  - expected pair count `6`
- raw product-box plans:
  - left/right object kind `:cartesian_raw_product_box_plan_3d`
  - left/right source mode dims `(5, 5, 5)`
  - left/right axis-local coefficient shapes `((5, 5), (5, 5), (5, 5))`
- raw PQS summaries:
  - left/right representation `:orthogonal_raw_product_box`
  - left/right boundary selected count `98`
- product/doside unit:
  - kind `:product_doside`
  - support count `25`
  - support state count `25`
  - coefficient matrix shape `(25, 25)`
- pair inventory:
  - every pair uses source-box algorithmic policy
  - source-box algorithmic pair count `6`

Source-plan/Ham behavior:

- The pass-138 source-plan payload now summarizes
  `raw_box_route_payload_status`.
- Default source-plan blocker remains `:missing_parent_axis_bundle_object`.
- Probe-enabled source-plan blocker remains
  `:missing_diatomic_complete_core_shell_source_realization_contract`.
- Existing Ham readiness blocker remains
  `:missing_diatomic_complete_core_shell_source_plan_producer`.

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

Validation:

- `julia --project=. test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl`
  passed.
  - precompile dominated with about 56 s.
  - focused test summaries:
    - `Be2 PQS Ham payload readiness fingerprint | Pass 116 Total 116 Time 12.7s`
    - `Be2 PQS probe-enabled Ham readiness fingerprint | Pass 115 Total 115 Time 32.0s`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  passed with `load ok`.
- `git diff --check` passed.

Git status:

```text
## main...origin/main [ahead 1]
```

Deletion/shrinkage report:

- deleted: removed the direct/ad hoc raw-box producer fingerprint assertions
  from the focused test.
- simplified: tests now prefer the route-owned private
  `diatomic_raw_box_route_payload`; pass-141 direct producer checks were reduced
  into compact payload summary assertions.
- quarantined: raw product-box probe output, source-plan/final-basis/H1/H1-J/Ham
  materialization, RHF/SCF/Fock, WL payloads, public APIs, exports, artifacts,
  hfdmrg, and CR2 remain outside this payload.
- not deleted because: support-window, source-plan, and Ham readiness payloads
  remain active compatibility blockers and private route-driver surfaces.
- exact remaining caller/blocker: no source-realization materializer yet turns
  the raw-box payload plus support-order/permutation facts into an honest
  `:pqs_multilayer_shell_source_plan` consumer shape.

-- repo-doer@macmini
