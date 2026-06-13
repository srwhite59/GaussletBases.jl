Pass 141 response - fingerprint Be2 raw-box route producer

Implemented and committed:

- commit: `9b801052 Fingerprint Be2 raw box producer`

Files changed:

- `test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl`

Exact producer/helper called:

- `_PQS_BE2_CCPM._pqs_pqs_product_raw_box_route_producer`
- `_PQS_BE2_CCPM._pqs_raw_product_box_plan_view`
- local test-only helper `_pqs_be2_axis_metrics`
- existing focused fixture ingredients from
  `_pqs_be2_ham_payload_fingerprint_components`

Fingerprint result:

- producer object kind: `:pqs_pqs_product_raw_box_route_producer`
- producer status: `:private_shadow_only`
- descriptor object kind: `:pqs_pqs_product_safe_term_route_descriptor`
- descriptor roles: `(:pqs_left, :pqs_right, :product)`
- descriptor retained dimension: `221`
- descriptor expected pair count: `6`

Raw product-box plans:

- available for `:pqs_left` and `:pqs_right`
- both have `object_kind = :cartesian_raw_product_box_plan_3d`
- both have `source_mode_dims = (5, 5, 5)`
- both carry three `axis_local_coefficients` matrices with shape `(5, 5)`

Raw PQS/source coefficient objects:

- available through `producer.raw_pqs_plans.pqs_left` and
  `producer.raw_pqs_plans.pqs_right`
- `_pqs_raw_product_box_plan_view(...)` reports
  `representation = :orthogonal_raw_product_box`
- both boundary selectors retain `98` PQS boundary modes
- both carry three axis-local coefficient arrays

Product/doside unit:

- available as `producer.product_unit`
- `kind = :product_doside`
- support index count: `25`
- support state count: `25`
- coefficient matrix shape: `(25, 25)`

Pair inventory:

- available as `producer.all_pairs_inventory`
- `every_pair_uses_source_box_algorithmic_policy = true`
- `source_box_algorithmic_pair_count = 6`

Private/shadow/probe-only status:

- `producer.diagnostics.private_shadow_only = true`
- raw product-box plans and retained rules are built by this private producer,
  but dense raw pair matrices are not materialized by it.
- validation-only dense raw pair matrix flag remains true.
- no shell projection, Lowdin cleanup, support-local PQS oracle,
  support-coefficient matrix, retained PQS weights, IDA weight division,
  packet adoption, fixed-block routing, QWHamiltonian consumption, or public
  default consumption is claimed.

Validation:

- `julia --project=. test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl`
  passed.
  - `Be2 PQS Ham payload readiness fingerprint | Pass 94 Total 94 Time 11.7s`
  - `Be2 PQS probe-enabled Ham readiness fingerprint | Pass 73 Total 73 Time 17.8s`
  - `Be2 PQS raw-box producer fingerprint | Pass 46 Total 46 Time 12.4s`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  passed with `load ok`.
- `git diff --check` passed.

Git status:

```text
## main...origin/main [ahead 1]
```

Deletion/shrinkage report:

- deleted: none; this was a test-only fingerprint pass.
- simplified: the focused Be2 fixture helper now exposes components for
  producer-level fingerprints while preserving the existing assembly helper.
- quarantined: raw-box producer remains private/shadow-only and is not wired
  into `cartesian_assembly(...)`; source-plan/final-basis/H1/H1-J/Ham,
  RHF/SCF/Fock, WL payloads, public APIs, exports, artifacts, hfdmrg, and CR2
  remain out of scope.
- not deleted because: the producer is still evidence for constructible raw-box
  route objects, not route authority.
- exact remaining caller/blocker: no route-owned diatomic source-realization
  payload/materializer yet wraps these raw product-box plans, PQS coefficient
  objects, product unit, pair inventory, and support-order permutation
  provenance into an honest source-plan consumer shape.

-- repo-doer@macmini
