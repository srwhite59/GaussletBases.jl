Pass 204 review - accepted

Reviewed the split of `src/cartesian_contracted_parent_metrics/source_box_route_shadow.jl`.

Result:
- Accepted.
- The mixed 8601-line implementation file was deleted.
- `src/CartesianContractedParentMetrics.jl` now includes three explicit private
  files for the former source-box route-shadow code:
  - `product_staged_metric_fallbacks.jl`
  - `source_box_pair_shadow.jl`
  - `legacy_source_box_fixtures.jl`

What changed:
- No public API changes.
- No intended numerical behavior changes.
- No H1/H1-J/RHF/PQS multilayer physics changes.
- No new tests or public route concepts.

Resulting line counts checked by manager:
```text
product_staged_metric_fallbacks.jl  1305
source_box_pair_shadow.jl           3069
legacy_source_box_fixtures.jl       4224
CartesianContractedParentMetrics.jl   68
```

Validation run by manager:
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - passed.
- `julia --project=. -e 'using Test, LinearAlgebra, SparseArrays, GaussletBases; include("test/nested/cartesian_contracted_parent_metric_packet_runtests.jl")'`
  - passed, 29/29.
- `julia --project=. test/nested/pqs_component_route_report_adapter_runtests.jl`
  - passed, 605/605.
- `git diff --check`
  - passed.

Line budget:
- Deleted old `source_box_route_shadow.jl`: 8601 lines.
- Added new split files: 1305 + 3069 + 4224 = 8598 lines.
- Wrapper changed by +3/-4.
- Net source/test/bin line budget: -4.

Deletion-candidate audit:
- `_product_doside_source_box_shadow_blocks` looks like the clearest deletion
  candidate after a stale-test audit: no live source caller outside the legacy
  fixture file was reported, with direct pressure from
  `pqs_projected_q_shell_local_layer_integration_runtests.jl`.
- The rest of the audited product/doside source-box family still has callers in
  `current_route_metadata_export.jl`, component report expectations, route
  skeleton tests, high-order source-construction integration tests, or internal
  fixture paths.

Deletion/shrinkage assessment:
- deleted: mixed `source_box_route_shadow.jl` implementation file.
- simplified: the former junk-drawer include is now split by purpose:
  metric fallback, PQS source-box pair-plan prototype, and legacy fixtures.
- quarantined: product/doside fixture/oracle code and contact-cap/
  outer-mismatch/atom-box route shadows now live in
  `legacy_source_box_fixtures.jl`.
- not deleted because: product/doside helpers still have live source/report/test
  callers; deletion needs a dedicated stale-caller audit.
- exact remaining caller/blocker: `current_route_metadata_export.jl`,
  `pqs_projected_q_shell_local_layer_integration_runtests.jl`,
  `pqs_component_route_report_adapter_runtests.jl`,
  `pqs_source_box_route_skeleton_runtests.jl`, and high-order
  source-construction integration tests still preserve product/doside
  source-box fixture vocabulary.

Manager note:
- This satisfies the requested first simplification of the next hotspot.
- The natural next cleanup, if requested, is not another split. It is a
  stale-caller audit/deletion pass for the product/doside fixture family,
  starting with `_product_doside_source_box_shadow_blocks`.

-- repo-manager@macmini
