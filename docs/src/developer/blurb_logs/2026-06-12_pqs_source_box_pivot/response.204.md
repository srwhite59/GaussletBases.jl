Pass 204 response - split source_box_route_shadow

Implementation:
- Deleted the mixed implementation file:
  - `src/cartesian_contracted_parent_metrics/source_box_route_shadow.jl`
- Created three behavior-preserving private include files:
  - `src/cartesian_contracted_parent_metrics/product_staged_metric_fallbacks.jl`
  - `src/cartesian_contracted_parent_metrics/source_box_pair_shadow.jl`
  - `src/cartesian_contracted_parent_metrics/legacy_source_box_fixtures.jl`
- Updated `src/CartesianContractedParentMetrics.jl` to include the three new files directly.
- Preserved private names and behavior. No public API, numerical, H1/H1-J/RHF, route, or physics changes.

Resulting include line counts:
```text
product_staged_metric_fallbacks.jl  1305
source_box_pair_shadow.jl           3070
legacy_source_box_fixtures.jl       4226
CartesianContractedParentMetrics.jl   68
```

Function-family placement:
- `product_staged_metric_fallbacks.jl`
  - staged axis utilities;
  - product-staged retained metric blocks;
  - PQS/product low-order support-local metric references;
  - fallback staged separable/kinetic/x2 metric block helpers;
  - support-local retained entries, staged unit entries, and resolved-payload metric fallback helpers;
  - linear-vector helpers used by metric packet construction.
- `source_box_pair_shadow.jl`
  - `_pqs_raw_product_box_*` source-mode, selector, structural, operator, and realization plans;
  - `_pqs_product_source_box_pair_plan`;
  - `_pqs_pqs_source_box_pair_plan`;
  - PQS/product and PQS/PQS source-box block-from-factor helpers;
  - PQS/product/PQS/PQS source-box density, local Gaussian, nuclear, and reference helpers tied to the pair-plan prototype path.
- `legacy_source_box_fixtures.jl`
  - product/doside source-box shadow/reference block family;
  - product/doside density-density, local Gaussian, raw-weighted density, and nuclear fixture blocks;
  - homonuclear raw-box geometry fixture producers;
  - PQS/PQS/product all-pairs inventory and route descriptor fixtures;
  - contact-cap, outer-mismatch, atom-box, and support-local oracle comparison helpers;
  - route-retained-unit fact audit helpers.

Scoped line budget:
- `git diff --numstat -- src test bin tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl` currently reports tracked-file changes only:
  - `3 4 src/CartesianContractedParentMetrics.jl`
  - `0 8601 src/cartesian_contracted_parent_metrics/source_box_route_shadow.jl`
- New include files are untracked at response time, so their line counts are accounted separately:
  - `1305 product_staged_metric_fallbacks.jl`
  - `3070 source_box_pair_shadow.jl`
  - `4226 legacy_source_box_fixtures.jl`
- source/test/bin scoped lines added: 8604
- source/test/bin scoped lines deleted: 8605
- net: -1

Deletion-candidate audit:
- `_product_doside_source_box_pair_plan`
  - Source callers: internal product/doside reference path.
  - Test/doc callers: legacy/source-box policy docs; direct behavior coverage is in `test/nested/pqs_projected_q_shell_local_layer_integration_runtests.jl`.
  - Status: not deletion-ready; feeds product/doside shadow/reference helpers.
- `_product_doside_source_box_reference_block`
  - Source callers: `legacy_source_box_fixtures.jl`, `source_box_pair_shadow.jl`, and `current_route_metadata_export.jl`.
  - Test callers: `pqs_projected_q_shell_local_layer_integration_runtests.jl`, high-order source-construction integration tests, and route skeleton/report expectations.
  - Status: not deletion-ready; still used by current-route metadata/export and fixture/oracle comparisons.
- `_product_doside_source_box_shadow_blocks`
  - Source callers: no live source caller outside legacy fixture code.
  - Test callers: `pqs_projected_q_shell_local_layer_integration_runtests.jl`.
  - Status: deletion candidate after a separate stale-test audit.
- `_product_doside_source_box_local_gaussian_sum_block`
  - Source callers: internal nuclear/local Gaussian fixture helpers.
  - Test callers: no direct test call found outside the integration fixture path.
  - Status: not standalone deletion-ready; can be audited with the product/doside nuclear fixture family.
- `_product_doside_source_box_density_density_interaction_block`
  - Source callers: `legacy_source_box_fixtures.jl` and `current_route_metadata_export.jl`.
  - Test callers: `pqs_projected_q_shell_local_layer_integration_runtests.jl`, `pqs_source_box_route_skeleton_runtests.jl`, and component report adapter expectations.
  - Status: not deletion-ready; route/report metadata still names and calls it.
- `_product_doside_source_box_nuclear_attraction_by_center`
  - Source callers: `current_route_metadata_export.jl`.
  - Test callers: `pqs_projected_q_shell_local_layer_integration_runtests.jl` and component report adapter expectations.
  - Status: not deletion-ready; still used by current-route component metadata.

Validation:
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - passed; printed `load ok` after precompilation.
- `julia --project=. -e 'using Test, LinearAlgebra, SparseArrays, GaussletBases; include("test/nested/cartesian_contracted_parent_metric_packet_runtests.jl")'`
  - passed, 29/29.
- `julia --project=. test/nested/pqs_component_route_report_adapter_runtests.jl`
  - passed, 605/605.
- Not run by default:
  - `julia --project=. test/nested/pqs_projected_q_shell_local_layer_integration_runtests.jl`
  - Reason: it is explicitly marked integration/slow and is the direct legacy fixture/pair-plan coverage file; the pass requested not running a >60s-only legacy test by default.
- `git diff --check`
  - passed.
- `git diff --cached --check`
  - passed; no staged diff at response time.
- `rg -n "source_box_route_shadow\\.jl" src test bin examples docs --glob '!docs/src/developer/blurb_logs/**'`
  - no matches, exit 1 expected.

Git status:
```text
## main...origin/main
 M src/CartesianContractedParentMetrics.jl
 D src/cartesian_contracted_parent_metrics/source_box_route_shadow.jl
?? src/cartesian_contracted_parent_metrics/legacy_source_box_fixtures.jl
?? src/cartesian_contracted_parent_metrics/product_staged_metric_fallbacks.jl
?? src/cartesian_contracted_parent_metrics/source_box_pair_shadow.jl
```

Deletion/shrinkage report:
- deleted: the mixed `source_box_route_shadow.jl` implementation file.
- simplified: `CartesianContractedParentMetrics.jl` now includes explicitly named metric fallback, source-box pair shadow, and legacy fixture files.
- quarantined: product/doside fixture/oracle code and contact-cap/outer-mismatch/atom-box route shadows are isolated in `legacy_source_box_fixtures.jl`.
- not deleted because: the audited product/doside helpers still have source callers in `current_route_metadata_export.jl` and live fixture/report tests; deleting them safely requires a separate stale-caller audit.
- exact remaining caller/blocker: `current_route_metadata_export.jl`, `test/nested/pqs_projected_q_shell_local_layer_integration_runtests.jl`, `test/nested/pqs_component_route_report_adapter_runtests.jl`, `test/nested/pqs_source_box_route_skeleton_runtests.jl`, and high-order source-construction integration tests still preserve product/doside source-box fixture vocabulary.

-- repo-doer@macmini
