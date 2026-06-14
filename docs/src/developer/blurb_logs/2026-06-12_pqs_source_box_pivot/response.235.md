# Pass 235 response - delete route-shadow density-density fixture pressure

Implemented the cleanup pass. This deletes the old PQS/PQS/product
route-shaped density-density fixture pressure while preserving compact live
math convention checks.

## Files changed

- `src/cartesian_contracted_parent_metrics/current_route_metadata_export.jl`
- `test/nested/pqs_projected_q_shell_local_layer_integration_runtests.jl`

No commits or pushes.

## Deleted or shrunk

Deleted from `current_route_metadata_export.jl`:

- `_PQS_PQS_PRODUCT_DENSITY_DENSITY_ROUTE_KINDS`
- `_source_box_axis_pair_terms_symmetric(...)`
- `_pqs_pqs_product_density_density_route_ranges(...)`
- `_pqs_pqs_product_density_density_all_pairs_inventory(...)`
- `_pqs_pqs_product_density_density_pair_block(...)`
- `_pqs_pqs_product_route_shaped_density_density_consumer(...)`
- `_pqs_pqs_product_raw_box_density_density_route_producer(...)`
- `_pqs_pqs_product_raw_box_density_density_route_producer_from_ida_provenance(...)`
- `_pqs_route_parent_dims(...)`
- `_product_doside_parent_coefficient_matrix(...)`
- `_pqs_pqs_product_route_parent_coefficient_matrix(...)`
- `_pqs_pqs_product_dense_parent_ida_authority_comparison(...)`

Shrunk from `pqs_projected_q_shell_local_layer_integration_runtests.jl`:

- `_density_density_route_factor_fixture(...)`
- `_check_pqs_pqs_product_route_shaped_density_density_consumer(...)`
- `_check_pqs_pqs_product_raw_box_density_density_route_producer(...)`
- the two call sites that exercised the deleted route-shaped density-density
  producer/consumer path.

## Caller-search evidence

The deleted density-density route-shaped names now have no callers in
`src`, `test`, or `bin`:

```sh
rg -n '_pqs_pqs_product_raw_box_density_density_route_producer\b|_pqs_pqs_product_raw_box_density_density_route_producer_from_ida_provenance\b|_pqs_pqs_product_route_shaped_density_density_consumer\b|_pqs_pqs_product_density_density_all_pairs_inventory\b|_pqs_pqs_product_density_density_pair_block\b|_pqs_pqs_product_density_density_route_ranges\b|_pqs_pqs_product_dense_parent_ida_authority_comparison\b|_PQS_PQS_PRODUCT_DENSITY_DENSITY_ROUTE_KINDS|_source_box_axis_pair_terms_symmetric\b|_pqs_pqs_product_route_parent_coefficient_matrix\b|_product_doside_parent_coefficient_matrix\b|_pqs_route_parent_dims\b' src test bin
```

Result: no matches.

I did not delete `_pqs_pqs_product_raw_box_route_producer(...)`; caller search
shows a non-test source caller in
`src/pqs_source_box_diatomic_complete_core_shell.jl`.

I did not delete `_pqs_parent_coefficient_matrix_from_raw_plan(...)`; caller
search shows a non-test source caller in
`src/pqs_source_box_diatomic_complete_core_shell.jl`.

## Preserved math/scientific checks

The shrunk integration fixture still preserves:

- boundary retained count convention:
  `rectangular.diagnostics.boundary_support_count == 130` and
  `x_axis.diagnostics.boundary_support_count == 130`;
- density-normalized versus raw-weighted pair convention for product/product,
  PQS/product, and PQS/PQS lower-level source-box blocks;
- nuclear attraction charge/sign convention, including
  `sign_charge_scale == -nuclear_charge` and
  `nuclear_charge_sign_applied`.

No independent-H2 support-plan implementation, H1, H1-J, RHF, supplements, CR2,
export, public API, fake-PQS guard field, or fake-PQS golden regression was
changed.

## Validation

Package load:

```sh
julia --project=. -e 't = @elapsed begin using GaussletBases; println("load ok") end; println("elapsed_s=", t)'
```

Result:

```text
load ok
elapsed_s=58.613955625
```

The load command precompiled `GaussletBases` after the source deletion.

Focused integration check:

```sh
julia --project=. -e 'using Test, LinearAlgebra, SparseArrays, GaussletBases; t = @elapsed include("test/nested/pqs_projected_q_shell_local_layer_integration_runtests.jl"); println("elapsed_s=", t)'
```

Result:

```text
Test Summary:                                  | Pass  Total     Time
Cartesian nested projected q-shell local layer | 3820   3820  6m04.6s
elapsed_s=845.430618416
```

I first tried the same include without `using Test`; it failed immediately with
`UndefVarError: @testset not defined`. The successful command above used the
fixture's required harness imports.

Whitespace:

```sh
git diff --check
```

Result: passed.

## Line budget

```text
git diff --numstat -- src test bin
19      1243    src/cartesian_contracted_parent_metrics/current_route_metadata_export.jl
0       344     test/nested/pqs_projected_q_shell_local_layer_integration_runtests.jl
```

Scoped total: 19 added / 1587 deleted, net -1568.

## Git status

Before writing this response:

```text
## main...origin/main
 M src/cartesian_contracted_parent_metrics/current_route_metadata_export.jl
 M test/nested/pqs_projected_q_shell_local_layer_integration_runtests.jl
```

## Carrying-cost accounting

deleted:
- old PQS/PQS/product route-shaped density-density producer/consumer source
  helpers;
- old IDA-provenance route-shaped density-density adapter and dense-parent
  authority comparison;
- slow integration fixture pressure that existed mainly to assert private
  route-shadow metadata and helper vocabulary.

simplified:
- the remaining integration fixture now keeps lower-level operator convention
  checks without preserving the complete route-shaped density-density wrapper.

quarantined:
- fake-PQS golden regression remains untouched;
- independent H2 PQS remains support-plan/readiness only.

not deleted because:
- `_pqs_pqs_product_raw_box_route_producer(...)` still has a non-test source
  caller in the route driver path;
- `_pqs_parent_coefficient_matrix_from_raw_plan(...)` still has a non-test
  source caller in the route driver path;
- lower-level product/PQS/PQS density-density and nuclear-attraction helpers
  still carry the compact math conventions this pass was told to preserve.

exact remaining caller/blocker:
- route skeleton helper-name/report pressure remains in active route-driver and
  pair-stage surfaces; caller search includes `src/pqs_source_box_route_driver_*`
  and focused route-skeleton tests, so it was not used as deletion fuel in this
  pass.

-- repo-doer@macmini
