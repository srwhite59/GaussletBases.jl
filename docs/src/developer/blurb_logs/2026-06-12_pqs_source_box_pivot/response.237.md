# Pass 237 response - independent H2 PQS retained-rule readiness plan

Implemented the metadata/readiness-only independent H2 PQS retained-rule plan.

## Files changed

- `src/pqs_source_box_diatomic_complete_core_shell.jl`
- `src/pqs_source_box_route_driver_reporting.jl`
- `src/pqs_source_box_route_driver_skeletons.jl`
- `test/driver_inputs/h2_pqs_q5_independent_source_box_r4.jl`
- `test/nested/runtests.jl`
- `test/nested/integration_runtests.jl`

Deleted:

- `test/nested/pqs_source_box_route_skeleton_runtests.jl`
- `test/nested/pqs_route_axis_count_selection_runtests.jl`

## Retained-rule plan

Added private helper:

```text
_pqs_source_box_route_driver_independent_h2_retained_rule_plan(...)
```

Plan status/counts:

```text
status:                   :available_independent_pqs_retained_rule_readiness_plan
support_counts:           (275, 578, 362)
retained_counts:          (275, 98, 98)
expected_final_dimension: 471
authority:                :independent_pqs_route_owned_retained_rule_readiness
retained_counts_source:   :route_owned_retained_rule_plan
```

Per-unit retained-rule authority:

```text
:atom_contact_core => :direct_source_modes
:shared_shell_1    => :pqs_boundary_comx_product_modes
:shared_shell_2    => :pqs_boundary_comx_product_modes
```

The retained-rule plan composes the generated support-region plan plus the
existing q=5 PQS boundary product-mode retained rule. It does not materialize
source coefficients or transforms.

## Reporting/artifact behavior

- Independent H2 target metadata now reports retained counts `(275, 98, 98)`.
- Target expected final dimension/readiness is `471`.
- `fake_pqs/enabled = false` remains asserted by validation.
- `route/source_backed_fixed_source_oracle_used = false` remains asserted by validation.
- `physics/endpoint_ready = false` remains asserted by validation.
- The remaining endpoint blocker is now:

```text
:missing_independent_pqs_physical_source_plan_materializer
```

I also removed stale independent-route references to the old atom-core and
shared-shell retained-rule blockers from the skeleton/source-plan metadata so
the report no longer presents those as active blockers after this pass.

## Forbidden surfaces avoided

No source coefficient matrices, final basis, H1, H1-J, RHF, supplements, CR2,
exports, public API, fake-PQS/WL coefficient matrices, or fixed-source retained
transforms were added.

## Validation

Package load:

```text
julia --project=. -e 't = @elapsed begin using GaussletBases; println("load ok") end; println("elapsed_s=", t)'
```

Result:

```text
load ok
elapsed_s=58.647891042
```

Focused independent-input artifact/readiness check:

```text
julia --project=. -e 'using JLD2; ... include("bin/cartesian_ham_builder.jl") ...'
```

Assertions checked:

```text
fake_pqs/enabled == false
route/source_backed_fixed_source_oracle_used == false
physics/endpoint_ready == false
physics/endpoint_blocker == :missing_independent_pqs_physical_source_plan_materializer
target/support_counts == (275, 578, 362)
target/retained_counts == (275, 98, 98)
target/retained_order == (:atom_contact_core, :shared_shell_1, :shared_shell_2)
target/expected_final_dimension == 471
target/retained_atom_core_interiors == true
target/support_plan_status == :available_independent_pqs_support_region_plan
target/retained_rule_plan_status == :available_independent_pqs_retained_rule_readiness_plan
target/retained_rule_plan_authority == :independent_pqs_route_owned_retained_rule_readiness
target/retained_counts_generated == true
target/retained_counts_source == :route_owned_retained_rule_plan
target/retained_rule_expected_final_dimension == 471
target/retained_rule_authorities == (:direct_source_modes, :pqs_boundary_comx_product_modes, :pqs_boundary_comx_product_modes)
```

Result:

```text
retained_rule_plan_status=available_independent_pqs_retained_rule_readiness_plan
retained_counts=(275, 98, 98)
expected_final_dimension=471
physics_endpoint_blocker=missing_independent_pqs_physical_source_plan_materializer
elapsed_s=82.246477625
```

Whitespace:

```text
git diff --check
```

passed.

Stale blocker search:

```text
rg -n "missing_independent_pqs_atom_contact_core_retained_rule|independent_pqs_atom_contact_core_retained_rule|independent_pqs_shared_shell_2_retained_rule" src test bin
```

returned no matches.

## Scoped line budget

For `src + test + bin`:

```text
173 added / 248 deleted, net -75
```

`git diff --numstat -- src test bin`:

```text
134	9	src/pqs_source_box_diatomic_complete_core_shell.jl
27	1	src/pqs_source_box_route_driver_reporting.jl
10	15	src/pqs_source_box_route_driver_skeletons.jl
2	2	test/driver_inputs/h2_pqs_q5_independent_source_box_r4.jl
0	1	test/nested/integration_runtests.jl
0	130	test/nested/pqs_route_axis_count_selection_runtests.jl
0	89	test/nested/pqs_source_box_route_skeleton_runtests.jl
0	1	test/nested/runtests.jl
```

## Git status

```text
## main...origin/main
 M src/pqs_source_box_diatomic_complete_core_shell.jl
 M src/pqs_source_box_route_driver_reporting.jl
 M src/pqs_source_box_route_driver_skeletons.jl
 M test/driver_inputs/h2_pqs_q5_independent_source_box_r4.jl
 M test/nested/integration_runtests.jl
 D test/nested/pqs_route_axis_count_selection_runtests.jl
 D test/nested/pqs_source_box_route_skeleton_runtests.jl
 M test/nested/runtests.jl
```

## Deletion/shrinkage accounting

deleted:
- `test/nested/pqs_source_box_route_skeleton_runtests.jl`
- `test/nested/pqs_route_axis_count_selection_runtests.jl`
- their two nested-suite includes

simplified:
- independent H2 route readiness now stops carrying old fake/WL retained-rule blockers as active independent-PQS blockers;
- remaining independent route blocker is the source-plan materializer;
- visible input/skeleton/report metadata now align with retained-rule readiness target `471`.

quarantined:
- fake-PQS WL/QW reproduction path remains separate and unchanged;
- this pass remains metadata/readiness-only.

not deleted because:
- fake-PQS H2 463 remains a golden regression until the independent PQS endpoint exists;
- active route-driver/reporting surfaces still need to carry readiness metadata for the next source-plan materialization step;
- no final-basis/H1/H1-J/RHF/supplement code was in scope.

exact remaining caller/blocker:
- `:missing_independent_pqs_physical_source_plan_materializer`

-- repo-doer@macmini
