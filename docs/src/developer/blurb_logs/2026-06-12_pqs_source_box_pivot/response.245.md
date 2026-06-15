Pass 245 complete.

Files changed:
- `test/nested/cartesian_unit_pairs_contract_runtests.jl`
- `test/nested/cartesian_pair_operator_plans_contract_runtests.jl`

Exact lines added/deleted:
- `git diff --numstat -- src test bin`
  - `0  63  test/nested/cartesian_unit_pairs_contract_runtests.jl`
  - `0  37  test/nested/cartesian_pair_operator_plans_contract_runtests.jl`
- Total `src + test + bin`: `0` added / `100` deleted, net `-100`.

Helper/test blocks deleted or collapsed:
- Deleted `_unit_pairs_blocked_retained_plan(...)`.
- Deleted the synthetic `"CartesianUnitPairs missing RouteCore sidecars"` testset.
- Collapsed pair-operator metadata-only/test-vocabulary assertions by deleting:
  - compatibility no-transform-contract constructor mirror assertions;
  - `route_core_sidecars = false` not-requested inventory assertion cloud;
  - `:metadata_only_not_materialized` materialization count assertion;
  - missing-transform blocker count assertions for `:missing_left_transform_contract` and `:missing_right_transform_contract`;
  - duplicate-transform blocker count and duplicate-key metadata assertions.

Active contract preserved:
- Unit-pair tests still cover unavailable summary plus active metadata inventory:
  - retained unit count;
  - pair count and upper-triangular pair order;
  - unit pair index table order;
  - route-core pair inventory availability/count for the normal metadata path;
  - pair family counts;
  - non-materialized summary flags.
- Pair-operator tests still cover:
  - unavailable summary;
  - active pair-operator plan shape/counts;
  - source/operator path counts;
  - transform and realization path matching;
  - final block path summaries;
  - compact missing/duplicate transform-contract blocked summaries.

Validation and timings:
- `julia --project=. test/nested/cartesian_unit_pairs_contract_runtests.jl`
  - passed
  - Julia-reported testsets:
    - `CartesianUnitPairs unavailable summary`: 18 passed, 0.0s
    - `CartesianUnitPairs metadata inventory`: 27 passed, 0.7s
- `julia --project=. test/nested/cartesian_pair_operator_plans_contract_runtests.jl`
  - passed
  - Julia-reported testsets:
    - `CartesianPairOperatorPlans unavailable summary`: 22 passed, 0.0s
    - `CartesianPairOperatorPlans metadata inventory`: 23 passed, 2.1s
    - `CartesianPairOperatorPlans transform-contract authority`: 8 passed, 0.2s
    - `CartesianPairOperatorPlans transform-contract blockers`: 8 passed, 0.5s
- `git diff --check`
  - passed

Post-edit stale-symbol audit:
- `rg` found no remaining occurrences in the two edited files for:
  - `_unit_pairs_blocked_retained_plan`
  - `missing RouteCore sidecars`
  - `blocked_missing_route_core_final_units`
  - `synthetic_missing_sidecar`
  - `metadata_only_not_materialized`
  - `duplicate_transform_contract_unit_keys`
  - `missing_left_transform_contract`
  - `missing_right_transform_contract`
  - `duplicate_retained_unit_transform_contract`

Git status:
```text
## main...origin/main
 M test/nested/cartesian_pair_operator_plans_contract_runtests.jl
 M test/nested/cartesian_unit_pairs_contract_runtests.jl
```

Deletion/shrinkage accounting:
- deleted: synthetic unit-pair missing-sidecar helper/testset and pair-operator count-only blocker/materialization assertions.
- simplified: pair-operator blocked-path tests now keep compact status/count/final-block readiness checks without preserving private blocker-count vocabulary.
- quarantined: none.
- not deleted because: active happy-path unit-pair shape/order/family checks and pair-operator shape/path/transform matching checks are still live module contracts.
- exact remaining blocker/caller: `_pair_ops_count(...)` remains because active path/final-block summary assertions still use it; no remaining caller for `_unit_pairs_blocked_retained_plan(...)`.

-- repo-doer@macmini
