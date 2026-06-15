# Pass 255 response - shrink atom-growth report-stage RouteCore mirrors

Cleanup-only pass completed. I edited only:

```text
test/nested/cartesian_report_stage_low_order_policy_runtests.jl
```

No source, H2/PQS route behavior, RHF, H1-J, supplements, CR2/export, public
API, driver input, or manifest files were changed.

## Deleted assertion families

Deleted stale atom-growth report-stage RouteCore mirror assertions:

```text
low_order_route_core_pair_family_counts == summary.route_core_pair_family_counts
low_order_route_core_pair_operator_preflight.route_core_final_unit_count == 8
low_order_route_core_pair_operator_preflight.route_core_pair_count == 36
!low_order_route_core_pair_operator_preflight.operator_blocks_materialized
low_order_route_core_pair_operator_plan.route_core_final_unit_count == 8
low_order_route_core_pair_operator_plan.route_core_pair_count == 36
!low_order_route_core_pair_operator_plan.operator_blocks_materialized
!low_order_route_core_pair_operator_plan.hamiltonian_matrices_materialized
```

## Coverage kept

Compact atom-growth RouteCore readiness coverage remains:

```text
atom-growth route selected, not terminal, not legacy
materialization required/deferred and no Hamiltonian/operator/pair blocks
pair inventory known/available
pair count = 36
family-count assertion on the compact summary
RouteCore summary/inventory availability/status
RouteCore final-unit count = 8
RouteCore pair count = 36
pair-operator ready/status/blocker at summary/report level
preflight availability/status/blocker at report level
plan availability/status/blocker at report level
LW complete-shell enumeration checks
```

## Line budget

Scoped `src + test + bin`:

```text
0 insertions, 13 deletions
net -13
```

Deletion/shrinkage report:

```text
deleted:
- stale nested preflight/plan object field mirrors
- stale flat report-to-summary family-count mirror

simplified:
- atom-growth report-stage RouteCore coverage now keeps compact readiness facts
  without asserting private nested object shapes

quarantined:
- nested RouteCore preflight/plan object details remain internal test-only
  scaffolding, not route authority

not deleted because:
- compact atom-growth route selection, pair inventory, pair count, preflight,
  plan, and LW enumeration checks remain active

exact remaining caller/blocker:
- no remaining caller requires the deleted nested preflight/plan shape mirrors
```

## Validation

Commands run:

```text
julia --project=. -e 'Meta.parseall(read("test/nested/cartesian_report_stage_low_order_policy_runtests.jl", String)); println("parse ok")'
git diff --check
git status --short --branch
```

Results:

```text
parse smoke: passed
git diff --check: passed
```

I did not run the broad report-stage test; it is not required for this
stale-mirror cleanup and has been treated in this loop as a poor per-pass gate.

Final git status:

```text
## main...origin/main
 M test/nested/cartesian_report_stage_low_order_policy_runtests.jl
```

-- repo-doer@macmini
