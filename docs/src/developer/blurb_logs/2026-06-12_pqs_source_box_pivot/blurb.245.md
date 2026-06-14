# Pass 245 blurb - shrink synthetic RouteCore sidecar blocker tests

Role: repo-doer.

Read before editing:

- `AGENTS.md`
- `BlurbStyle.md`
- `docs/src/developer/pqs_manager_running_log.md`
- `docs/src/developer/old_flat_cartesian_retirement_audit_2026-06-14.md`
- `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.244.md`
- `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/review.244.md`

Task type: cleanup/test shrink.

Purpose:

Shrink synthetic missing-sidecar and metadata-only blocker vocabulary in active
unit-pair/pair-operator contract tests. These tests should protect current
module contracts such as pair shape, pair order, transform contract matching,
and compact summary readiness, not preserve every stale private
`route_core_sidecar_status` or metadata-only blocker field.

Deletion candidate:

```text
files:
  test/nested/cartesian_unit_pairs_contract_runtests.jl
  test/nested/cartesian_pair_operator_plans_contract_runtests.jl
```

Subagent/read-only audit summary:

```text
risk class:
  green/yellow.

expected line savings:
  about 60-120 test lines.

obsolete pressure:
  synthetic missing-sidecar blocker blocks and metadata-only materialization
  count clouds, including:
    :blocked_missing_route_core_final_units
    route_core_sidecar_status = :blocked
    route_core_sidecars = false
    metadata-only materialization count clouds

replacement/current authority:
  active unit-pair and pair-operator module summaries;
  RouteCore typed inventory when sidecars are requested.
```

Prioritize:

1. In `cartesian_unit_pairs_contract_runtests.jl`, delete the synthetic
   `"CartesianUnitPairs missing RouteCore sidecars"` testset if it only proves
   stale missing-sidecar blocker vocabulary.
2. Remove helper code used only by that testset, such as
   `_unit_pairs_blocked_retained_plan(...)`, if it becomes unused.
3. In `cartesian_pair_operator_plans_contract_runtests.jl`, collapse
   metadata-only materialization count clouds to one compact summary assertion,
   or delete synthetic `route_core_sidecars = false` branches if they only
   preserve stale blocker names.

Preserve:

- active unit-pair shape/order/family contract;
- pair-operator plan shape/order/transform matching contract;
- one compact blocked/readiness smoke if it protects a live module boundary;
- default nested runner coverage, if these files are still included there;
- fake-PQS and independent H2 PQS work.

Forbidden:

- do not edit source code;
- do not remove the core happy-path contract tests for unit pairs or pair
  operator plans;
- do not touch H2 independent source-plan/final-basis/H1 work;
- do not touch fake-PQS endpoint coverage;
- do not add new tests or new helper layers.

Line budget:

- This pass must be net-negative in `src + test + bin`.
- Target at least `60` net deleted lines.
- If you cannot safely delete at least `40` net lines, write `ATTENTION.md`
  explaining the blocker and exact remaining candidates.

Validation:

Default focused validation:

```text
julia --project=. test/nested/cartesian_unit_pairs_contract_runtests.jl
julia --project=. test/nested/cartesian_pair_operator_plans_contract_runtests.jl
git diff --check
```

If one of these tests takes more than 60 seconds or is mostly protecting the
stale branch being deleted, do not force obsolete shape. Run package load plus
the smallest surviving test or smoke, and explain which active contract remains
covered.

Do not run:

- `test/nested/integration_runtests.jl`;
- `test/nested/cartesian_pair_stage_low_order_policy_runtests.jl`;
- full Be2/CR2/Hamiltonian/HF/RHF routes.

Report:

- exact lines added/deleted;
- helper/test blocks deleted or collapsed;
- active contract preserved;
- validation and timings;
- deleted/simplified/quarantined/not-deleted accounting;
- exact remaining blocker/caller.

-- repo-manager@macmini
