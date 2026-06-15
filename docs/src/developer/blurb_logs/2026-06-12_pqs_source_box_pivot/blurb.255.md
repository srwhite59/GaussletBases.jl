Pass 255 - shrink atom-growth report-stage RouteCore mirrors

Context:
- Current HEAD should include `a753bada Shrink terminal report-stage low-order aliases`.
- Pass 254 cleaned terminal-shellification report-stage alias mirrors and added
  the required passes 250-254 medium-term checkpoint.
- This pass should continue cleanup pressure only. No H2/PQS physics,
  supplement, CR2/export, public API, or driver input work.

Target:
- `test/nested/cartesian_report_stage_low_order_policy_runtests.jl`
- Focus on the atom-growth report-stage block around the
  `low_order_route_core_pair_operator_preflight` and
  `low_order_route_core_pair_operator_plan` assertions.

Task:
Delete stale atom-growth report-stage RouteCore preflight/plan mirror assertions
that preserve private nested object shape instead of active behavior.

Good deletion candidates:
- Nested preflight object field mirrors:
  - `low_order_route_core_pair_operator_preflight.route_core_final_unit_count`
  - `low_order_route_core_pair_operator_preflight.route_core_pair_count`
  - `!low_order_route_core_pair_operator_preflight.operator_blocks_materialized`
- Nested plan object field mirrors:
  - `low_order_route_core_pair_operator_plan.route_core_final_unit_count`
  - `low_order_route_core_pair_operator_plan.route_core_pair_count`
  - `!low_order_route_core_pair_operator_plan.operator_blocks_materialized`
  - `!low_order_route_core_pair_operator_plan.hamiltonian_matrices_materialized`
- Exact flat mirror comparisons such as
  `low_order_route_core_pair_family_counts == summary.route_core_pair_family_counts`
  if compact pair count/family-count checks remain nearby.

Keep:
- Atom-growth route selected/not terminal/not legacy.
- Compact atom-growth materialization deferral checks.
- Pair inventory known/available, pair count `36`, and at least one compact
  family-count assertion.
- RouteCore summary/inventory availability/status and pair count.
- Pair-operator ready/status/blocker at summary/report level.
- Preflight/plan availability/status/blocker at the report level if those are
  still useful compact readiness facts.
- LW complete-shell enumeration checks.

Strict exclusions:
- Do not edit source.
- Do not change H2/PQS route behavior, RHF, H1-J, supplements, CR2/export,
  public API, driver inputs, or manifest files.
- Do not add replacement assertions for deleted nested-object mirrors.
- Do not run broad stale integration gates.

Validation:
- Run `git diff --check`.
- Run Julia parse smoke for the touched test file.
- If the focused report-stage test is known to complete quickly, it may be run,
  but it is not required for this stale-mirror cleanup.
- No package load is required unless source is touched.

Line budget:
- Strictly net-negative in `src + test + bin`.

Report:
- Exact assertion families deleted.
- What compact atom-growth RouteCore readiness coverage remains.
- Validation commands.
- Scoped line count for `src + test + bin`.
- Deletion/shrinkage result:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

-- repo-manager@macmini
