# Pass 237 blurb - independent H2 PQS retained-rule readiness plan

Role: repo-doer.

Read before editing:

- `AGENTS.md`
- `BlurbStyle.md`
- `docs/src/developer/pqs_manager_running_log.md`
- `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.236.md`
- `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/review.236.md`

Task:

Add a compact retained-rule/readiness plan for the independent H2 PQS route.
This is metadata/readiness only. Do not materialize source coefficients, final
basis, H1, H1-J, RHF, supplements, CR2, export, or public API.

Manager decision from pass 236:

Independent PQS should not force the fake/WL retained counts `(251, 98, 114)`.
Use the route-owned authorities that currently exist:

```text
support counts:  (275, 578, 362)
retained counts: (275, 98, 98)
final target dimension/readiness: 471

:atom_contact_core => :direct_source_modes
:shared_shell_1    => :pqs_boundary_comx_product_modes
:shared_shell_2    => :pqs_boundary_comx_product_modes
```

Allowed implementation:

- Add a compact private helper/object near the independent H2 target payload,
  for example:

```julia
_PQSIndependentH2CoreShellRetainedRulePlan
_pqs_source_box_route_driver_independent_h2_retained_rule_plan(...)
```

- Compose existing authority:
  - support plan from pass 234;
  - direct source-mode rule for `:atom_contact_core`;
  - `CartesianRawProductSources.pqs_boundary_product_mode_retained_rule(...)`
    or the existing terminal-lowering retained-rule summary for each q=5 shared
    shell.
- Report status, blocker, support order/counts, retained order/counts,
  expected final dimension, per-unit retained-rule authority, and missing
  objects.

Required semantics:

- retained-rule/readiness plan status is available only if support plan is
  generated and counts/rules match `(275, 98, 98)`.
- retained counts must be `(275, 98, 98)`.
- expected final dimension/readiness target must be `471`.
- artifact/reporting must make clear this is independent PQS retained-rule
  readiness, not fake-PQS/WL reproduction.
- `fake_pqs/enabled = false` and
  `source_backed_fixed_source_oracle_used = false` remain true/false as before.
- physics endpoint remains false.

Forbidden:

- no coefficient matrices;
- no final basis;
- no H1, H1-J, RHF, supplements, CR2, export, public API;
- no fake-PQS/WL coefficient matrices or fixed-source retained transforms;
- no claim that old counts `(251, 98, 114)` are independent PQS.

Line-budget:

- Keep `src + test + bin` net-negative if possible.
- Use same-surface cleanup first: stale fake/WL retained-count aliases,
  duplicate blocker fields, or route-skeleton helper-name/report pressure that
  becomes obsolete after the retained-rule plan exists.
- If honest net-negative implementation is blocked, write `ATTENTION.md` with
  exact diff totals and deletion candidates. Do not delete scientific endpoint
  tests or fake-PQS guard fields.

Validation:

- Run package load.
- Run the smallest focused independent-input artifact/readiness check.
- Run `git diff --check`.
- Report timings for any Julia command over 60 seconds.

Report:

- retained-rule plan status/counts/authority;
- expected final dimension/readiness target;
- forbidden surfaces avoided;
- validation and timings;
- scoped line-budget totals;
- deleted/simplified/quarantined/not-deleted accounting.

-- repo-manager@macmini
