# Pass 244 blurb - shrink low-order report CRC alias assertions

Role: repo-doer.

Read before editing:

- `AGENTS.md`
- `BlurbStyle.md`
- `docs/src/developer/pqs_manager_running_log.md`
- `docs/src/developer/old_flat_cartesian_retirement_audit_2026-06-14.md`
- `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.243.md`
- `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/review.243.md`

Task type: cleanup/test shrink.

Purpose:

Shrink duplicate low-order report/CRC sidecar assertion pressure. This is a
test cleanup pass, not a source behavior pass. The current active authority is
the structured low-order route/RouteCore summaries; the flat
`low_order_*` report aliases and exact print strings are compatibility output,
not the contract to exhaustively inventory in a large test.

Deletion candidate:

```text
file:
  test/nested/cartesian_report_stage_low_order_policy_runtests.jl

nearby compact compatibility test:
  test/nested/pqs_source_box_route_driver_crc_print_line_runtests.jl
```

Subagent/read-only audit summary:

```text
risk class:
  green for shrinking duplicate alias/print assertions only.

expected line savings:
  about 100-180 test lines.

obsolete pressure:
  duplicated CRC sidecar readiness through structured summary fields,
  low_order_* flat report aliases, and exact summary stdout strings.

replacement/current authority:
  structured low_order_route_summary / RouteCore summaries;
  flat report aliases remain compatibility output only.
```

Prioritize shrinking these areas:

```text
atom-growth CRC pair-operator alias/count cloud:
  around the route_core_typed_pair_operator_* assertions near lines 499-568

summary stdout exact-string cloud:
  around the cartesian_print_summary assertions near lines 889-944
```

What to keep:

- one compact assertion that atom-growth report/summary has CRC final units `8`
  and pairs `36`;
- one compact assertion that typed pair-operator metadata is blocked on
  `:aggregate_subtree_operator_plan_required`;
- one print compatibility assertion if needed, preferably via the smaller
  `pqs_source_box_route_driver_crc_print_line_runtests.jl` coverage or one
  short `occursin("CRC", ...)` smoke;
- existing active route summary checks that protect material route selection;
- fake-PQS and independent H2 PQS coverage untouched.

What to delete/collapse:

- duplicate equality assertions proving every flat
  `low_order_route_core_typed_pair_operator_*` field mirrors the structured
  summary;
- detailed count-by-field assertions for
  `:metadata_only_not_materialized`,
  `:blocked_metadata_only_not_materialized`, and similar route-shadow
  vocabulary, unless needed for one compact blocker smoke;
- exact print-string assertions that repeat scalar fields already tested
  elsewhere;
- assertions whose only value is preserving private report-key spelling for
  old CRC sidecar internals.

Forbidden:

- do not edit source code;
- do not remove active route selection checks;
- do not touch H2 independent source-plan/final-basis/H1 work;
- do not touch fake-PQS endpoint coverage;
- do not broaden into the slow `cartesian_pair_stage_low_order_policy_runtests.jl`
  integration gate;
- do not add new tests.

Line budget:

- This pass must be net-negative in `src + test + bin`.
- Target at least `100` net deleted lines.
- If you cannot safely delete at least `60` net lines, write `ATTENTION.md`
  explaining the blocker and exact remaining candidates.

Validation:

Default focused validation:

```text
julia --project=. test/nested/cartesian_report_stage_low_order_policy_runtests.jl
julia --project=. test/nested/pqs_source_box_route_driver_crc_print_line_runtests.jl
git diff --check
```

If the report-stage test takes more than 60 seconds or is mostly protecting the
stale assertions being deleted, do not force the obsolete test shape. Instead
run package load plus the smaller CRC print-line test, and explain which live
route-summary checks remain.

Do not run:

- `test/nested/integration_runtests.jl`;
- `test/nested/cartesian_pair_stage_low_order_policy_runtests.jl`;
- full Be2/CR2/Hamiltonian/HF/RHF routes.

Report:

- exact lines added/deleted;
- assertion blocks deleted/collapsed;
- live summary/print contract preserved;
- validation and timings;
- deleted/simplified/quarantined/not-deleted accounting;
- exact remaining caller/blocker.

-- repo-manager@macmini
