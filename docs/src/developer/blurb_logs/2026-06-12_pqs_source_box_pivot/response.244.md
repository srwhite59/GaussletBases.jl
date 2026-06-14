Pass 244 complete.

Files changed:
- `test/nested/cartesian_report_stage_low_order_policy_runtests.jl`

Exact lines added/deleted:
- `git diff --numstat -- src test bin`
  - `4  121  test/nested/cartesian_report_stage_low_order_policy_runtests.jl`
- Net `src + test + bin`: `-117`.

Assertion blocks deleted/collapsed:
- Deleted `_cartesian_report_stage_count_by_field`, which only supported the stale count-by-field assertion cloud.
- Collapsed the atom-growth CRC typed pair-operator block by deleting:
  - duplicate per-field count assertions for `:metadata_only_not_materialized`,
    `:blocked_metadata_only_not_materialized`, blocker counts, and source-path counts;
  - duplicate family-count/materialization assertions;
  - flat `low_order_route_core_typed_pair_operator_*` report alias mirror assertions;
  - detailed materialization-readiness count mirrors.
- Preserved compact typed pair-operator blocker coverage:
  - `route_core_typed_pair_operator_plan_inventory_status == :blocked_route_core_pair_operator_plan_inventory`
  - blocker `:aggregate_subtree_operator_plan_required`
  - materialization readiness blocker `:aggregate_subtree_operator_plan_required`
- Collapsed exact `cartesian_print_summary` stdout assertions by deleting scalar report-field string checks and the full CPB/PQS prototype print-string assertion.

Live summary/print contract preserved:
- The test still checks atom-growth RouteCore summary/report final-unit and pair counts:
  - final units `8`
  - pairs `36`
- The print smoke now checks only compact compatibility substrings:
  - `CRC sidecars: final units 8, pairs 36`
  - `CRC pair families: `
  - `CRC pair-operator plan: ready metadata plan`
  - `CRC typed pair-operator inventory:`
  - `:aggregate_subtree_operator_plan_required`

Validation and timings:
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - passed, printed `load ok`
- `julia --project=. test/nested/pqs_source_box_route_driver_crc_print_line_runtests.jl`
  - passed: `1` test, `0.0s` reported by Julia
- `git diff --check`
  - passed
- Attempted default report-stage gate:
  - `julia --project=. test/nested/cartesian_report_stage_low_order_policy_runtests.jl`
  - exceeded the 60s threshold and exposed stale assertion shape rather than a live compact contract.
  - In the edited target block, the live object retains the aggregate-subtree blocker but no longer supports the old available/count cloud; I removed those stale target assertions.
  - It also failed on unrelated pre-existing terminal/report exact-field assertions outside this pass target, including missing `terminal_shellification_scaffold` and terminal lowering-contract fields on the current assembly object.
  - Per the blurb fallback, I did not force that obsolete test shape. I terminated the long focused process with `pkill -f cartesian_report_stage_low_order_policy_runtests.jl` after it continued streaming stale failures.

Git status:
```text
## main...origin/main
 M test/nested/cartesian_report_stage_low_order_policy_runtests.jl
```

Deletion/shrinkage accounting:
- deleted: count-by-field helper; duplicate CRC typed pair-operator count/materialization/source-path assertion cloud; exact summary stdout scalar/string cloud.
- simplified: retained one compact blocked typed-pair metadata smoke and compact CRC print substrings.
- quarantined: none.
- not deleted because: active route selection, RouteCore final-unit/pair-count checks, terminal route checks outside the requested target area, fake-PQS coverage, and independent H2 PQS coverage were left untouched.
- exact remaining caller/blocker: `test/nested/cartesian_report_stage_low_order_policy_runtests.jl` itself remains a stale broad report-stage gate; unrelated terminal-shellification exact-field assertions still block using it as focused validation.

-- repo-doer@macmini
