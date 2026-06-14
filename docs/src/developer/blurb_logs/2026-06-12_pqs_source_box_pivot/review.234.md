# Pass 234 manager review - accepted

Accepted under the manager-approved line-budget exception.

The pass generated the independent H2 PQS support-region fingerprint from route
geometry via `CartesianShellification.raw_terminal_geometry(...)`. It replaced
the pass-232 target-constant blocker with generated support-plan authority while
keeping all retained-transform and physics surfaces blocked.

Key checks:

- `target/support_plan_status =
  :available_independent_pqs_support_region_plan`.
- `target/support_plan_authority =
  :cartesian_shellification_route_geometry`.
- `target/support_counts_generated = true`.
- `target/support_counts_source =
  :cartesian_shellification_route_geometry`.
- Support counts are generated in target order: `(275, 578, 362)`.
- Coverage/disjointness counters are clean:
  `coverage_complete = true`, duplicate/missing/outside counts all zero.

Manager validation:

- Reviewed the source diff and response.
- Ran `git diff --check`.
- Ran package load:
  `julia --project=. -e 't = @elapsed begin using GaussletBases; println("load ok") end; println("elapsed_s=", t)'`
  with `elapsed_s=0.648348375`.
- Did not rerun the focused 66.6s driver artifact check; doer reported it
  passed.

Line-budget exception:

- Scoped diff was `118` added / `10` deleted, net `+108`.
- This is within the approved pass-234 exception after the earlier
  `ATTENTION.md`.
- Future passes should use the old-flat-path audit to pay down this line debt.

Guardrail:

- This is support-region authority only. Source plan, retained transforms,
  final basis, H1, H1-J, RHF, supplements, CR2, export, and public API remain
  blocked.

Next:

- Use the old-flat-path deletion audit for pass 235. The best target is the
  route-shadow PQS/PQS/product density-density fixture path and its slow
  integration-test assertions.

-- repo-manager@macmini
