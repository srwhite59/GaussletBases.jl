Review 029 - manager recovery

No `.agent_handoffs/response.029.md` appeared after an extended wait. The
visible worktree change was narrow and matched the blurb's cleanup target, so
this review treats pass 029 as a manager recovery from the live diff rather
than a normal doer handoff.

Accepted change:

- `test/nested/bond_aligned_diatomic_high_order_recipe_opt_in_source_construction_integration_runtests.jl`
  deletes 37 helper-vocabulary assertions from the current-route safe-term
  oracle section.

Why this is acceptable:

- The test still exercises `_pqs_current_route_safe_term_matrices(...)` as a
  private diagnostic/oracle path.
- It still checks matrix dimensions, finite matrices, matrix agreement with the
  oracle, debug-only status, fixed-block authority comparison, term count, and
  max authority error.
- The removed assertions were mostly detailed no-claim flags, timing fields,
  pair-count subfields, and a negative unsupported-term throw check. Those
  preserved old helper vocabulary more than live behavior.

Validation:

- `julia --project=. -e 'using GaussletBases; println("load ok")'` passed.
- `git diff --check` passed.
- Direct execution of the edited integration file is not valid because it
  depends on top-level harness fixtures.
- The proper slow nested harness was attempted with:
  `TMPDIR=/private/tmp GAUSSLETBASES_TEST_GROUPS=nested GAUSSLETBASES_SLOW_TESTS=1 ... test/runtests.jl`.
  It failed before reaching the edited file with a stack overflow in
  `test/nested/pqs_source_box_route_driver_report_runtests.jl`, at
  `_pqs_route_driver_check_materialization_status(...)`. The stack trace shows
  specialization on a giant route-report `NamedTuple`, which is unrelated to
  the assertion-deletion diff and is now the next live validation blocker.

Deletion/shrinkage:

- Deleted 37 assertions.
- Added no new tests.
- Added no new source, docs, metadata, adapters, or artifacts.
- Remaining stale surface: the current-route safe-term helper remains as an
  oracle/debug path. It is not yet deletable because the integration test still
  uses it for old-route authority comparison.

Next pass:

- Fix the route-driver test helper stack overflow without adding new tests.
  The likely target is preventing specialization on giant route-report
  `NamedTuple` arguments in `_pqs_route_driver_check_materialization_status(...)`
  or an adjacent helper.

-- repo-manager@macmini
