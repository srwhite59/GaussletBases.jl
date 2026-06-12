Purpose:

Fix the slow-harness stack overflow that blocked validation of pass 029.

Why now:

Pass 029 only removed helper-vocabulary assertions, but the proper nested slow
harness could not reach the edited file. It failed earlier in
`test/nested/pqs_source_box_route_driver_report_runtests.jl` with a
`StackOverflowError` at `_pqs_route_driver_check_materialization_status(...)`.
The stack trace shows Julia specializing on an enormous route-report
`NamedTuple`. This is exactly the report-shape/compile-pressure issue we have
been trying to reduce.

Exact task:

Make the route-driver status test/helper avoid stack overflow from giant
route-report specialization.

Start with these surfaces:

- `test/nested/pqs_source_box_route_driver_report_runtests.jl`
- `_pqs_route_driver_check_materialization_status(pqs_report, white_lindsey_report)`
- the call at the `"Route-driver standard unit inventory report"` testset

Likely fix:

- Add a local no-specialization boundary such as `@nospecialize pqs_report
  white_lindsey_report` in the helper, or split the helper so the giant route
  report is passed through compact summary/status values instead of being a
  fully specialized argument.
- Prefer the smallest test-side fix that makes the existing validation pass.
- Do not move this into production source unless inspection proves the same
  stack overflow is in a live source path.

Do not:

- add new tests;
- broaden the integration runner;
- add metadata fields;
- change source-box/PQS route behavior;
- touch final-basis H1 implementation;
- implement IDA, density-density, RHF, drivers, exports, or artifacts;
- replace the giant report with a new giant struct.

Validation:

- Run the focused failing route-driver test file through a valid harness if it
  can be run directly; otherwise run the narrowest top-level command that
  includes it.
- If you must run the slow nested harness, use `TMPDIR=/private/tmp` to avoid
  the known artifact-path assertion failure:

```text
TMPDIR=/private/tmp GAUSSLETBASES_TEST_GROUPS=nested GAUSSLETBASES_SLOW_TESTS=1 julia --project=. test/runtests.jl
```

- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git diff --check`

Deletion/shrinkage report required:

- whether the fix deletes/simplifies helper assertions or just adds a
  no-specialization boundary;
- whether any old test/report surface became unnecessary;
- if nothing was deleted, why no existing surface became obsolete;
- whether any new test was avoided and why existing coverage is enough;
- remaining report-shape or giant-`NamedTuple` surfaces to retire next.

Report back:

- write `.agent_handoffs/response.030.md.tmp`, then atomically rename to
  `.agent_handoffs/response.030.md`;
- also write the curated copy to
  `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.030.md`;
- include files changed;
- include validation run;
- include deletion/shrinkage report;
- sign `-- repo-doer@macmini`.

-- repo-manager@macmini
