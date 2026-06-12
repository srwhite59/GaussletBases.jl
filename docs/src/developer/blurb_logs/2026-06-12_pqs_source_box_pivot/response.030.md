Response 030 - manager recovery implementation

No doer response appeared for pass 030, and the worktree stayed clean. The
manager recovered the pass directly because the validation blocker was precise
and did not require a design decision.

Changed:

- `test/nested/pqs_source_box_route_driver_report_runtests.jl`
  - Added `@nospecialize pqs_report white_lindsey_report` at the start of
    `_pqs_route_driver_check_materialization_status(...)`.

Why:

- The slow nested harness failed before reaching the pass 029 edited file with
  a `StackOverflowError` while compiling/specializing
  `_pqs_route_driver_check_materialization_status(...)`.
- The stack trace showed specialization on a giant route-driver report
  `NamedTuple`.
- The helper only needs to inspect fields at runtime; it does not benefit from
  full specialization on the complete report type.

Validation:

- `julia --project=. test/nested/pqs_source_box_route_driver_report_runtests.jl`
  passed: 1189 checks, about 13m43s.
- `julia --project=. -e 'using GaussletBases; println("load ok")'` passed.
- `git diff --check` passed.

Deletion/shrinkage:

- No code/test surface became obsolete from this one-line fix.
- No new tests were added; the existing route-driver test was the right
  validation surface.
- The fix reduces compile/specialization pressure without creating new
  metadata, adapters, structs, or report fields.

Remaining cleanup:

- Giant route-report `NamedTuple` payloads remain a known pressure source.
  This pass only stopped one validation blocker; it did not redesign the
  route-report shape.
- The pass 029 edited integration file still needs validation through the slow
  nested harness now that the route-driver stack overflow is fixed.

-- repo-manager@macmini
