Pass 077 complete.

Active H1 path:

- `test/nested/pqs_direct_retained_final_h1_runtests.jl` now builds the active `plan` through:
  - `CartesianShellification.shellify(...)`;
  - `CartesianTerminalLowering.lower_terminal_regions(..., PQSLowering(q = 5))`;
  - `pqs_multilayer_shell_region_plan(...)`;
  - `pqs_multilayer_shell_source_plan(bundles, region_plan; ...)`.
- The active final-basis construction, support states, support overlap/kinetic/nuclear matrices, and H1 solve now consume that region-plan-backed `plan`.

Remaining explicit-box bridge comparison:

- The old explicit-box source plan is still built as `explicit_box_plan`.
- It is used only for compact bridge/oracle checks:
  - status;
  - layer count;
  - core support count;
  - shell support count;
  - shell final retained count;
  - final retained count through `explicit_box_final_basis`.
- It is no longer the active object for final-basis construction, support operators, or H1.

Validation run:

- `julia --project=. -e 't = @elapsed include("test/nested/pqs_direct_retained_final_h1_runtests.jl"); println("elapsed_s=", t)'`
  - passed: 44 tests;
  - reported `elapsed_s=5.942399084`.
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - passed.
- `git diff --check`
  - passed.

Deletion/shrinkage report:

- Explicit-box/test authority removed or downgraded:
  - `fixture.plan` is no longer the explicit-box source plan;
  - active final-basis construction and all active H1 support-operator paths now use the shellification/lowering-backed source plan.
- Explicit-box bridge remaining:
  - `explicit_box_plan` remains because pass 076 intentionally kept the explicit-box entry point working as a bridge;
  - it is only used to verify compact equivalence against the active region-backed path.
- Tests:
  - no new test file was added;
  - the existing H1 gate was minimally changed, not broadened;
  - the test did not get smaller in line count, but its authority moved away from explicit PQS box arithmetic.
- Remaining stale/duplicate surface to retire next:
  - the explicit-box source-plan entry point remains a bridge and still owns duplicate box-depth/layer-box arithmetic internally;
  - after another source-plan implementation pass, the explicit-box bridge can be quarantined more clearly or reduced to a compatibility wrapper.

-- repo-doer@macmini
