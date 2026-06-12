Purpose:
  Shrink the remaining explicit-box authority in the PQS H1 gate. Pass 076
  added a shellification/lowering-backed region plan and source-plan entry
  point, but the active H1 operator path in the test still uses the explicit
  `pqs_multilayer_shell_source_plan(bundles, core_box, outer_box; ...)` plan.

Task:
  Make the region-plan source plan the active fixture path in
  `test/nested/pqs_direct_retained_final_h1_runtests.jl`.

Implementation guidance:
  - Build shellification and terminal lowering as in pass 076.
  - Build `region_plan = pqs_multilayer_shell_region_plan(...)`.
  - Build the active `plan` via
    `pqs_multilayer_shell_source_plan(bundles, region_plan; ...)`.
  - Keep the explicit-box source plan only as a compact bridge/oracle
    comparison, not as the object used for final-basis construction,
    support operators, or H1.
  - Compare only compact facts between the active region-backed plan and the
    explicit-box bridge: status, support counts, layer count, shell retained
    count, and final retained count or final overlap identity if already cheap.

Support-space guardrail:
  Keep existing support-space dense overlap/kinetic/nuclear helpers limited to
  this H1 seam. Do not extend them toward density, RHF, driver route, or a
  scalable PQS operator algorithm in this pass.

Do not:
  - add new physics values, fixture-rule policy, IDA, density-density, RHF, or
    driver wiring;
  - add exports or artifacts;
  - add a broad new test file;
  - assert broad metadata vocabulary;
  - delete the explicit-box bridge entry point yet;
  - request UI escalation. In unattended baton mode, write
    `.agent_handoffs/ATTENTION.md` and stop if permission is genuinely needed.

Validation:
  - `julia --project=. -e 't = @elapsed include("test/nested/pqs_direct_retained_final_h1_runtests.jl"); println("elapsed_s=", t)'`;
  - `julia --project=. -e 'using GaussletBases; println("load ok")'`;
  - `git diff --check`.

Report:
  - confirm active H1 path now uses the region-plan source plan;
  - describe the remaining explicit-box bridge comparison;
  - validation run;
  - deletion/shrinkage report:
      - what explicit-box/test authority was removed or downgraded;
      - what explicit-box bridge remains and why;
      - whether the test got smaller or only changed authority.

-- repo-manager@macmini
