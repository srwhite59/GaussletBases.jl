Purpose:
  Implement the small shellification/lowering-backed boundary identified in
  pass 075, so multi-layer PQS source planning no longer has to be driven only
  by explicit `core_box` / `outer_box` arithmetic.

Target:
  Add a compact PQS region-plan seam that consumes:

  ```text
  CartesianShellification.ShellificationPlan
  CartesianTerminalLowering.lower_terminal_regions(..., PQSLowering(q = ...))
  ```

  and then add a new `pqs_multilayer_shell_source_plan(bundles, region_plan; ...)`
  entry point that reproduces the current explicit-box result for the existing
  complete core/shell H1 fixture.

Implementation guidance:
  1. Add the smallest module-owned object/helper needed. A name like
     `pqs_multilayer_shell_region_plan(...)` is fine.
  2. The region plan should carry only geometry/lowering facts:
       - direct core terminal region and core box;
       - ordered complete-shell layer records and their outer/inner boxes;
       - selected PQS lowering contracts/source CPBs for shell layers;
       - support coverage and duplicate/disjointness fingerprints;
       - provenance from `CartesianShellification` and
         `CartesianTerminalLowering`.
  3. It must not carry PQS descriptors, shell projection/Lowdin matrices,
     support operator blocks, H1, IDA, density-density, RHF, final-basis
     transfer, driver data, exports, or artifacts.
  4. Add `pqs_multilayer_shell_source_plan(bundles, region_plan; ...)` as a
     consumer. It may delegate to the existing explicit-box implementation for
     now, but the public authority for the boxes should be the region plan, not
     ad hoc test/probe arithmetic.
  5. Keep the existing explicit-box entry point working as a bridge; do not
     delete it yet.

Fixture to compare:
  Use the same geometry as `test/nested/pqs_direct_retained_final_h1_runtests.jl`:

  ```text
  parent count 7
  current_box = (1:7, 1:7, 1:7)
  inner_box = (2:6, 2:6, 2:6)
  raw source dims = 5 x 5 x 5
  complete final retained count = 223
  ```

  Build a shellification/lowering-backed region plan for the one-center case
  and verify that the new entry point matches the explicit-box plan on compact
  facts such as support counts, layer count, shell retained count, final
  retained count after `pqs_multilayer_complete_core_shell_final_basis`, and
  final overlap identity. Do not assert broad metadata vocabulary.

Support-space guardrail:
  The existing support-space dense kinetic/nuclear helpers are accepted as
  H1 seam/oracle machinery. Do not generalize them into the scalable PQS
  operator algorithm in this pass. This pass is about geometry/lowering
  authority, not operator scaling.

Do not:
  - add H1/RHF/IDA/density-density behavior;
  - add fixture-rule policy for `Z`, `d`, `s`, radius, q, or shell depth;
  - add driver wiring, exports, or acceptance gates;
  - add broad tests or helper-vocabulary tests;
  - remove the explicit-box bridge entry point yet;
  - request UI escalation. In unattended baton mode, write
    `.agent_handoffs/ATTENTION.md` and stop if permission is genuinely needed.

Test policy:
  Add or update one compact module/fixture check only. Prefer extending the
  existing H1 gate with a small comparison that the region-plan entry point
  reproduces the explicit-box plan, rather than creating a new broad test file.

Validation:
  - focused H1 gate or smaller focused test if you add one;
  - `julia --project=. -e 'using GaussletBases; println("load ok")'`;
  - `git diff --check`.

Report:
  - object/helper name and exact ownership boundary;
  - comparison result versus the explicit-box bridge;
  - validation run;
  - deletion/shrinkage report:
      - what responsibility moved out of ad hoc PQS box arithmetic;
      - what explicit-box/probe code remains as bridge only;
      - whether any tests were added, shrunk, or only minimally extended.

-- repo-manager@macmini
