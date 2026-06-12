Purpose:
  Feed the pass-091 PQS H1/J diagnostic driver slot from real route-owned
  complete core/shell inputs, or stop with the exact missing input boundary.
  Do not add more report fields.

Context:
  Pass 091 added a compact driver-facing H1/J diagnostic payload slot in
  `cartesian_assembly(...)` / `cartesian_report(...)`. The active driver state
  correctly remains blocked on missing complete core/shell inputs:
  region plan, source plan, final basis, H1 payload, axis weights, raw pair
  numerator terms, and Coulomb expansion.

Existing route-owned construction sequence to reuse, not duplicate:
  `test/nested/pqs_direct_retained_final_h1_runtests.jl` already demonstrates
  the complete one-center path:

  - build a reviewed parent axis bundle;
  - shellify with `CartesianShellification.shellify(...)`;
  - lower with `CartesianTerminalLowering.lower_terminal_regions(..., PQSLowering(...))`;
  - build `pqs_multilayer_shell_region_plan(...)`;
  - build `pqs_multilayer_shell_source_plan(bundles, region_plan; ...)`;
  - build `pqs_multilayer_complete_core_shell_final_basis(...)`;
  - build `pqs_multilayer_complete_core_shell_h1_payload(...)`.

Target:
  Add a small internal driver-stage helper, if feasible, that obtains those
  route-owned inputs from the existing driver objects (`parent`, `shells`,
  `units`, `transforms`, `pairs`, `recipe`) and passes them to the pass-091
  H1/J diagnostic payload helper.

Implementation guidance:
  1. First audit where the complete core/shell ingredients already live in the
     driver spine. In particular inspect:
       - `parent.parent_axis_bundle_object`;
       - terminal shellification/lowering data in `shells`, `units`, and
         `transforms`;
       - whether axis weights, raw pair-factor terms, and Coulomb expansion are
         already available from the same parent/axis bundle objects or need a
         named route-input seam.
  2. If all required ingredients are already present and the code is small,
     build a private helper that returns the complete core/shell H1/J diagnostic
     route payload and let `cartesian_assembly(...)` pass it into the pass-091
     payload slot.
  3. If any required ingredient is not already route-owned in the driver spine,
     do not synthesize it from a test fixture and do not hard-code side/core
     boxes. Instead stop with a precise blocker list: which driver stage should
     own each missing input and why.
  4. If implementing, keep the payload diagnostic/internal and blocked unless
     all route-owned inputs are present. Do not promote a physics fixture or
     acceptance gate.

Do not:
  - add RHF, SCF, Fock construction, density iteration, GTO, exports, artifacts,
    q/side ladder policy, fixture-rule policy, or acceptance status;
  - hard-code the `pqs_direct_retained_final_h1_runtests.jl` fixture into the
    driver;
  - use the explicit-box compatibility bridge as active driver authority;
  - use WL/fixed-block data as active authority;
  - add another report-only placeholder or more report fields;
  - expand broad route-driver/report/materialization tests;
  - run slow broad tests as routine validation.

Test policy:
  - Prefer a compact dry-run smoke or focused H1/J gate only.
  - If a new test is necessary, it must replace or shrink existing probe/test
    pressure and protect the live driver-owned input seam.
  - Do not add helper-vocabulary assertions.
  - Do not use `test/nested/pqs_source_box_route_driver_report_runtests.jl` or
    the assembly/report low-order policy files as per-pass validation unless
    you first justify why the touched contract requires those slow gates.

Validation:
  - Run the focused smoke/probe for this seam.
  - Run `julia --project=. test/nested/pqs_direct_retained_final_h1_runtests.jl`
    if H1/final-basis construction is touched.
  - Run `julia --project=. -e 'using GaussletBases; println("load ok")'`.
  - Run `git diff --check`.

Report back:
  - files changed;
  - whether the H1/J driver slot is now fed by real route-owned inputs;
  - if not, exact missing input boundary by driver stage;
  - whether any test/probe glue became less necessary;
  - validation run;
  - deletion/shrinkage report.

Deletion/shrinkage report required:
  - what old code, test, metadata, or compatibility path became unnecessary;
  - what was deleted or simplified;
  - if nothing was deleted, why no existing surface was made obsolete;
  - whether any new test replaces/shrinks older coverage or is genuinely new
    live-contract coverage;
  - any remaining stale or duplicate surfaces to retire next.

No escalation:
  Do not request UI escalation. If a command needs permission or an external
  condition blocks progress, write `.agent_handoffs/ATTENTION.md` with the
  blocker and stop.

Write the response to:
  - `.agent_handoffs/response.092.md`
  - `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.092.md`

-- repo-manager@macmini
