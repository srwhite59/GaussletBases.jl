Purpose:
  Feed the next PQS H1/J driver-slot input seam: build the complete core/shell
  PQS region/source-plan inputs from the selected terminal shellification and
  lowering already carried by the driver.

Context:
  Pass 093 made one-center `:pqs_source_box` select terminal shellification.
  The dry-run now reports terminal shellification and terminal lowering
  available, while H1/J remains blocked on:

  - `pqs_multilayer_shell_region_plan`
  - `pqs_multilayer_shell_source_plan`
  - final basis
  - H1 payload
  - axis weights
  - raw pair terms
  - Coulomb expansion

Target:
  Add, if small and clean, a transform/unit-stage seam that builds:

  - `pqs_multilayer_shell_region_plan(shellification_plan, lowering_plan; ...)`
  - `pqs_multilayer_shell_source_plan(parent_axis_bundle_object, region_plan; ...)`

  from route-owned driver objects. The likely source is
  `terminal_route_state.shellification_plan`,
  `terminal_route_state.lowering_plan`, and `parent.parent_axis_bundle_object`
  (or the same bundle object propagated through the driver).

Implementation guidance:
  1. First audit whether `cartesian_transforms(...)` or another existing
     transform-stage helper already has access to both:
       - the selected terminal shellification/lowering plans; and
       - the parent axis bundle object needed by `pqs_multilayer_shell_source_plan`.
  2. If access is already present and the code is small, materialize a compact
     diagnostic/internal complete core/shell source-plan payload:
       - region plan;
       - source plan;
       - compact status/blocker/summary;
       - no final-basis/H1/J expansion yet unless the existing H1/J driver slot
         can consume the payload without adding fixture policy.
  3. If the parent axis bundle object is not available at the transform stage,
     do not rebuild it from scratch and do not hard-code the test fixture.
     Stop with the exact stage that should carry it.
  4. Keep the explicit-box bridge reference-only. The active source plan must
     use `pqs_multilayer_shell_source_plan(bundles, region_plan; ...)`.

Do not:
  - add RHF, SCF, Fock construction, density iteration, GTO, exports, artifacts,
    acceptance status, q/side ladder policy, or fixture tuning;
  - build final basis or H1 unless it is a tiny consequence of feeding an
    already available route-owned source plan into existing helpers;
  - use the explicit-box bridge as active authority;
  - use WL/fixed-block data as active authority;
  - add more report-only fields;
  - grow `test/nested/pqs_direct_retained_final_h1_runtests.jl`;
  - run broad CPBM/report/materialization tests as routine validation.

Test policy:
  - Prefer a compact dry-run smoke proving the region/source-plan status.
  - Do not add a permanent test unless it replaces/shrinks existing probe
    pressure and protects this live seam.
  - No helper-vocabulary assertions.

Validation:
  - Run the compact dry-run smoke/probe for the source-plan seam.
  - Run `julia --project=. -e 'using GaussletBases; println("load ok")'`.
  - Run `git diff --check`.
  - Avoid slow broad tests unless explicitly justified.

Report back:
  - files changed;
  - whether driver-owned region/source plans are now available;
  - exact source-plan status/blocker;
  - whether the H1/J diagnostic slot advanced or remains blocked;
  - what remains missing for final-basis/H1/J;
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
  - `.agent_handoffs/response.094.md`
  - `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.094.md`

-- repo-manager@macmini
