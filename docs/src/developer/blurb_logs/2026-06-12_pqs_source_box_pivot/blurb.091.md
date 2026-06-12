Purpose:
  Make the driver consume the compact PQS H1/J seam before the seam grows.
  Pass 090 found the correct insertion point but did not add a placeholder
  report hook because the driver spine does not yet own a complete core/shell
  diagnostic route payload.

Target:
  Add, if it is small and mechanical, a driver-facing complete core/shell PQS
  diagnostic route payload that can be carried by `cartesian_assembly(...)` and
  summarized by `cartesian_report(...)`.

Existing route-owned seams to consume:
  - `pqs_multilayer_shell_region_plan(...)`
  - `pqs_multilayer_shell_source_plan(bundles, region_plan; ...)`
  - `pqs_multilayer_complete_core_shell_final_basis(...)`
  - `pqs_multilayer_complete_core_shell_h1_payload(...)`
  - `pqs_multilayer_support_weights(...)`
  - `pqs_multilayer_support_pair_raw_numerator_matrix(...)`
  - `pqs_multilayer_complete_core_shell_h1_j_payload(...)`

Implementation guidance:
  1. First audit the existing driver helper surfaces in
     `src/pqs_source_box_route_driver_helpers.jl`, especially
     `cartesian_assembly`, `cartesian_report`, and `cartesian_materialization`.
  2. If the seam is small, add a compact internal diagnostic route payload
     object/helper. It should represent the complete core/shell H1/J diagnostic
     route flow, not merely a report metadata placeholder.
  3. Let `cartesian_assembly(...)` carry the payload or an explicit blocked
     status. Let `cartesian_report(...)` expose only a compact summary:
     status, blocker, final dimension, H1 energy, self-Coulomb value, density
     gauge, and the key nonclaim flags.
  4. Keep `cartesian_materialization(...)` unchanged unless a very small
     summary-only handoff is already the established local pattern.
  5. If a real driver-facing route object is not small or the insertion point
     is ambiguous, stop and write the exact missing object fields and driver
     stage where it should enter. Do not add a placeholder summary.

Do not:
  - add RHF, SCF, Fock construction, density iteration, or an RHF acceptance
    route;
  - add fixture-rule policy, q/side ladder policy, side13 promotion, GTO, PQS
    density-density beyond the existing H1/J diagnostic seam, exports, or
    artifacts;
  - make WL/fixed-block data active authority;
  - set `driver_route_materialized = true` unless the driver actually owns and
    carries the complete core/shell diagnostic route payload;
  - add broad tests or helper-vocabulary assertions;
  - expand the large CPBM contract test as a development notebook;
  - turn `pqs_multilayer_complete_core_shell_h1_j_payload(...)` into a private
    SCF route.

Test policy:
  - Prefer updating one compact driver/route-stage test only if it protects the
    new driver-owned payload.
  - Do not add a test for a report-only placeholder.
  - A focused H1/J probe or existing compact PQS final-H1 gate is acceptable if
    needed to verify the payload, but avoid slow nested harnesses unless a
    touched contract genuinely requires them.

Validation:
  - Run the focused test/probe touched by the change, if any.
  - Run `julia --project=. -e 'using GaussletBases; println("load ok")'`.
  - Run `git diff --check`.

Report back:
  - files changed;
  - whether a real driver-facing diagnostic route payload was added, or why it
    was not;
  - where the payload enters the driver spine;
  - compact report fields exposed, if any;
  - what remains private/diagnostic;
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
  - `.agent_handoffs/response.091.md`
  - `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.091.md`

-- repo-manager@macmini
