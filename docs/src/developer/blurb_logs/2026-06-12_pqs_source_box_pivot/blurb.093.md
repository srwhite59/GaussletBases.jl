Purpose:
  Start feeding the PQS H1/J driver slot at the first missing owner: shell-stage
  selection of a complete core/shell PQS shellification plan.

Context:
  Pass 092 found that the parent stage can provide a one-center parent axis
  bundle, but `:pqs_source_box` dry-run currently reports:

  - `low_order_shellization_source = :not_applicable`
  - `low_order_terminal_shellification_selected = false`

  Therefore the driver cannot yet produce `PQSLowering`,
  `pqs_multilayer_shell_region_plan(...)`, or
  `pqs_multilayer_shell_source_plan(...)`.

Target:
  Add, if small and clean, a shell-stage complete core/shell PQS selection seam
  for one-center `:pqs_source_box` diagnostics. It should make
  `cartesian_shells(...)` own or precisely block the
  `CartesianShellification.shellify(...)` output needed by the complete
  core/shell PQS path.

Implementation guidance:
  1. First inspect the existing shell-stage selection policy in
     `src/pqs_source_box_route_driver_helpers.jl`, especially
     `cartesian_shells(...)` and the low-order shellization helpers.
  2. Reuse existing `CartesianShellification` policy objects where possible.
     For the current one-center diagnostic, the reference fixture uses
     `OneCenterShellification(core_side = 5, q = 5)` over parent index axes.
  3. If the driver already has a route-input policy knob for terminal
     shellification, use that. If it does not, add only a narrow internal
     diagnostic policy/status for complete core/shell PQS; do not invent a broad
     fixture-rule framework.
  4. The shell stage should output either:
       - an available shellification plan with compact provenance/status, or
       - a precise blocker explaining the missing parent/index/policy input.
  5. Do not proceed into lowering, source-plan realization, final basis, H1/J,
     or RHF in this pass unless the shell-stage change is already trivial and
     explicitly local. The primary target is shell-stage ownership.

Do not:
  - add RHF, SCF, density iteration, Fock construction, GTO, exports, artifacts,
    acceptance status, q/side ladder policy, or fixture tuning;
  - hard-code the entire `pqs_direct_retained_final_h1_runtests.jl` fixture into
    the driver;
  - use the explicit-box bridge as active shellification authority;
  - use WL/fixed-block data as active authority;
  - add more H1/J report fields;
  - grow `test/nested/pqs_direct_retained_final_h1_runtests.jl` unless it
    replaces or shrinks older pressure;
  - run the broad CPBM contract, broad route-driver report file, or
    assembly/report low-order policy files as routine validation.

Test policy:
  - Prefer a compact dry-run smoke proving the shell-stage status changed from
    `:not_applicable` to either available or a more precise blocker.
  - Do not add a permanent test unless it protects this live shell-stage seam
    and replaces/shrinks existing probe pressure.
  - No helper-vocabulary assertions.

Validation:
  - Run the compact dry-run smoke/probe for the shell-stage seam.
  - Run `julia --project=. -e 'using GaussletBases; println("load ok")'`.
  - Run `git diff --check`.
  - Do not run slow broad tests unless the edit genuinely touches their
    contract and you explain why.

Report back:
  - files changed;
  - whether `:pqs_source_box` shell stage now selects or blocks a complete
    core/shell PQS shellification plan;
  - exact shell-stage status and blocker;
  - what remains missing for lowering/source/final-basis/H1/J;
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
  - `.agent_handoffs/response.093.md`
  - `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.093.md`

-- repo-manager@macmini
