Purpose:
  Use the pass-094 complete core/shell source-plan payload to build the
  driver-owned final basis and H1 payload, without entering J/RHF.

Context:
  Pass 094 made `cartesian_assembly(...)` build a local complete core/shell
  source-plan payload from terminal shellification/lowering and the parent axis
  bundle. The H1/J missing-input list now excludes region plan, source plan, and
  Coulomb expansion. It remains blocked on final basis, H1 payload, axis
  weights, and raw pair terms.

Target:
  Add, if small and clean, assembly-stage construction of:

  - `pqs_multilayer_complete_core_shell_final_basis(source_plan; ...)`
  - `pqs_multilayer_complete_core_shell_h1_payload(source_plan; final_basis, ...)`

  using the route-owned source plan and parent/route-owned center/Coulomb data.
  Feed those into the existing H1/J diagnostic payload helper so its missing
  list advances to the density/J inputs only.

Implementation guidance:
  1. Reuse the local source-plan payload from pass 094. Do not recompute the
     region/source plan through a second path.
  2. Audit how the driver already represents nuclear center records. Use those
     if present. If no route-owned center records are available, stop with a
     precise blocker instead of hard-coding the H1 test fixture.
  3. For electron-nuclear factors, prefer existing parent axis bundle data or
     the already-supported centered axis-layer route. Do not invent a new GTO
     or quadrature path.
  4. Keep the dense final-basis/H1 payload internal to assembly or summarized
     compactly. Do not add broad report fields or retain dense payloads in the
     public report unless an existing local pattern requires it.
  5. Do not build axis weights/raw pair numerator terms or self-Coulomb J in
     this pass unless those inputs are already route-owned and connecting them
     is a trivial local extension. The intended primary result is H1.

Do not:
  - add RHF, SCF, Fock construction, density iteration, GTO, exports, artifacts,
    acceptance status, q/side ladder policy, or fixture tuning;
  - hard-code `pqs_direct_retained_final_h1_runtests.jl` fixture centers,
    boxes, or bundle construction into the driver;
  - use the explicit-box bridge as active authority;
  - use WL/fixed-block data as active authority;
  - add more report-only fields;
  - grow `test/nested/pqs_direct_retained_final_h1_runtests.jl`;
  - run broad CPBM/report/materialization tests as routine validation.

Test policy:
  - Prefer a compact dry-run smoke proving final basis and H1 payload are
    available, with H1/J still blocked only on density/J inputs.
  - Do not add a permanent test unless it replaces/shrinks existing probe
    pressure and protects this live driver-owned seam.
  - No helper-vocabulary assertions.

Validation:
  - Run the compact dry-run smoke/probe for the final-basis/H1 seam.
  - Run `julia --project=. test/nested/pqs_direct_retained_final_h1_runtests.jl`
    if final-basis/H1 construction helpers are touched.
  - Run `julia --project=. -e 'using GaussletBases; println("load ok")'`.
  - Run `git diff --check`.
  - Avoid slow broad tests unless explicitly justified.

Report back:
  - files changed;
  - whether driver-owned final basis and H1 payload are now available;
  - H1/J diagnostic status and remaining missing inputs;
  - exact blocker if final basis/H1 could not be built;
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
  - `.agent_handoffs/response.095.md`
  - `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.095.md`

-- repo-manager@macmini
