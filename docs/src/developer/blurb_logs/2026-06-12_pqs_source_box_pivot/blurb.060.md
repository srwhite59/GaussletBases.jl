Purpose:
  Implement the narrow route-owned multi-layer PQS shell source plan discovered
  in pass 059, then run only a final-basis/H1 smoke for the WL-aligned side13
  fixture.

Context:
  Pass 059 showed repeated legal one-cell shell descriptors can represent the
  WL-aligned side13 three-shell PQS geometry:

  ```text
  parent count: 13
  mapping:      AsinhMapping(c = 0.1, s = 1.0, tail = 10)
  core:         (4:10)^3
  shell 1:      (3:11)^3 / (4:10)^3
  shell 2:      (2:12)^3 / (3:11)^3
  shell 3:      (1:13)^3 / (2:12)^3
  ```

  The collapsed shell sector had support count `1854`, retained count `1206`,
  final dimension `1549`, and final overlap identity error about `5.51e-13`.

Architecture guardrail:
  The public route driver spine is the executable sketch of the interface:

  ```text
  system -> recipe -> parent -> shells -> units -> transforms -> pairs
  -> assembly -> report -> materialization
  ```

  The new multi-layer PQS seam should be shaped so it can later enter the
  shells/transforms/final-basis stages. Do not build another private physics
  lane or special driver. Lower modules should not call upward into driver
  helpers.

Exact task:
  Add one narrow production helper/object for a multi-layer PQS shell source
  plan. Naming can follow the codebase, but the concept should be clear, e.g.

  ```text
  pqs_multilayer_shell_source_plan(...)
  ```

  It should:
  - accept an ordered core box and outer/current box or explicit shell layer
    boxes;
  - build repeated legal one-cell projected q-shell descriptors from the core
    outward;
  - call the existing shell-realization plan for each layer;
  - validate shell support disjointness and core/shell disjointness;
  - validate combined support covers the intended parent/current box;
  - concatenate shell support rows and block-diagonal shell final coefficients
    into one collapsed shell sector;
  - return a compact route object plus summary fields needed by the complete
    core/shell final-basis helper.

  Use the existing complete core/shell final-basis helper after the plan is
  materialized. Do not generalize that helper unless the implementation proves
  it is necessary.

Smoke probe:
  Run a side13 final-basis/H1 smoke using the new route-owned plan:

  ```text
  parent count: 13
  AsinhMapping(c = 0.1, s = 1.0, tail = 10)
  core:         (4:10)^3
  outer box:    (1:13)^3
  shell layers: 3
  ```

  Report final dimension, overlap identity error, Z=1 H1, Z=2 H1, and timing.
  Do not run density-density or RHF in this pass.

Spacing caveat:
  Record only a caveat, not a new research task: PQS may not follow the same
  q/ns accuracy intuition as WL, and the relation among `Z`, core spacing `d`,
  mapping parameter `s`, radius, and shell depth remains provisional.

Do not:
  - run RHF;
  - add GTO, driver wiring, exports, or artifacts;
  - promote q=9/q=11 or the side13 H1 smoke as acceptance;
  - make fixed-block matrices route authority;
  - add broad tests or metadata-vocabulary assertions.

Tests:
  Prefer one compact module-contract test only if needed to protect the new
  plan object. It should check counts, disjoint support, parent coverage, and
  final-basis availability for a small fixture. Avoid asserting every report
  field. If a tmp/work smoke is enough for this pass, do not add a permanent
  test yet.

Docs:
  Update the PQS near-term plan with the implemented seam and side13 H1 smoke
  result if successful. Keep the note concise and non-acceptance.

Validation:
  - run the new focused test if added, otherwise the side13 H1 smoke probe;
  - `julia --project=. -e 'using GaussletBases; println("load ok")'`;
  - `git diff --check`.

Required response:
  - files edited and helper/object names;
  - final-basis/H1 smoke result and timing;
  - whether complete core/shell final-basis helper stayed unchanged;
  - whether the seam is shaped for driver-stage consumption;
  - validation run and result;
  - deletion/shrinkage report:
      - what old/fallback/oracle surface became less necessary;
      - what was deleted or simplified, if anything;
      - if nothing was deleted, why not;
      - whether any new test replaces/shrinks older coverage or is genuinely
        new live-contract coverage;
      - remaining stale or duplicate surfaces to retire next.

Continue the baton loop after writing `response.060.md`. Do not stop after one
pass unless blocked by a real design decision or a failure that cannot be
resolved without manager input. Do not request UI escalation; write
`.agent_handoffs/ATTENTION.md` if you are blocked.

-- repo-manager@macmini
