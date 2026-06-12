Purpose:
  Align the next PQS He fixture with the White-Lindsey He test philosophy before
  running more physics probes or promoting a gate.

Context:
  The q=5/q=7/q=9/q=11 PQS ladder established route/scaling behavior, but every
  point used one surrounding shell and fixed `R ~= 8`. That is not enough for a
  physical convergence claim. The next tests should be similar in spirit to the
  WL He tests: box size, local spacing/distortion, and shell count must move
  together.

Exact task:
  Do a read-only/probe-only fixture-design audit. No production source changes
  and no permanent tests.

  Compare the current PQS complete core/shell fixture knobs against the WL He
  fixture knobs recorded in:

  ```text
  docs/src/developer/numerical_contracts.md
  docs/src/developer/pqs_near_term_final_basis_realization_plan.md
  test/nested/cartesian_wl_gausslet_he_atom_acceptance_runtests.jl
  ```

  Answer:

  1. What are the WL-side physical fixture controls for He?
     Include `R`/box extent, spacing `d`, distortion/mapping parameter `s`,
     `ns`/shell depth, and when GTO supplement data enters.

  2. What are the corresponding PQS controls?
     Include parent axis count, `AsinhMapping(a, s, tail_spacing)`, current box,
     inner/core box, raw source dimensions, surrounding shell-layer count, and
     final dimension.

  3. Does the current complete core/shell PQS construction support more than
     one surrounding shell by changing `current_box`/`inner_box`/raw source
     dimensions, or is there a missing route ingredient?

  4. What is the first WL-aligned PQS He fixture to try next?
     It should keep an adequate box, use reviewed central spacing/distortion,
     and include enough shell layers to be a real fixture-quality probe.

Optional cheap probe:
  If it is easy and does not require production source changes, run a
  final-basis-only or H1-only smoke probe for a multi-shell PQS fixture to
  determine whether the route accepts it. Do not run a long RHF job in this
  pass.

Artifacts:
  Write ignored artifacts:

  ```text
  tmp/work/pqs_wl_aligned_fixture_design_audit.jl
  tmp/work/pqs_wl_aligned_fixture_design_audit_summary.txt
  ```

Docs:
  Update the PQS near-term plan only if the result is concise and useful. The
  docs should say the next PQS physical tests should be WL-aligned fixture
  probes, not q-only scaling probes.

Do not:
  - promote q=9 or q=11 as a gate;
  - add permanent tests;
  - run a new long RHF ladder;
  - add GTO, driver wiring, exports, or artifacts;
  - make broad source refactors.

Validation:
  - run the audit/probe if you create one;
  - `julia --project=. -e 'using GaussletBases; println("load ok")'`;
  - `git diff --check`.

Required response:
  - files edited;
  - WL fixture controls vs PQS fixture controls;
  - whether multi-shell PQS is currently supported or blocked;
  - recommended first WL-aligned PQS He fixture;
  - any cheap probe result and timing;
  - validation run and result;
  - deletion/shrinkage report:
      - what old/fallback/oracle surface became less necessary;
      - what was deleted or simplified, if anything;
      - if nothing was deleted, why not;
      - whether any new test replaces/shrinks older coverage or is genuinely
        new live-contract coverage;
      - remaining stale or duplicate surfaces to retire next.

Continue the baton loop after writing `response.058.md`. Do not stop after one
pass unless blocked by a real design decision or a failure that cannot be
resolved without manager input. Do not request UI escalation; write
`.agent_handoffs/ATTENTION.md` if you are blocked.

-- repo-manager@macmini
