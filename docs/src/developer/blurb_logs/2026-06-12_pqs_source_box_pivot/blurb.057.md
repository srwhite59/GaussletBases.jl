Purpose:
  Fixture-quality audit before promoting any complete core/shell PQS He RHF
  gate. Do not promote q=9 yet.

Correction:
  The previous q=9 gate-promotion direction is superseded. Large `q` is only
  meaningful if the physical box, mapping/distortion parameters, local spacing,
  and number of shells are also appropriate. The current q=5/q=7/q=9/q=11
  ladder used one surrounding shell and was primarily a route/scaling probe.

  If you already started from the earlier pass-057 blurb, stop that promotion
  work and switch to this audit. Do not add a permanent q=9 gate in this pass.

Context:
  The ladder results are strong route evidence:

  ```text
  q=5   dim 223   RHF -2.7213372828531668
  q=7   dim 561   RHF -2.810068050134403
  q=9   dim 1115  RHF -2.8499091618019303
  q=11  dim 1933  RHF -2.8559475204289022
  ```

  But this does not by itself prove that the fixture follows the physically
  intended WL/PQS convergence pattern. In particular, ask whether:

  ```text
  box radius is large enough;
  central spacing is fine enough;
  mapping/distortion parameters are appropriate;
  one surrounding shell is enough;
  increasing q alone is a valid accuracy ladder.
  ```

Exact task:
  Do a read-only/probe-only fixture audit of the q=5/q=7/q=9/q=11 complete
  core/shell PQS He ladder. No production source changes and no permanent test.

  For each q point, report:
  - current box and inner box;
  - raw source dimensions;
  - physical coordinate endpoints and approximate radius;
  - central spacing / minimum local spacing near the origin;
  - spacing range over the inner/core region if cheap;
  - mapping parameters used by `_pqs_density_test_bundle` / test bundle;
  - number of surrounding shell layers actually present;
  - final dimension;
  - H1/RHF result already measured;
  - whether the point should be interpreted as route/scaling validation or
    physical accuracy evidence.

  Also answer explicitly:

  ```text
  Do we have enough shells for a physical convergence claim?
  Is q=9 or q=11 appropriate as a permanent gate now?
  What fixture change would be needed next: more shells, different d/s,
  larger R, smaller central spacing, or some combination?
  ```

Artifacts:
  Write ignored artifacts:

  ```text
  tmp/work/pqs_complete_core_shell_fixture_quality_audit.jl
  tmp/work/pqs_complete_core_shell_fixture_quality_audit_summary.txt
  ```

Docs:
  Update `docs/src/developer/pqs_near_term_final_basis_realization_plan.md`
  or `docs/src/developer/numerical_contracts.md` only if there is already a
  concise place for this status. The note should say:

  ```text
  q ladder is route/scaling validation;
  not yet accepted physical convergence;
  one-shell fixture is not enough to claim final He accuracy;
  q=9/q=11 gate promotion is deferred pending fixture review.
  ```

Do not:
  - add or promote a q=9/q=11 permanent test;
  - add a ladder of tests;
  - modify production source;
  - run new long physics ladders;
  - add GTO, driver wiring, exports, or artifacts.

Validation:
  - run the fixture-quality audit if you create a script;
  - `julia --project=. -e 'using GaussletBases; println("load ok")'`;
  - `git diff --check`.

Required response:
  - files edited;
  - concise table of q=5/q=7/q=9/q=11 fixture geometry and physics status;
  - direct answer on whether there are enough shells;
  - direct answer on whether q=9/q=11 should be promoted now;
  - recommended next fixture direction;
  - validation run and result;
  - deletion/shrinkage report:
      - what old/fallback/oracle surface became less necessary;
      - what was deleted or simplified, if anything;
      - if nothing was deleted, why not;
      - whether any new test replaces/shrinks older coverage or is genuinely
        new live-contract coverage;
      - remaining stale or duplicate surfaces to retire next.

Continue the baton loop after writing `response.057.md`. Do not stop after one
pass unless blocked by a real design decision or a failure that cannot be
resolved without manager input. Do not request UI escalation; write
`.agent_handoffs/ATTENTION.md` if you are blocked.

-- repo-manager@macmini
