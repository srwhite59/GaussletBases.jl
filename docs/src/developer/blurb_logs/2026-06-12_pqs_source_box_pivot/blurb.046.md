Purpose:
  Audit and, only if small, fix the one-body operator convention for the
  complete `223`-function PQS H1 route.

Context:
  Pass 045 materialized the complete core/surrounding-shell H1 route, but the
  Z=1 hydrogen energy was:

  ```text
  H1 = -2.0638461028784776
  exact H reference = -0.5
  ```

  That is below the variational lower bound for a correct Z=1 one-electron
  Hamiltonian. The center charge in the probe was `1.0`, so this is not just a
  mislabeled Z=2 run. Treat the H1 value as a convention/fixture blocker, not a
  physics result.

Exact task:
  Run a focused one-body convention audit for the complete core/shell H1 route.
  Do not add IDA/RHF/GTO or new route features.

Required checks:
  1. Recompute the same support-space overlap, kinetic, and uncharged nuclear
     matrices used by the probe.
  2. Solve the full support-space generalized eigenproblem

     ```text
     H_support c = E S_support c
     ```

     on the combined core+shell support rows, if numerically possible. If this
     already gives an energy below `-0.5`, the operator matrices are wrong
     before final-basis transfer.
  3. For the final-basis lowest vector, report:

     ```text
     kinetic expectation
     uncharged nuclear expectation
     charged nuclear expectation for Z=1
     total
     ```

  4. Verify the center charge, nuclear sign convention, and charge-application
     stage.
  5. Compare the overlap/kinetic/nuclear construction against an existing
     trusted one-body oracle for the same geometry if one can be reached
     without making old fixed-block matrices the active route. If mapping or
     representation mismatch prevents comparison, record that blocker
     explicitly.
  6. Check whether the issue is a small convention bug, such as:
     - applying charge twice;
     - wrong sign in uncharged by-center nuclear input;
     - using raw Gaussian factor matrices where projected/operator-normalized
       factors are required;
     - inconsistent support-row/operator metric convention.

Implementation rule:
  - If the bug is small and local, fix it.
  - If the audit points to a larger missing operator convention, do not paper it
    over. Report the exact blocker and leave the H1 route marked non-acceptance.

Do not:
  - use `_pqs_current_route_safe_term_matrices(...)` as the active path;
  - use old fixed-block matrices as active authority;
  - use a generalized-overlap final solve as the final-basis acceptance route;
  - add IDA weights, density-density, RHF, GTO, driver wiring, exports, or
    artifacts;
  - add broad tests or expand CPBM contract tests.

Probe:
  Update or create ignored artifacts:

  ```text
  tmp/work/pqs_complete_core_shell_h1_convention_probe.jl
  tmp/work/pqs_complete_core_shell_h1_convention_probe_summary.txt
  ```

  Include:
  - support generalized H1 value or blocker;
  - final H1 value;
  - kinetic/nuclear decomposition;
  - charge/sign convention facts;
  - oracle comparison or exact oracle-blocker;
  - whether a source fix was made;
  - nonclaim flags.

Validation:
  - focused convention probe;
  - if source changed, rerun the H1 probe;
  - `julia --project=. -e 'using GaussletBases; println("load ok")'`;
  - `git diff --check`.

Required response:
  - files edited and exact function/object names changed;
  - audit findings and whether a fix was made;
  - H1/decomposition numbers after the audit/fix;
  - validation run and result;
  - deletion/shrinkage report:
      - what old/fallback/oracle surface became less necessary;
      - what was deleted or simplified, if anything;
      - if nothing was deleted, why not;
      - whether any new test replaces/shrinks older coverage or is genuinely
        new live-contract coverage;
      - remaining stale or duplicate surfaces to retire next.

Continue the baton loop after writing `response.046.md`. Do not stop after one
pass unless blocked by a real design decision or a failure that cannot be
resolved without manager input. Do not request UI escalation; write
`.agent_handoffs/ATTENTION.md` if you are blocked.

-- repo-manager@macmini
