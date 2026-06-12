Purpose:
  Define and probe the smallest complete one-center PQS/final-basis H1 target,
  replacing the boundary-shell-only physical fixture before any IDA/RHF work.

Manager decision:
  After pass 040, choose option 2: first define a complete one-center
  PQS/final basis so H/He+ H1 is physically meaningful. Do not proceed to final
  IDA weights on the boundary-shell-only route as the next scientific target.

Context:
  Pass 040 showed the direct-retained final-basis route is mechanically correct
  but physically incomplete:

  - the active route used only the PQS boundary shell, retained/final count 98;
  - H / Z=1 H1 was `-0.08171962129085239` vs exact `-0.5`;
  - route agreement with the shell-support oracle was roundoff;
  - therefore the problem is basis completeness, not one-body route mechanics.

  The likely missing physical ingredient is the complete one-center final basis,
  not IDA. The preferred first physical target is a `5 x 5 x 5` inner core plus
  one surrounding shell. A single all-source-mode `5 x 5 x 5` retained rule can
  be used as a diagnostic fallback if the core-plus-shell route is not yet
  expressible, but it is not the preferred complete one-center target.

Exact task:
  1. Audit the current PQS/source/final-basis surfaces and identify the
     smallest complete one-center final-basis route that can be tested without
     a broad framework.

     Consider explicitly, in this priority order:
     - a `5 x 5 x 5` inner core plus one surrounding shell;
     - the direct-core plus boundary-shell retained/final-basis split needed to
       represent that core-plus-shell fixture;
     - an all-source-mode retained rule over one `5 x 5 x 5` source box
       (`125` retained modes) only as a diagnostic fallback;
     - whether current shell-realization/Lowdin inputs can represent the
       complete core-plus-shell final basis.

  2. If an existing route can express the complete one-center final basis,
     create an ignored `tmp/work` probe and run H / Z=1 H1 through the same
     direct-retained final-basis path:

     source/final basis -> direct retained overlap/kinetic/nuclear -> final
     one-electron Hamiltonian -> ordinary symmetric H1 solve.

     Compare to exact `-0.5` and to a shell-support oracle if available.

  3. If existing production objects cannot express it cleanly, do not invent a
     broad route. Report the exact missing object/surface, for example:

     - missing all-source-mode retained rule;
     - missing combined direct-core + boundary-shell final-basis realization;
     - missing route-owned shell projection/Lowdin for complete one-center
       source modes;
     - missing direct retained operator blocks for a non-boundary retained rule.

  4. Preserve route boundaries:

     - no `_pqs_current_route_safe_term_matrices(...)` active path;
     - no generalized-overlap final solve;
     - no IDA, density-density, RHF, driver, GTO, export, or artifact work.

  5. Write a concise ignored artifact, for example:

     `tmp/work/pqs_complete_one_center_h1_probe_summary.txt`

     Include:
     - selected complete-basis strategy;
     - retained/final count;
     - physical extent;
     - H1 value and exact-reference error if runnable;
     - oracle delta if available;
     - direct-path flags;
     - exact blocker if not runnable;
     - coarse timings if cheap.

Do not:
  - add a permanent test in this pass;
  - add source behavior unless a tiny, obviously scoped helper is required and
    the response justifies why a probe-only construction was not enough;
  - add IDA, density-density, RHF, driver wiring, exports/artifacts, or GTO
    changes;
  - preserve the 98-function boundary-shell H1 as a physical acceptance target;
  - treat the all-125 single-source-box case as the preferred target if the
    `5 x 5 x 5` inner core plus one-shell route is reachable.

Validation:
  - run the probe or audit script;
  - `julia --project=. -e 'using GaussletBases; println("load ok")'`;
  - `git diff --check`.

Required response:
  - complete-basis strategy selected or exact blocker, with explicit discussion
    of the preferred `5 x 5 x 5` inner core plus one-shell target;
  - probe artifact path if runnable;
  - retained/final dimension and physical extent;
  - H1 value, exact-reference error, and oracle delta if available;
  - direct-path flags and nonclaim flags;
  - validation run and result;
  - deletion/shrinkage report:
      - whether anything was deleted or shrunk;
      - if no deletion, why this was probe/audit-only;
      - exact next implementation target for complete one-center PQS H1.

Continue the baton loop after writing `response.041.md`. Do not stop after one
pass unless blocked by a real design decision or a failure that cannot be
resolved without manager input. Do not request UI escalation; write
`.agent_handoffs/ATTENTION.md` if you are blocked.

-- repo-manager@macmini
