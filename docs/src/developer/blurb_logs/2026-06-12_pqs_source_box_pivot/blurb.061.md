Purpose:
  Run the WL-style pre-RHF density/J diagnostic for the new multi-layer PQS
  side13 final basis. Do not run RHF yet.

Context:
  Pass 060 implemented `pqs_multilayer_shell_source_plan(...)` and the side13
  final-basis/H1 smoke:

  ```text
  parent count: 13
  mapping:      AsinhMapping(c = 0.1, s = 1.0, tail = 10)
  core:         (4:10)^3
  shell layers: 3
  final dim:    1549
  Z=2 H1:      -1.975561823201342
  ```

  This is comparable in purpose to the WL side13 H1/J readiness checks. Before
  RHF, verify the density interaction and self-Coulomb convention.

Exact task:
  Build an ignored `tmp/work` probe for the multi-layer PQS side13 pre-RHF J
  diagnostic.

  Use:
  - `pqs_multilayer_shell_source_plan(...)`;
  - the existing complete core/shell final-basis helper;
  - the existing final one-body/H1 helpers;
  - the existing pre-final density-interaction seam if it can consume the
    multi-layer plan's support states and weights.

  Compute:
  - Z=2 H1 lowest orbital energy;
  - support raw pair numerator and support weights for the multi-layer support;
  - pre-final density interaction;
  - H1 self-Coulomb J using
    `c_prefinal = combined_lowdin_cleanup * c_final`;
  - comparison to the hydrogenic reference `5Z/8 = 1.25`;
  - comparison to the WL side13 diagnostic if easy:
    `H1 ~= -1.9748150892830352`, `J ~= 1.2158294767735702`.

Artifacts:
  Write ignored artifacts:

  ```text
  tmp/work/pqs_multilayer_shell_side13_j_probe.jl
  tmp/work/pqs_multilayer_shell_side13_j_probe_summary.txt
  ```

Decision rule:
  If H1 and J are coherent and the density interaction is finite/symmetric with
  positive pre-final weights, report that RHF is a reasonable next probe.
  If the J diagnostic is far off or the density gauge fails, do not run RHF;
  report the exact blocker.

Do not:
  - run RHF in this pass;
  - add permanent tests;
  - add GTO, driver wiring, exports, or artifacts;
  - use signed final weights or raw no-division density;
  - use fixed-block pair data as active authority;
  - turn the provisional spacing rule into a new study.

Docs:
  Update the PQS near-term plan only if the J diagnostic yields a concise
  status note. Keep it non-acceptance.

Validation:
  - run the J probe;
  - if source changes, run the focused H1 smoke/probe or relevant compact test;
  - `julia --project=. -e 'using GaussletBases; println("load ok")'`;
  - `git diff --check`.

Required response:
  - files edited;
  - H1 and J diagnostic numbers;
  - density-interaction weight and symmetry diagnostics;
  - comparison to WL side13 if available;
  - whether RHF is a reasonable next probe or blocked;
  - validation run and result;
  - deletion/shrinkage report:
      - what old/fallback/oracle surface became less necessary;
      - what was deleted or simplified, if anything;
      - if nothing was deleted, why not;
      - whether any new test replaces/shrinks older coverage or is genuinely
        new live-contract coverage;
      - remaining stale or duplicate surfaces to retire next.

Continue the baton loop after writing `response.061.md`. Do not stop after one
pass unless blocked by a real design decision or a failure that cannot be
resolved without manager input. Do not request UI escalation; write
`.agent_handoffs/ATTENTION.md` if you are blocked.

-- repo-manager@macmini
