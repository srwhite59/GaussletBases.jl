Purpose:
  Correct the next PQS direction before more source-plan/H1 helper code
  accumulates. The recent multi-layer PQS work is useful, but
  `pqs_multilayer_shell_source_plan(...)` is currently doing shellification-like
  geometry/coverage work itself from explicit `core_box` / `outer_box` inputs.
  That should not harden into the route architecture.

Current concern:
  The clean driver/module split is:

  ```text
  cartesian_shells / CartesianShellification
      parent geometry + nuclei + route policy
      -> owned terminal regions / core / ordered shell layers / coverage

  PQS lowering/source planning
      shellification-owned regions
      -> repeated one-cell PQS source descriptors / collapsed shell sector

  CartesianFinalBasisRealization
      core support + shell sector
      -> final orthonormal basis and final operator transfer
  ```

  The current `pqs_multilayer_shell_source_plan(...)` is a tactical bridge. It
  accepts explicit boxes, requires symmetric equal shell depth, computes layer
  boxes, support coverage, duplicate counts, and collapsed shell sector data
  internally. That is acceptable as a prototype, but it should become a
  consumer of shellification/lowering ownership data, not a private
  shellification substitute.

Task:
  Do a focused audit/design pass for this boundary.

  1. Read the current `CartesianShellification` ownership contract and the
     current `pqs_multilayer_shell_source_plan(...)` implementation.
  2. Identify exactly which facts in `pqs_multilayer_shell_source_plan(...)`
     are really shellification/geometry ownership facts, and which facts are
     PQS lowering/source-plan facts.
  3. Propose the smallest concrete intermediate object or helper boundary that
     lets PQS multi-layer source planning consume shellification/lowering
     output. It should record at least:
       - core owned region / core box;
       - ordered one-cell shell layer regions or layer boxes;
       - outer support coverage;
       - disjointness / duplicate-support checks;
       - provenance that the geometry came from shellification/lowering rather
         than ad hoc PQS box arithmetic.
  4. Say whether this can be a small mechanical implementation next pass, or
     whether existing `CartesianShellification` output is missing a necessary
     fact.

Allowed edits:
  - Prefer no source edits in this pass.
  - A concise docs update is allowed if it records the boundary clearly in
    `docs/src/developer/pqs_near_term_final_basis_realization_plan.md`.
  - Add the tracked response log.

Do not:
  - add another special case to `pqs_multilayer_shell_source_plan(...)`;
  - add or change H1/HF/RHF/IDA/density-density behavior;
  - add fixture-rule policy for `Z`, `d`, `s`, radius, q, or shell depth;
  - add exports, driver wiring, or acceptance gates;
  - add broad tests or helper-vocabulary tests;
  - request UI escalation. In unattended baton mode, write
    `.agent_handoffs/ATTENTION.md` and stop if permission is genuinely needed.

Validation:
  - If docs-only: `git diff --check` is enough.
  - If you unexpectedly touch executable code, run the smallest focused test
    that covers it plus `julia --project=. -e 'using GaussletBases; println("load ok")'`.

Report:
  - boundary recommendation;
  - exact source files/functions inspected;
  - whether a small implementation pass can follow;
  - any docs edited;
  - validation run;
  - deletion/shrinkage report:
      - what special-case responsibilities should move out of
        `pqs_multilayer_shell_source_plan(...)`;
      - what old probe/test glue would become unnecessary after that move;
      - if nothing was deleted, why this audit did not obsolete a concrete
        surface yet.

-- repo-manager@macmini
