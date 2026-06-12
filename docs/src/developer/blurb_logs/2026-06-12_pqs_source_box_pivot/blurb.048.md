Purpose:
  Stabilize the complete `223`-function PQS H1 route as the durable compact
  one-body gate, and reduce pressure from the old boundary-shell-only H1 gate.

Context:
  Pass 047 resolved the bad H1 convention in the probe. The correct nuclear
  factor source for this PGDG/numerical-reference bundle is:

  ```text
  pgdg_intermediate.gaussian_factor_terms
  ```

  Corrected result:

  ```text
  complete PQS final dimension: 223
  complete PQS H1:             -0.48047934800387226
  same-geometry fixed oracle:  -0.48047920531279725
  exact H reference:           -0.5
  ```

  The previous `98`-function boundary-shell H1 gate was useful mechanically,
  but it is not a physical one-center acceptance basis.

Exact task:
  Update the permanent PQS H1 test coverage so the live gate protects the
  complete core/surrounding-shell route, not the boundary-shell-only physical
  non-target.

Implementation guidance:
  - Prefer updating `test/nested/pqs_direct_retained_final_h1_runtests.jl`
    rather than adding another broad test file.
  - The durable gate should use:

    ```text
    current_box = (1:7, 1:7, 1:7)
    inner_box = (2:6, 2:6, 2:6)
    raw_source_dims = (5, 5, 5)
    final dimension = 223
    nuclear factor source = pgdg_intermediate.gaussian_factor_terms
    ordinary symmetric final solve
    ```

  - Keep assertions compact:
    - core/shell counts and final dimension;
    - final overlap identity;
    - Hamiltonian finite/symmetric;
    - H1 near `-0.48047934800387226`;
    - agreement with same-geometry fixed-block oracle within a small tolerance
      if cheap;
    - no fixed-block matrix authority in the active route;
    - no `_pqs_current_route_safe_term_matrices`;
    - no generalized final solve, IDA, density-density, RHF, GTO, driver,
      export, or artifact claims.
  - Shrink or relabel the `98`-function boundary-shell-only checks so they are
    not presented as the physical H1 gate. Delete them if the complete route
    now covers the live contract.
  - Update the near-term PQS docs or numerical contracts note with the complete
    H1 result and the nuclear-factor convention. Keep it concise.

Do not:
  - add IDA weights, density-density, RHF, GTO, driver wiring, exports, or
    artifacts;
  - add a second broad H1 test alongside the old one unless it replaces/shrinks
    older coverage;
  - add metadata-vocabulary assertions beyond the compact live contract.

Validation:
  - run the updated focused PQS H1 test;
  - `julia --project=. -e 'using GaussletBases; println("load ok")'`;
  - `git diff --check`;
  - docs ASCII scan if docs are edited.

Required response:
  - files edited;
  - H1 value and oracle delta from the permanent test/probe;
  - what old boundary-shell-only test pressure was removed, shrunk, or relabeled;
  - validation run and result;
  - deletion/shrinkage report:
      - what old/fallback/oracle surface became less necessary;
      - what was deleted or simplified;
      - if nothing was deleted, why not;
      - whether any new test replaces/shrinks older coverage or is genuinely
        new live-contract coverage;
      - remaining stale or duplicate surfaces to retire next.

Continue the baton loop after writing `response.048.md`. Do not stop after one
pass unless blocked by a real design decision or a failure that cannot be
resolved without manager input. Do not request UI escalation; write
`.agent_handoffs/ATTENTION.md` if you are blocked.

-- repo-manager@macmini
