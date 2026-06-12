Purpose:
  Promote the direct-retained PQS final-basis H1 seam into one compact durable
  gate, while shrinking old slow integration coverage that preserves
  `_pqs_current_route_safe_term_matrices(...)` helper vocabulary.

Context:
  Pass 038 proved the explicit final-basis H1 seam in an ignored probe:

  - `5 x 5 x 5` source modes, source count `125`;
  - boundary/final retained count `98 / 98`;
  - direct retained overlap/kinetic/by-center nuclear all used;
  - active retained path did not materialize raw source one-body blocks;
  - final Hamiltonian matched shell-support oracle to about `1.3e-15`;
  - H1 lowest eigenvalue matched oracle to about `3.8e-16`;
  - ordinary symmetric eigensolve, not generalized overlap solve.

Exact task:
  1. Add one compact durable test for the direct-retained PQS final-basis H1
     seam, for example:

     `test/nested/pqs_direct_retained_final_h1_runtests.jl`

     Base it on `tmp/work/pqs_final_basis_h1_probe.jl`, but keep the permanent
     test compact and readable.

  2. The new gate should check only the live workflow contract:

     - source dims/count `5 x 5 x 5 / 125`;
     - boundary/final retained count `98 / 98`;
     - final overlap identity error is small;
     - retained overlap/kinetic/nuclear direct-boundary flags are true;
     - active retained overlap/kinetic/nuclear did not materialize raw source
       one-body blocks;
     - final Hamiltonian is finite/symmetric;
     - final Hamiltonian and H1 eigenvalue match the shell-support oracle;
     - ordinary symmetric eigensolve is used, not a generalized overlap solve;
     - no IDA, density-density, RHF, driver, export, or artifact claim.

  3. Add the new test to `test/nested/runtests.jl` only if it remains compact
     and reasonably fast. If it is too slow for the default nested runner,
     explain and place it in the appropriate integration runner instead.

  4. Shrink the old slow integration section in
     `test/nested/bond_aligned_diatomic_high_order_recipe_opt_in_source_construction_integration_runtests.jl`
     around these calls:

     - `CCPM._pqs_current_route_safe_term_matrices(...)`
     - `CCPM._pqs_current_route_safe_term_authority_comparison(...)`

     Remove or reduce the detailed assertions over whole-route safe-term
     matrices, oracle matrices, per-term helper fields, and authority comparison
     metadata that are now helper-vocabulary pressure. Keep only a minimal
     historical/private-diagnostic smoke if needed.

  5. Do not delete the private CCPM helpers in this pass. They may remain
     oracle/debug code until a source deletion pass has a clear caller-driven
     deletion condition.

Do not:
  - add IDA, density-density, RHF, driver wiring, exports/artifacts, or GTO
    changes;
  - expand CPBM or CFBR contract tests with duplicate H1 coverage;
  - copy the full probe artifact into tracked docs;
  - preserve broad old helper-vocabulary assertions just because they already
    exist.

Validation:
  - run the new H1 gate directly;
  - run the edited integration test if feasible; if it is too slow, report that
    explicitly and run the narrowest available syntax/load validation;
  - `julia --project=. -e 'using GaussletBases; println("load ok")'`;
  - `git diff --check`.

Required response:
  - files created/edited;
  - H1 gate result and runtime;
  - exact old integration assertions/blocks removed or retained;
  - validation run and result;
  - deletion/shrinkage report:
      - net test-line/assertion change if easy to report;
      - whether the H1 gate replaces old coverage rather than duplicating it;
      - private CCPM helper surfaces left for future deletion/quarantine.

Continue the baton loop after writing `response.039.md`. Do not stop after one
pass unless blocked by a real design decision or a failure that cannot be
resolved without manager input. Do not request UI escalation; write
`.agent_handoffs/ATTENTION.md` if you are blocked.

-- repo-manager@macmini
