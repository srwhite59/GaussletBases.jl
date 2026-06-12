Purpose:
  Perform the first mechanical module-boundary extraction for PQS final-basis
  realization, based on the pass 033 audit.

Context:
  PQS has crossed into a real final-basis H1 seam. The current final-basis
  realization/operator-transfer helpers live inside
  `CartesianPairBlockMaterialization`, but that module should not become the
  owner of later PQS final-basis concepts. This pass should create the narrow
  `CartesianFinalBasisRealization` module and move only the helpers that can
  live there without depending on CPBM result types.

Manager decision from review 033:
  Move now:
    - `pqs_source_shell_realization_final_basis`
    - `pqs_source_shell_projected_one_body_matrix`
    - `pqs_source_shell_final_one_body_from_boundary_matrix`

  Leave in `CartesianPairBlockMaterialization` for now:
    - `pqs_source_shell_final_electron_nuclear_by_center_from_boundary_block`
    - `pqs_source_shell_final_one_electron_hamiltonian`

  Reason:
    The by-center nuclear helper conceptually belongs with final-basis operator
    transfer, but its current input type is CPBM-owned
    `PairBlockMaterializationResult`. Moving it now would force a signature
    redesign or module cycle. Keep this pass mechanical.

Exact task:
  1. Create:

     `src/cartesian_final_basis_realization/CartesianFinalBasisRealization.jl`
     `src/cartesian_final_basis_realization/pqs_source_shell_final_basis.jl`

  2. Move the three approved functions and their private helper functions into
     the new module:

     - `pqs_source_shell_realization_final_basis`
     - `pqs_source_shell_projected_one_body_matrix`
     - `pqs_source_shell_final_one_body_from_boundary_matrix`

  3. Include the new module from `src/GaussletBases.jl` after
     `CartesianRawProductSources` and before `CartesianPairBlockMaterialization`.

  4. Make `CartesianPairBlockMaterialization` import/use the new module and keep
     existing CPBM-qualified calls working by direct aliases or reexports only.
     Do not add broad compatibility adapters.

  5. Leave the nuclear-by-center final transfer and final Hamiltonian assembly
     in CPBM. If they need validation helpers that moved, either use public
     moved functions/fields or keep a small CPBM-local validator. Do not create
     a module cycle.

  6. Update comments/docs only where they would otherwise say the moved file is
     CPBM-owned. Keep wording concise.

Do not:
  - add IDA, RHF, density-density, direct retained-boundary kernels, driver
    wiring, exports/artifacts, or GTO changes;
  - redesign result shapes or introduce structs in this pass;
  - grow `test/nested/cartesian_pair_block_materialization_contract_runtests.jl`
    with new helper-vocabulary assertions;
  - add a broad test file unless needed to preserve the moved module contract.

Test policy:
  Prefer updating existing references and running focused validation. If you add
  a new test file, it must replace or shrink the existing CPBM final-basis
  section rather than duplicate it. A CPBM alias smoke is acceptable if small.

Suggested validation:
  - focused CPBM/final-basis contract section or compact moved-module test;
  - `julia --project=. -e 'using GaussletBases; println("load ok")'`;
  - `git diff --check`.

Required response:
  - files moved/created/edited;
  - exact functions now owned by `CartesianFinalBasisRealization`;
  - exact functions intentionally left in CPBM and why;
  - compatibility choice for CPBM-qualified calls;
  - validation run and result;
  - deletion/shrinkage report:
      - what CPBM no longer owns;
      - what was deleted, moved, or simplified;
      - whether tests were moved/shrunk or why not;
      - remaining stale/duplicate surfaces to retire next.

Continue the baton loop after writing `response.034.md`. Do not stop after one
pass unless blocked by a real design decision or a failure that cannot be
resolved without manager input. Do not request UI escalation; write
`.agent_handoffs/ATTENTION.md` if you are blocked.

-- repo-manager@macmini
