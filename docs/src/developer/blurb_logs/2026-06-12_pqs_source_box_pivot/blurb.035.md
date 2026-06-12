Purpose:
  Clean up test ownership after the `CartesianFinalBasisRealization` extraction
  before adding new PQS kernels.

Context:
  Pass 034 moved final-basis realization/operator-transfer implementation out
  of `CartesianPairBlockMaterialization` into `CartesianFinalBasisRealization`.
  The implementation boundary is now better, but the main behavioral checks for
  the moved functions still live in the large CPBM contract test around the
  synthetic `final_basis_source_box` fixture.

Exact task:
  1. Add a compact module-contract test for `CartesianFinalBasisRealization`,
     preferably:

     `test/nested/cartesian_final_basis_realization_contract_runtests.jl`

  2. Move the behavior checks for these module-owned functions out of
     `test/nested/cartesian_pair_block_materialization_contract_runtests.jl`
     and into the new compact test:

     - `pqs_source_shell_realization_final_basis`
     - `pqs_source_shell_projected_one_body_matrix`
     - `pqs_source_shell_final_one_body_from_boundary_matrix`

  3. Keep the new test focused. It should check:

     - `q=5/L=5` source dims/count `5 x 5 x 5 / 125`;
     - boundary retained count `98`;
     - final retained count `98`;
     - final overlap identity for the simple identity fixture;
     - shell-support oracle projection for one symmetric shell operator;
     - retained-boundary overlap/kinetic final transfer;
     - no IDA/RHF/driver/artifact claims only where they are critical boundary
       claims, not every report flag.

  4. Update `test/nested/runtests.jl` to include the new compact test near the
     other Cartesian module-contract tests.

  5. Shrink the CPBM contract section. CPBM may keep:

     - a tiny alias smoke showing CPBM-qualified names point to the CFBR
       functions;
     - the CPBM-owned by-center nuclear final transfer and final one-electron
       Hamiltonian checks.

     It should not duplicate the moved module's detailed final-basis behavior
     checks.

Do not:
  - add source behavior, new helper fields, IDA, RHF, density-density, direct
    retained-boundary kernels, driver wiring, or exports/artifacts;
  - expand the CPBM contract test with new helper-vocabulary assertions;
  - add a broad integration test.

Validation:
  - `julia --project=. test/nested/cartesian_final_basis_realization_contract_runtests.jl`
  - `julia --project=. test/nested/cartesian_pair_block_materialization_contract_runtests.jl`
  - `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - `git diff --check`

Required response:
  - files created/edited;
  - which checks moved from CPBM to CFBR;
  - what remains in CPBM and why;
  - validation run and result;
  - deletion/shrinkage report:
      - CPBM test lines/assertions removed or simplified;
      - whether the new test replaces old coverage rather than duplicating it;
      - remaining test bloat or oracle surfaces to retire next.

Continue the baton loop after writing `response.035.md`. Do not stop after one
pass unless blocked by a real design decision or a failure that cannot be
resolved without manager input. Do not request UI escalation; write
`.agent_handoffs/ATTENTION.md` if you are blocked.

-- repo-manager@macmini
