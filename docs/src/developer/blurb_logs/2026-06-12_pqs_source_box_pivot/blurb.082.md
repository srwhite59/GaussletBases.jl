Purpose:
  Implement the narrow PQS H1 assembly payload recommended by pass 081, so the
  H1 gate stops manually stitching together support kinetic, support nuclear,
  final one-body transfer, Hamiltonian assembly, and H1 solve.

Task:
  Add a small internal helper, likely near the multi-layer PQS source/final
  basis helpers, with a name like:

  ```julia
  pqs_multilayer_complete_core_shell_h1_payload(
      plan;
      final_basis,
      coulomb_expansion,
      center_records,
      axis_layers = nothing,
      gaussian_factor_terms_by_center = nothing,
      metadata = (;),
  )
  ```

  Shape can be adjusted if a nearby file has a clearer local naming pattern.

Required behavior:
  - Require an available `pqs_multilayer_shell_source_plan`.
  - Consume the complete core/shell final basis, or build it only if that is
    already a local pattern and does not duplicate test work.
  - Build support kinetic with `pqs_multilayer_support_kinetic_matrix(plan)`.
  - Build separated support electron-nuclear matrices with
    `pqs_multilayer_support_electron_nuclear_by_center_matrices(...)`.
  - Transfer kinetic and each separated nuclear matrix through
    `CartesianFinalBasisRealization.pqs_complete_core_shell_final_one_body_matrix`.
  - Assemble the final one-electron Hamiltonian with existing
    `CartesianFinalBasisRealization` helper.
  - Solve ordinary final-basis H1 with existing helper.
  - Return a compact payload containing final kinetic, final nuclear records,
    final Hamiltonian, H1 solve, and summary/nonclaim fields.

Nonclaims:
  The helper must explicitly not materialize IDA, density-density, RHF, GTO,
  driver route adoption, exports, artifacts, or fixture-rule policy.

Test update:
  Update `test/nested/pqs_direct_retained_final_h1_runtests.jl` so the active
  H1 path uses the new payload helper. Keep the small axis-layer/origin-factor
  nuclear convention comparison if it remains the only focused check for that
  path. Keep the fixed-block oracle energy comparison.

Do not:
  - add new physics values or change the accepted H1 value;
  - add IDA, density-density, RHF, GTO, driver wiring, exports, artifacts, or
    fixture-rule policy;
  - add a broad new test file;
  - assert every payload field;
  - request UI escalation. In unattended baton mode, write
    `.agent_handoffs/ATTENTION.md` and stop if permission is genuinely needed.

Validation:
  - focused H1 gate;
  - load check;
  - `git diff --check`.

Report:
  - helper name and file;
  - what test-local glue was removed;
  - validation run;
  - deletion/shrinkage report:
      - manual assembly removed from test;
      - remaining test-local oracle/convention pieces;
      - any support-space helper pressure still remaining.

-- repo-manager@macmini
