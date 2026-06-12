Purpose:
  Apply test diet and support-space discipline to the PQS H1 gate. The current
  gate still builds a dense support-space overlap matrix only to transfer it
  back into the final basis and recheck identity, while
  `pqs_multilayer_complete_core_shell_final_basis(...)` already computes and
  checks the final overlap.

Task:
  Audit and, if safe, remove the redundant support-overlap transfer from
  `test/nested/pqs_direct_retained_final_h1_runtests.jl`.

Expected cleanup if safe:
  - Remove the direct call to `_pqs_multilayer_support_product_matrix(...)` from
    the H1 test.
  - Remove `final_overlap = pqs_complete_core_shell_final_one_body_matrix(...)`
    for `term = :overlap`.
  - Keep the final overlap identity check through
    `fixture.final_basis.final_overlap_identity_error`.
  - Keep kinetic, by-center nuclear, final Hamiltonian, ordinary H1 solve, and
    fixed-block oracle comparison.

If you find that the overlap transfer check protects a live contract not
covered elsewhere, keep it and explain precisely what that contract is. Do not
add a replacement broad metadata test.

Guardrail:
  Dense support-space one-body helpers are allowed for this H1 seam, but should
  not be preserved in tests merely as helper vocabulary. The durable gate
  should protect physics/workflow behavior and the region-backed route
  authority, not every intermediate support operator.

Do not:
  - change source code unless the test reveals a real small cleanup need;
  - add new physics values, IDA, density-density, RHF, driver wiring, exports,
    artifacts, or fixture-rule policy;
  - add tests;
  - request UI escalation. In unattended baton mode, write
    `.agent_handoffs/ATTENTION.md` and stop if permission is genuinely needed.

Validation:
  - focused H1 gate;
  - `git diff --check`;
  - load check only if source code changes.

Report:
  - whether the overlap transfer was removed or retained;
  - validation run;
  - deletion/shrinkage report:
      - test lines/helper calls removed;
      - what live contract remains covered;
      - any remaining helper-vocabulary checks to consider later.

-- repo-manager@macmini
