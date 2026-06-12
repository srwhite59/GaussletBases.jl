Pass 082 complete.

Helper added:

- `pqs_multilayer_complete_core_shell_h1_payload(...)`
- File: `src/pqs_multilayer_shell_source_plan.jl`

Helper behavior:

- Requires an available `pqs_multilayer_shell_source_plan`.
- Requires an available complete core/shell final basis.
- Builds support kinetic with `pqs_multilayer_support_kinetic_matrix(plan)`.
- Builds separated support electron-nuclear matrices with `pqs_multilayer_support_electron_nuclear_by_center_matrices(...)`.
- Transfers kinetic and each separated nuclear matrix through `CartesianFinalBasisRealization.pqs_complete_core_shell_final_one_body_matrix`.
- Assembles the final one-electron Hamiltonian with `CartesianFinalBasisRealization.pqs_complete_core_shell_final_one_electron_hamiltonian`.
- Solves ordinary final-basis H1 with `CartesianFinalBasisRealization.pqs_complete_core_shell_final_h1_solve`.
- Returns final kinetic, final by-center nuclear records, final Hamiltonian, H1 solve, and a compact summary/nonclaim payload.

Test-local glue removed:

- The active H1 test no longer manually:
  - builds support kinetic;
  - builds explicit-origin support nuclear for the active path;
  - transfers active support kinetic to final basis;
  - transfers active support nuclear to final basis;
  - assembles the active H1 Hamiltonian;
  - calls the active H1 solve directly.
- `PQSH1CFBR` was removed from the test because active final assembly is now behind the payload helper.

Remaining test-local oracle/convention pieces:

- The axis-layer/origin-factor nuclear comparison remains as a focused convention check.
- The fixed-block oracle energy comparison remains.
- Compact region-plan-vs-explicit-box bridge checks remain.

Validation run:

- `julia --project=. test/nested/pqs_direct_retained_final_h1_runtests.jl`
  - passed: 46 tests.
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - passed.
- `git diff --check`
  - passed.

Deletion/shrinkage report:

- Manual H1 assembly was removed from the test and replaced by one payload call.
- The test still has a small convention/oracle section for axis-layer nuclear factors because that remains the focused check for the origin-factor path.
- Support-space helper pressure remaining:
  - support kinetic and support nuclear are still used inside the H1 payload;
  - these helpers remain scoped to the H1 seam and are not promoted as scalable PQS operator algorithms.
- Net source lines increased because the new payload is now the module-owned seam. The test’s active route logic is simpler and no longer owns final H1 assembly.

-- repo-doer@macmini
