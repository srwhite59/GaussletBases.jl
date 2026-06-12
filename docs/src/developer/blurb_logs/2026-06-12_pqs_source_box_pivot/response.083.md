Pass 083 complete.

Final payload fields after trimming:

- `object_kind`
- `status`
- `blocker`
- `final_kinetic`
- `final_nuclear_by_center`
- `final_hamiltonian`
- `h1`
- `summary`
- `metadata`

Removed from the payload return shape:

- `source_plan`
- `final_basis`
- `support_kinetic`
- `support_nuclear_by_center`

Convention comparison support matrix:

- The active H1 path still uses `pqs_multilayer_complete_core_shell_h1_payload(...)`.
- The axis-layer/origin-factor convention comparison now builds its explicit support nuclear matrix separately with `pqs_multilayer_support_electron_nuclear_by_center_matrices(...)`.
- That support matrix is not exposed through the H1 payload.

Validation run:

- `julia --project=. test/nested/pqs_direct_retained_final_h1_runtests.jl`
  - passed: 46 tests.
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - passed.
- `git diff --check`
  - passed.

Deletion/shrinkage report:

- Removed broad returned payload fields that carried consumed inputs or support-space intermediates.
- Kept final assembly outputs only: final kinetic, final by-center nuclear records, final Hamiltonian, and H1 solve.
- The H1 test still contains a small support-space convention/oracle comparison, but that is now separate from the active payload contract.
- Remaining support-space helper pressure:
  - support kinetic and nuclear helpers remain internal implementation details of the H1 payload;
  - the explicit support nuclear helper remains in the test only for the axis-layer/origin-factor convention check.

-- repo-doer@macmini
