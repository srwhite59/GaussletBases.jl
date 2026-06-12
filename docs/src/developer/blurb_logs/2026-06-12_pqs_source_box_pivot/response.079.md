Pass 079 complete.

Overlap transfer result:

- Removed the redundant support-space overlap transfer from `test/nested/pqs_direct_retained_final_h1_runtests.jl`.
- Removed:
  - direct `_pqs_multilayer_support_product_matrix(...)` overlap construction;
  - `pqs_complete_core_shell_final_one_body_matrix(...; term = :overlap)`;
  - the two assertions comparing transferred overlap to `fixture.final_basis.final_overlap` and identity.
- Retained the live final-overlap contract through:
  - `fixture.final_basis.final_overlap_identity_error < 1.0e-10`.

Live contract still covered:

- region-plan-backed source path remains active;
- complete core/shell final basis still checks final overlap identity;
- support kinetic and by-center nuclear construction remain covered for H1;
- final Hamiltonian symmetry/finite checks remain covered;
- ordinary H1 solve and fixed-block oracle comparison remain covered.

Validation run:

- `julia --project=. test/nested/pqs_direct_retained_final_h1_runtests.jl`
  - passed: 42 tests.
- `git diff --check`
  - passed.
- Load check was not run because no source code changed in this pass.

Deletion/shrinkage report:

- Removed one direct support-overlap helper call.
- Removed one final overlap transfer call.
- Removed two overlap-transfer assertions.
- The test is smaller and no longer preserves overlap helper vocabulary merely to recheck an identity already produced by `pqs_multilayer_complete_core_shell_final_basis(...)`.
- Remaining helper-vocabulary checks to consider later:
  - the support-space kinetic and electron-nuclear helpers are still live H1 seam/oracle machinery;
  - the explicit-box bridge comparison remains compact, but should be revisited once explicit-box source planning is quarantined further.

-- repo-doer@macmini
