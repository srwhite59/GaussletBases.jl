Review 079: accepted.

The pass removes the redundant support-space overlap transfer from the PQS H1
gate. The test no longer calls `_pqs_multilayer_support_product_matrix(...)`
for overlap, no longer transfers that overlap through
`pqs_complete_core_shell_final_one_body_matrix(...; term = :overlap)`, and no
longer asserts the transferred overlap against identity.

The live contract remains covered by `fixture.final_basis.final_overlap_identity_error`
plus the active region-plan-backed source path, kinetic/nuclear final transfer,
Hamiltonian symmetry/finite checks, ordinary H1 solve, and fixed-block oracle
comparison.

Manager validation:

- `julia --project=. -e 't = @elapsed include("test/nested/pqs_direct_retained_final_h1_runtests.jl"); println("elapsed_s=", t)'`
  passed, 42/42, elapsed about 6.22s.
- `git diff --check` passed.

Deletion/shrinkage:

- removed one private support-overlap helper call;
- removed one final overlap transfer call;
- removed two overlap-transfer assertions;
- no replacement test was added.

Remaining pressure:

- support kinetic and by-center nuclear helpers remain live H1 seam machinery;
- explicit-box bridge comparison remains compact and should be revisited after
  explicit-box source planning is quarantined further.

-- repo-manager@macmini
