Implemented the first CPBM-owned shell-projected one-body operator seam for PQS.

Files changed:

- `src/cartesian_pair_block_materialization/CartesianPairBlockMaterialization.jl`
- `src/cartesian_pair_block_materialization/pqs_source_shell_final_basis.jl`
- `test/nested/cartesian_pair_block_materialization_contract_runtests.jl`

Implementation:

- Added and exported:

```julia
pqs_source_shell_projected_one_body_matrix(final_basis, shell_operator; term)
```

- The helper requires an available `:pqs_source_shell_realization_final_basis`.
- It accepts a caller-supplied real shell-support operator matrix.
- It validates shape and finite entries.
- For symmetric one-body terms, it checks shell-operator symmetry and reports the symmetry error.
- It computes:

```text
O_boundary = P' * O_shell_support * P
O_final = R' * O_shell_support * R
```

where `P = final_basis.shell_projection` and `R = final_basis.final_shell_coefficients`.

- It cross-checks the equivalent Lowdin form:

```text
O_final == L' * O_boundary * L
```

where `L = final_basis.lowdin_cleanup`.

- It returns the boundary operator, final operator, finite/symmetry/cross-check diagnostics, and explicit nonclaims:
  - no H1 solve;
  - no charge summing;
  - no IDA;
  - no density-density;
  - no RHF;
  - no Hamiltonian/driver/export/artifact materialization.

Test coverage:

- Extended the compact CPBM contract test with a supplied symmetric shell-support matrix.
- The test verifies:
  - materialized projected one-body status;
  - boundary and final operators have the expected values for the identity final-basis fixture;
  - Lowdin/boundary cross-check error is zero;
  - H1/charge/IDA/RHF/driver/artifact flags remain false;
  - no generated shell-support operator or current-route safe-term authority is claimed.

No real-fixture probe was added in this pass. Pass 020 already validated the real projected-q-shell final-basis object against the old shell plan; this pass only adds the algebraic projection seam for a caller-supplied shell-support operator.

No old route authority was adopted:

- No call to `_pqs_current_route_safe_term_matrices(...)`.
- No analytic kinetic/nuclear shell-support operator generation.
- No H1 solve.
- No old fixed-block matrix authority.
- No direct/full-parent fallback.

Exact next blocker:

```text
:missing_pqs_shell_support_one_body_operator_source
```

The algebraic projection seam now exists. The remaining constructive blocker is producing shell-support one-body operator matrices for overlap, kinetic, and by-center electron-nuclear terms through the intended PQS route, then feeding them through this projection seam.

Validation:

- `julia --project=. test/nested/cartesian_pair_block_materialization_contract_runtests.jl`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git diff --check`

Deletion/shrinkage report:

- No old code was deleted in this pass. The old source-box bridge/readiness and oracle surfaces remain useful until shell-support operator sources and a final PQS H1 probe exist.
- No compatibility adapter was added.
- The added test is new live-contract coverage for the algebraic projection seam; it does not preserve old current-route implementation vocabulary.
- Remaining stale surfaces to retire next: source-box bridge tests that only assert missing shell realization, old retained-source H1 probes, and current-route support-local safe-term oracle checks should shrink once shell-support one-body sources and final PQS H1 are available.

-- repo-doer@macmini
