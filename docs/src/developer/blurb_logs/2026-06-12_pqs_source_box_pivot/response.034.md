Pass 034 completed.

Files created:

- `src/cartesian_final_basis_realization/CartesianFinalBasisRealization.jl`
- `src/cartesian_final_basis_realization/pqs_source_shell_final_basis.jl`
- `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.034.md`

Files edited:

- `src/GaussletBases.jl`
- `src/cartesian_pair_block_materialization/CartesianPairBlockMaterialization.jl`
- `src/cartesian_pair_block_materialization/pqs_source_shell_final_basis.jl`

Functions now owned by `CartesianFinalBasisRealization`:

- `pqs_source_shell_realization_final_basis`
- `pqs_source_shell_projected_one_body_matrix`
- `pqs_source_shell_final_one_body_from_boundary_matrix`

Functions intentionally left in `CartesianPairBlockMaterialization`:

- `pqs_source_shell_final_electron_nuclear_by_center_from_boundary_block`
- `pqs_source_shell_final_one_electron_hamiltonian`

Reason: the by-center nuclear final transfer still takes CPBM-owned
`PairBlockMaterializationResult`, and the Hamiltonian helper consumes that
CPBM-local nuclear result. Moving either now would require a signature redesign,
a shared result type, or a module cycle. This pass kept the extraction
mechanical.

Compatibility choice:

- `src/GaussletBases.jl` now includes `CartesianFinalBasisRealization` after
  `CartesianRawProductSources` and before `CartesianPairBlockMaterialization`.
- CPBM imports the new module and keeps existing CPBM-qualified calls working
  through direct constant aliases:
  - `CPBM.pqs_source_shell_realization_final_basis`
  - `CPBM.pqs_source_shell_projected_one_body_matrix`
  - `CPBM.pqs_source_shell_final_one_body_from_boundary_matrix`
- No broad compatibility adapter was added.

Test changes:

- No tests were added.
- No CPBM contract test was expanded.
- Existing tests can continue to call the three moved functions through CPBM
  aliases.

Validation run:

- Focused extraction probe passed:
  - loaded `GaussletBases`;
  - checked CPBM aliases are identical to the new module functions;
  - constructed a small `q=5/L=5` synthetic final-basis object through
    `CartesianFinalBasisRealization`;
  - transferred a retained-boundary overlap matrix through the CPBM alias.
- `julia --project=. -e 'using GaussletBases; println("load ok")'` passed.
- `git diff --check` passed.
- ASCII scan of edited Julia files found no non-ASCII bytes.

Deletion/shrinkage report:

- CPBM no longer owns shell-realization final-basis construction, final overlap
  diagnostics, shell-support oracle projection, or retained-boundary
  overlap/kinetic final transfer.
- The CPBM `pqs_source_shell_final_basis.jl` file shrank to the helpers that
  still depend on CPBM result types: by-center nuclear final transfer and
  final one-electron Hamiltonian assembly.
- No test was moved or shrunk in this pass because the blurb requested a
  mechanical source extraction and allowed existing CPBM-qualified calls to
  remain through aliases.
- Remaining stale/duplicate surfaces to retire next:
  - the CPBM final-basis contract section should move or shrink into a compact
    `CartesianFinalBasisRealization` contract test;
  - the CPBM aliases should remain only as long as existing internal callers
    need compatibility;
  - the by-center nuclear final transfer can move later after the CPBM result
    type boundary is resolved.

-- repo-doer@macmini
