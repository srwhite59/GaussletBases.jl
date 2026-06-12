Review 033: accepted.

The audit stayed inside the requested boundary: no source, test, export, or
route behavior changed. It identified the useful module seam without trying to
solve every result-shape issue at once.

Manager decision for the next pass:

- Move now:
  - `pqs_source_shell_realization_final_basis`
  - `pqs_source_shell_projected_one_body_matrix`
  - `pqs_source_shell_final_one_body_from_boundary_matrix`
- Leave in `CartesianPairBlockMaterialization` for now:
  - `pqs_source_shell_final_electron_nuclear_by_center_from_boundary_block`
  - `pqs_source_shell_final_one_electron_hamiltonian`

Reason: the nuclear-by-center helper is conceptually final-basis operator
transfer, but its current input type is the CPBM-owned
`PairBlockMaterializationResult`. Moving it in the first pass would either
create an avoidable module cycle or force a signature redesign. Keep the first
move mechanical and extract only the pieces that can live cleanly above CPBM.

The proposed module name and include order are accepted:

```text
src/cartesian_final_basis_realization/CartesianFinalBasisRealization.jl
src/cartesian_final_basis_realization/pqs_source_shell_final_basis.jl
```

Include the new module after `CartesianRawProductSources` and before
`CartesianPairBlockMaterialization`. CPBM may keep compatibility aliases for
the moved functions, but should not grow broad adapter layers.

Next pass should perform that mechanical extraction and update the minimum
necessary tests/references. Do not add IDA, RHF, density-density, direct
retained kernels, result-type redesign, or a new broad test surface in the same
pass.

-- repo-manager@macmini
