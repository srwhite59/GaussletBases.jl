Implemented the narrow CPBM helper that converts a proven PQS retained-boundary one-body operator into the shell-realized final basis for overlap and kinetic.

Files changed:

- `src/cartesian_pair_block_materialization/CartesianPairBlockMaterialization.jl`
- `src/cartesian_pair_block_materialization/pqs_source_shell_final_basis.jl`
- `test/nested/cartesian_pair_block_materialization_contract_runtests.jl`

Implementation:

- Added and exported:

```julia
pqs_source_shell_final_one_body_from_boundary_matrix(
    final_basis,
    retained_boundary_result;
    term = nothing,
)
```

- The helper requires:
  - `final_basis.status == :available_pqs_shell_realization_final_basis`
  - `retained_boundary_result.object_kind == :pqs_retained_source_one_body_matrix`
  - `retained_boundary_result.matrix_materialized == true`
  - `retained_boundary_result.matrix_space == :retained_pqs_source_modes`
- It accepts only overlap and kinetic in this pass.
- It validates boundary-operator shape against `final_basis.boundary_source_mode_count`.
- It validates finite and symmetric boundary entries.
- It computes:

```text
O_final = L' * O_boundary * L
```

where `L = final_basis.lowdin_cleanup`.

- It returns the boundary operator, final operator, finite/symmetry diagnostics, and explicit provenance that the input is a retained boundary operator, not an arbitrary raw source operator.
- It keeps explicit nonclaims false for electron-nuclear, charge summing, H1, Hamiltonian assembly, IDA, density-density, RHF, driver route, exports, and artifacts.

Test coverage:

- Extended the compact CPBM contract test with synthetic retained-boundary overlap and kinetic result objects in the same shape as `pqs_retained_source_one_body_matrix(...)`.
- Verified final overlap/kinetic transformation through the identity final-basis fixture.
- Verified retained-boundary provenance and nonclaim fields.

No real-fixture probe was added in this pass. Pass 022 already proved on the real projected-q-shell fixture that retained-source overlap/kinetic are the shell-projected boundary operators to roundoff; this pass just adds the production transformation seam.

Exact next blocker:

```text
:missing_pqs_shell_boundary_electron_nuclear_operator_source
```

H1 is still not ready. Overlap and kinetic can now move from retained boundary operators to final shell-realized one-body matrices, but the by-center electron-nuclear boundary/source operator is still missing.

Validation:

- `julia --project=. test/nested/cartesian_pair_block_materialization_contract_runtests.jl`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git diff --check`

Deletion/shrinkage report:

- No old code was deleted in this pass. The final-from-boundary seam is now live, but H1 still needs electron-nuclear boundary/source materialization before old retained-source H1 probes can be retired.
- No compatibility adapter was added for arbitrary matrix inputs; the helper is deliberately tied to `pqs_retained_source_one_body_matrix(...)` results.
- The added test is live-contract coverage for the new seam and uses a compact synthetic retained-boundary result instead of broad old integration scaffolding.
- Remaining stale surfaces to retire next: retained-source overlap/kinetic readiness probes can shrink after a final overlap/kinetic probe consumes this helper; retained-source H1 probes and current-route support-local safe-term oracle checks remain until electron-nuclear and final H1 are available.

-- repo-doer@macmini
