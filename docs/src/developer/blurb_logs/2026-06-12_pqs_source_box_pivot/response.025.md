Status: complete.

Implemented the narrow CPBM final-basis by-center nuclear helper for PQS:

- `pqs_source_shell_final_electron_nuclear_by_center_from_boundary_block(final_basis, retained_boundary_result)`

Files changed:

- `src/cartesian_pair_block_materialization/CartesianPairBlockMaterialization.jl`
- `src/cartesian_pair_block_materialization/pqs_source_shell_final_basis.jl`
- `test/nested/cartesian_pair_block_materialization_contract_runtests.jl`

Implementation details:

- The helper requires an available `:pqs_source_shell_realization_final_basis`.
- The helper accepts the live retained-boundary input shape produced by CPBM: `PairBlockMaterializationResult` with term `:retained_source_electron_nuclear_by_center`.
- It validates:
  - retained PQS source-mode block space;
  - materialized block;
  - by-center metadata;
  - nuclear charge recorded;
  - nuclear charge not applied;
  - centers not summed;
  - uncharged by-center convention;
  - finite symmetric boundary block;
  - boundary block shape against `final_basis.boundary_source_mode_count`.
- It computes:

```text
V_final(center) = L' * V_boundary(center) * L
```

- It preserves center metadata:
  - `center_key`
  - `center_index`
  - `center_location`
  - `nuclear_charge`
  - `nuclear_charge_recorded = true`
  - `nuclear_charge_applied = false`
  - `centers_summed = false`
  - `uncharged_by_center_convention = true`

Nonclaims:

- no charge summing;
- no H1 solve;
- no Hamiltonian assembly;
- no IDA;
- no density-density;
- no RHF;
- no driver route;
- no exports/artifacts beyond exporting the new CPBM helper.

Test coverage:

- Added one compact CPBM module-contract check using a synthetic retained-source by-center nuclear block.
- The test verifies accepted term/space, final matrix transformation, center metadata preservation, charge-unapplied convention, centers-unsummed convention, and the explicit nonclaims.
- No real projected-q-shell integration test was added.
- No tmp probe was needed in this pass because pass 024 already established the real-fixture boundary-source equivalence.

Exact next blocker:

```text
:missing_pqs_final_one_electron_hamiltonian_assembly
```

The next pass should add the Hamiltonian-stage assembly that applies nuclear charges and sums separated by-center final nuclear matrices with kinetic. H1 remains unavailable until that stage exists.

Validation:

- `julia --project=. test/nested/cartesian_pair_block_materialization_contract_runtests.jl` passed.
- `julia --project=. -e 'using GaussletBases; println("load ok")'` passed.
- `git diff --check` passed.

Deletion/shrinkage report:

- No old production code, tests, metadata, or compatibility paths became fully obsolete in this pass because this was the first final-basis nuclear materialization seam.
- Nothing was deleted.
- I simplified the implementation before publishing by removing a generic adapter/blocked-result layer from the first draft; the committed working diff is scoped to the actual CPBM retained block result returned by the source helper.
- The new test is genuinely new live-contract coverage for the new helper. It does not replace older coverage yet.
- Remaining stale/duplicate surfaces to retire next: retained-source/PQS H1 scratch probes and current-route safe-term oracle vocabulary should be revisited after final one-electron Hamiltonian assembly exists.

-- repo-doer@macmini
