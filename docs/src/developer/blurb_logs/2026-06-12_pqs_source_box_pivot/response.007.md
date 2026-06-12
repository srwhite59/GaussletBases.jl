Pass 007 response: PQS source-space electron-nuclear supplied-factor block

Status: completed.

Files changed:

- `src/cartesian_pair_block_materialization/CartesianPairBlockMaterialization.jl`
- `src/cartesian_pair_block_materialization/pqs_source_safe_terms.jl`
- `test/nested/cartesian_pair_block_materialization_contract_runtests.jl`

Implemented:

- Added and exported:
  - `pqs_source_pair_electron_nuclear_by_center_block(record; coulomb_expansion, center_record, gaussian_factor_terms_1d)`
  - `pqs_source_pair_retained_electron_nuclear_by_center_block(record; ...)`
- The source helper requires a ready PQS/PQS source-pair record.
- It accepts caller-supplied term-first Gaussian factor arrays for x/y/z with
  shape `(nterms, left_axis_count, right_axis_count)`.
- It requires `length(coulomb_expansion.coefficients) == nterms`.
- It builds the raw product-source block as:

```text
sum_t (-coefficients[t]) * Gx[t] * Gy[t] * Gz[t]
```

- It returns a `PairBlockMaterializationResult` with:
  - term `:source_electron_nuclear_by_center`
  - block space `:raw_product_source_modes`
  - by-center metadata
  - nuclear charge recorded
  - nuclear charge not applied
  - centers not summed
  - uncharged by-center convention
  - no shell realization, Lowdin, IDA, Hamiltonian, driver, export, or artifact
    claim.
- The retained helper is only the requested direct wrapper:

```text
pqs_source_pair_retained_one_body_block(source_result)
```

No source-axis Gaussian factor generator, global matrix, H1 diagnostic, shell
projection, Lowdin, IDA, Hamiltonian assembly, driver path, or fallback route
was added.

Source/retained nuclear block evidence:

- Extended the existing compact PQS source-pair contract test with synthetic
  two-term Gaussian factor arrays and a two-term `CoulombGaussianExpansion`.
- The source block is checked against an independent manual two-term product
  contraction over the existing source-mode ordering.
- The test checks the live nuclear convention:
  - `nuclear_charge == 2.0`
  - `nuclear_charge_recorded == true`
  - `nuclear_charge_applied == false`
  - `centers_summed == false`
  - `uncharged_by_center_convention == true`
- The retained nuclear block is checked against:

```text
source_block[left_retained_columns, right_retained_columns]
```

- The retained wrapper is checked to produce the same block as explicit retained
  contraction.

Deletion/shrinkage report:

- No old code path became fully unnecessary yet because this pass only adds the
  supplied-factor PQS source-space nuclear block. The source-axis Gaussian
  factor generator and H1 diagnostic are still missing, so old CCPM local
  Gaussian helpers remain useful as oracle/reference surfaces for the next
  comparison.
- I did simplify one active metadata construction point: the repeated common
  raw product-source metadata in `_pqs_source_pair_product_result(...)` is now
  factored through `_pqs_source_pair_common_metadata(...)`, which the new
  nuclear helper also uses.
- No old test was deleted. The new assertions are genuinely new live-contract
  coverage for the supplied-factor source/retained nuclear block and do not
  preserve old CCPM helper vocabulary.
- Remaining stale/duplicate surfaces to retire later: CCPM
  `_pqs_pqs_source_box_*local_gaussian*` and
  `_pqs_pqs_source_box_*nuclear_attraction*` helpers should become oracle-only
  or shrink once the CPBM PQS source-axis Gaussian factor generator exists and
  has a direct equivalence check.

Validation:

- `julia --project=. test/nested/cartesian_pair_block_materialization_contract_runtests.jl` passed.
- `julia --project=. -e 'using GaussletBases; println("load ok")'` passed.
- `git diff --check` passed.

-- repo-doer@macmini
