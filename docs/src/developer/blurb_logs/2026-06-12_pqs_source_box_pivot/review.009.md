Accepted.

Pass 009 made the existing source-axis transform placeholder real enough for
the next PQS step. `RawProductBoxPlan` can now carry externally materialized
per-axis coefficient transforms, the retained-unit transform contract can pass
them through from PQS unit metadata, and CPBM preflight records expose the
left/right axis transform facts as compact objects.

The default path remains unchanged: no supplied transforms means three
`:not_materialized` axis facts with no coefficient matrices. No Gaussian
factor generation, nuclear assembly, shell realization, Lowdin, IDA,
Hamiltonian, driver, export, or artifact work was added.

Manager validation:

- `julia --project=. test/nested/cartesian_raw_product_sources_contract_runtests.jl`
- `julia --project=. test/nested/cartesian_retained_unit_transform_contracts_runtests.jl`
- `julia --project=. test/nested/cartesian_pair_block_materialization_contract_runtests.jl`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git diff --check`

Next target:

Use these materialized axis transforms to project supplied parent/source-axis
Gaussian term matrices into PQS source-mode Gaussian term arrays. This should
still avoid analytic Gaussian generation and H1 diagnostics; it is the narrow
projection layer between pass 009's axis facts and pass 007's supplied-factor
nuclear block.

-- repo-manager@macmini
