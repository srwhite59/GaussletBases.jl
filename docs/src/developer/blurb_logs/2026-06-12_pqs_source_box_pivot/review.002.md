Accepted.

The implementation added the first source-box-first retained rule in the right
place. `CartesianRawProductSources` now owns a compact
`PQSBoundaryProductModeRetainedRule`, and PQS retained-unit transform-contract
metadata now carries both the rule and a compact summary.

The live contract is clear:

```text
source modes (5, 5, 5) = 125
interior modes (3, 3, 3) = 27
boundary retained modes = 98
```

This is source-mode boundary selection. It does not materialize shell rows,
shell projection, Lowdin cleanup, coefficient maps, pair blocks, IDA data,
Hamiltonians, exports, or artifacts.

The added tests are justified live-contract coverage rather than broad metadata
coverage. They check the retained count, deterministic source-mode order,
boundary selector, nonclaims, and transform-contract attachment.

Manager validation:

- `julia --project=. test/nested/cartesian_raw_product_sources_contract_runtests.jl`
- `julia --project=. test/nested/cartesian_retained_unit_transform_contracts_runtests.jl`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git diff --check`

Next target:

Use this retained source-mode rule to build retained safe one-body blocks from
existing PQS/PQS raw source-space blocks:

```text
O_retained = T_left' * O_source * T_right
```

For the current selector rule this can be implemented as source-column
selection:

```text
O_retained = O_source[left_columns, right_columns]
```

Still do not implement shell realization, Lowdin, electron-nuclear, IDA,
Hamiltonian assembly, or drivers.

-- repo-manager@macmini
