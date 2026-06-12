Accepted.

This pass added the first retained PQS source-mode block contraction in the
right layer. The new helper consumes a materialized raw PQS/PQS source-space
one-body block and applies the source-mode retained selectors:

```text
O_retained = O_source[left_columns, right_columns]
```

For the current boundary selector rule this is the correct form of
`T_left' * O_source * T_right`. The implementation keeps the boundary clear:
the retained result is still `:retained_pqs_source_modes`, with no shell-row
realization, Lowdin cleanup, final shell-realized pair block, Hamiltonian, IDA,
export, or artifact claim.

The preflight change fixing the `UnitPairIndexTable` seam is also appropriate:
the plan now walks pair-operator records and unit pairs together instead of
assuming the lightweight unit-pair source is indexable like the old rich tuple
path.

The test additions are a little metadata-heavy, but acceptable for this pass
because they protect the new PQS drift boundary: retained source-mode blocks
must not silently become shell-realized/Lowdin/final blocks. Do not keep
expanding this assertion style in later passes; future tests should move toward
physics or route assembly contracts.

Manager validation:

- `julia --project=. test/nested/cartesian_pair_block_materialization_contract_runtests.jl`
- `julia --project=. test/nested/cartesian_retained_unit_transform_contracts_runtests.jl`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git diff --check`

Next target:

Wire the retained source-mode block contraction through the existing PQS
one-body selector/batch surface, still for safe one-body terms only. This makes
the retained block path the active consumer for the next source-box-first route
layer without jumping to shell realization, electron-nuclear, IDA, Hamiltonian,
or driver work.

-- repo-manager@macmini
