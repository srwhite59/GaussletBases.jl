Accepted.

Pass 016 added the first repo-owned PQS source-axis transform fact builder:

```text
pqs_source_axis_transform_facts_from_pgdg_axes(...)
```

It lives in CPBM, not `CartesianRawProductSources`, preserving CRPS as a
metadata-only module. The helper accepts x/y/z PGDG intermediates or mapped
ordinary gausslet bundles, calls the low-level doside/PGDG transform kernel
directly, and returns materialized CRPS `AxisSourceTransformFact`s.

Important boundary:

The implementation does not call `_cartesian_raw_product_box_plan(...)`, old
fixed-block matrices, CCPM wrappers, shell realization, Lowdin cleanup, IDA,
Hamiltonian, driver, exports, or artifacts.

The focused test compares the new helper to the old
`_cartesian_source_box_axis_transform_plan(...)` as an oracle on the pass-015
5 x 5 x 5 fixture. That is acceptable live-contract coverage for the new
transform-builder seam. Do not add another test for the next H1 probe.

Manager validation:

- `julia --project=. test/nested/cartesian_pair_block_materialization_contract_runtests.jl`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git diff --check`

Next target:

Rerun the pass-015 5 x 5 x 5 retained-source H1 probe using
`pqs_source_axis_transform_facts_from_pgdg_axes(...)` instead of the old nested
raw-product-box helper. This should be a `tmp/work` probe and should prove the
current path no longer needs the old helper for source-axis transforms.

-- repo-manager@macmini
