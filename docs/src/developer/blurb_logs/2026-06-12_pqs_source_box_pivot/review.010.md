Accepted.

Pass 010 added the narrow projection layer we needed:

```text
G_source[t] = C_left' * G_support[t] * C_right
```

It consumes materialized left/right `AxisSourceTransformFact`s from the CPBM
preflight record and supplied support-space Gaussian term matrices. It does not
generate analytic Gaussian factors, use CCPM wrappers, apply nuclear charge,
build nuclear blocks, or adopt a driver route.

The test stays in the existing PQS source-pair contract file and checks
projection against a direct manual matrix product, including a not-materialized
axis-transform failure. That is live-contract coverage for the new seam rather
than a new broad test surface.

Manager validation:

- `julia --project=. test/nested/cartesian_pair_block_materialization_contract_runtests.jl`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git diff --check`

Next target:

Generate the supplied support-space Gaussian term matrices from explicit
axis-layer objects, Coulomb exponents, and a center coordinate, then feed them
through the pass-010 projection helper. This should still avoid CCPM physical
nuclear wrappers and should not assemble an H1 diagnostic yet.

-- repo-manager@macmini
