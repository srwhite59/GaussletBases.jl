Pass 010 response: PQS source-mode Gaussian factor projection

Status: implemented.

Files changed:

- `src/cartesian_pair_block_materialization/CartesianPairBlockMaterialization.jl`
- `src/cartesian_pair_block_materialization/pqs_source_safe_terms.jl`
- `test/nested/cartesian_pair_block_materialization_contract_runtests.jl`

What changed:

- Added and exported:

```julia
pqs_source_pair_gaussian_factor_terms_1d(record; gaussian_factor_terms_axis)
```

- The helper requires a ready PQS/PQS source-pair preflight record with
  materialized left/right `AxisSourceTransformFact`s.
- It accepts supplied x/y/z term-first Gaussian matrices in source-support
  coordinates.
- It validates, per axis:
  - materialized left/right transform facts;
  - transform axes and source-mode dimensions;
  - real coefficient matrices;
  - source-support term array shape and common term count.
- It projects each term as:

```julia
G_source[t] = C_left' * G_support[t] * C_right
```

- It returns a named tuple `(x = ..., y = ..., z = ...)` with term-first
  arrays shaped `(nterms, left_source_mode_dim, right_source_mode_dim)`.
- It does not generate analytic Gaussian factors, call CCPM wrappers, apply
  nuclear charge, build nuclear blocks, materialize shell realization, run
  Lowdin cleanup, build IDA data, assemble Hamiltonians, adopt driver routes,
  export artifacts, or add fallback paths.

Projection evidence:

- The existing PQS source-pair contract test now uses synthetic support-space
  x/y/z Gaussian term arrays and synthetic source-axis transforms.
- It checks the projected arrays against manual
  `transpose(C_left) * G * C_right` for each axis.
- It checks projected shapes:
  - x: `(2, 3, 5)`
  - y: `(2, 3, 4)`
  - z: `(2, 3, 3)`
- It checks that a PQS source pair with default not-materialized transform
  facts throws instead of silently projecting.

Validation:

- `julia --project=. test/nested/cartesian_pair_block_materialization_contract_runtests.jl` passed.
- `julia --project=. -e 'using GaussletBases; println("load ok")'` passed.
- `git diff --check` passed.

Deletion/shrinkage report:

- No old code, test, metadata, or compatibility path became unnecessary yet.
  This pass adds the missing projection layer between pass 009's transform
  facts and pass 007's supplied-factor nuclear block, but it does not replace a
  live older CPBM implementation.
- Nothing was deleted or simplified.
- The new test coverage is live-contract coverage for the newly active
  projection helper. It extends the existing CPBM source-pair contract file and
  does not add a new test file or physics fixture.
- Remaining stale/duplicate surfaces to retire next:
  - CCPM centered local-Gaussian source-box wrappers remain oracle/reference
    only until CPBM owns the analytic source-axis Gaussian term generation.
  - A next narrow pass can generate source-support Gaussian term matrices from
    explicit axis layers/centers/exponents and feed this projection helper.

-- repo-doer@macmini
