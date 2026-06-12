Accepted.

Pass 011 added the missing analytic source-factor step:

```text
explicit x/y/z axis layers + Coulomb exponents + center
-> gaussian_factor_matrices(layer; exponents, center)
-> source-interval slices
-> pass-010 source-mode projection
```

The implementation stays in the intended PQS source-box lane. It does not call
the old CCPM source-box nuclear wrappers, does not apply nuclear charge, and
does not assemble electron-nuclear, H1, shell realization, Lowdin, IDA, or
driver data.

The test uses real small `UniformBasisSpec` axis layers and compares the helper
against direct low-level `gaussian_factor_matrices` slices projected through the
source-axis transforms. That is live-contract coverage for the new helper,
inside the existing CPBM contract file.

Manager validation:

- `julia --project=. test/nested/cartesian_pair_block_materialization_contract_runtests.jl`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git diff --check`

Next target:

Connect these generated projected Gaussian factors to the existing
`pqs_source_pair_electron_nuclear_by_center_block(...)` supplied-factor helper.
The next pass should produce an uncharged by-center source block, and optionally
the already-existing retained source-mode contraction, without starting a global
H1 route.

-- repo-manager@macmini
