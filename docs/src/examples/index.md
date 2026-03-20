# Examples

This page is the visible examples entry point for `GaussletBases.jl`.

The package has many runnable examples, but the useful way to approach them is
as a small curated sequence rather than as a raw directory listing.

From the repository root, each example runs as:

```bash
julia --project=. examples/NAME.jl
```

## Start with these

If you are new to the package, begin with:

1. `01_first_gausslet.jl`
2. `02_radial_basis.jl`
3. `03_radial_operators.jl`
4. `04_hydrogen_ground_state.jl`

That sequence is still the clearest public entry path.

## Then choose a branch

For the mature radial and atomic line, continue with:

- `15_atomic_hydrogen_ylm.jl`
- `16_atomic_ida_ingredients.jl`
- `19_atomic_ida_direct.jl`
- `20_atomic_ida_exchange.jl`
- `21_atomic_ida_fock.jl`
- `22_atomic_ida_uhf.jl`

For the newer ordinary Cartesian line, continue with:

- `23_cartesian_hydrogen_coulomb_expansion.jl`
- `24_mapped_cartesian_hydrogen.jl`
- `25_mapped_cartesian_hydrogen_backends.jl`
- `33_ordinary_cartesian_1s2_vee.jl`
- `29_hybrid_mapped_cartesian_hydrogen.jl`
- `34_hybrid_cartesian_1s2_vee.jl`
- `30_ordinary_sho_spectra.jl`

## Full curated guide

For the full grouped running guide, use the
[Example guide](../howto/example_guide.md).

For the advanced contraction, hierarchy, and prototype examples, use that full
guide rather than this shorter landing page.

For the workflow pages that explain these examples, also see:

- [Manual](../manual/index.md)
- [Current atomic branch](../explanations/current_atomic_branch.md)
- [Current ordinary branch](../explanations/current_ordinary_branch.md)
