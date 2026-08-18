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
- `38_qiu_white_reference_vee.jl`

`38_qiu_white_reference_vee.jl` is a slow reference example. Its nearest/GGT
path is the default public run; the MWG branch remains opt-in and experimental.

For the current public Cartesian base Hamiltonian path, run:

- `39_pqs_h2plus.jl`

It compares bounded PQS and White-Lindsey H2+ constructions through public API
only. See [Projected q-shells (PQS)](../manual/projected_q_shells.md) for scope
and the [Export layer](../reference/export.md#cartesian-base-hamiltonian) for
the facade contract.

The v0.2 public validation examples are:

- `40_screened_hartree_fixed_density.jl`
- `41_pqs_h2plus_table1.jl`

The first is a bounded supplied-field screened-Hartree assembly example. The
second is the slower fixed full-parent/PQS/White-Lindsey H2+ comparison used
for release validation. See [Reference-density Hartree screening](../manual/reference_density_hartree_screening.md)
and [Projected q-shells (PQS)](../manual/projected_q_shells.md).

The older 1D COMX-cleaned hybrid examples remain in `examples/` only as
legacy/internal experimental regressions and are intentionally omitted from
this public sequence.

## Full curated guide

For the full grouped running guide, use the
[Example guide](../howto/example_guide.md).

For the advanced contraction, hierarchy, and prototype examples, use that full
guide rather than this shorter landing page.

For the workflow pages that explain these examples, also see:

- [Manual](../manual/index.md)
- [Current atomic branch](../explanations/current_atomic_branch.md)
- [Current ordinary branch](../explanations/current_ordinary_branch.md)

For the lightweight point-cloud and path viewers used to inspect emitted 2D/3D
geometry datasets, also see [Visualization utilities](../howto/visualization.md).
