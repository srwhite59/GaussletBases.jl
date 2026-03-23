# Operators and diagnostics

If you want the narrative workflow before the API details, start with:

- [First radial workflow](../tutorials/first_radial_workflow.md)
- [Recommended atomic setup](../howto/recommended_atomic_setup.md)

## Small doctest

```jldoctest operators_and_diagnostics
julia> using GaussletBases

julia> rb = build_basis(RadialBasisSpec(:G10; rmax = 8.0, mapping = AsinhMapping(c = 0.1, s = 0.2), xgaussian_count = 0));

julia> grid = radial_quadrature(rb; accuracy = :medium);

julia> size(kinetic_matrix(rb, grid)) == (length(rb), length(rb))
true
```

This tiny doctest uses `xgaussian_count = 0` to keep the reference example
lightweight. The recommended public atomic setup still uses the default
`xgaussian` supplement unless you have a specific reason to turn it off.

## Diagnostics and quadrature

```@docs
basis_diagnostics
radial_quadrature
quadrature_points
quadrature_weights
```

## One-body radial operators

```@docs
overlap_matrix
kinetic_matrix
nuclear_matrix
centrifugal_matrix
multipole_matrix
```

`centrifugal_matrix(rb, grid; l)` is a single fixed-`l` radial block. It is
not yet the full explicit atomic `(l,m)` Hamiltonian. The usual next step is:

```julia
radial_ops = atomic_operators(rb, grid; Z = Z, lmax = 2)
atom = atomic_one_body_operators(radial_ops; lmax = 2)
```

which repeats the fixed-`l` radial blocks over the explicit `YlmChannel` list.

## Radial atomic bundle

```@docs
RadialAtomicOperators
atomic_operators
centrifugal
multipole
```
