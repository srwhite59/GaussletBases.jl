# Bases and mappings

If you want the workflow first, start with:

- [First radial workflow](../tutorials/first_radial_workflow.md)
- [Recommended atomic setup](../howto/recommended_atomic_setup.md)
- [Example guide](../howto/example_guide.md)

## Small doctest

```jldoctest bases_and_mappings
julia> using GaussletBases

julia> basis = build_basis(MappedUniformBasisSpec(:G10; count = 5, mapping = AsinhMapping(c = 0.1, s = 0.2)));

julia> length(basis)
5
```

## Basis recipes and construction

```@docs
UniformBasisSpec
MappedUniformBasisSpec
HalfLineBasisSpec
RadialBasisSpec
recommended_xgaussians
build_basis
```

## Coordinate mappings

```@docs
IdentityMapping
AsinhMapping
fit_asinh_mapping_for_extent
fit_asinh_mapping_for_strength
```

The old one-dimensional COMX-cleaned hybrid constructor is intentionally
omitted from the public reference. It remains only as legacy/internal
experimental code for surrogate comparisons and historical regression checks.
