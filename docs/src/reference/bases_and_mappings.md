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
white_lindsey_atomic_mapping
fit_asinh_mapping_for_extent
fit_asinh_mapping_for_strength
```

For the one-center atomic White-Lindsey-style path, the repo now exposes
`white_lindsey_atomic_mapping(Z=..., d=..., tail_spacing=...)` directly. This
route is `d`-driven: choose the physical core spacing `d` first, then resolve
any later `count` choice against that already-chosen map. Do not treat
`count -> s` as the front-door contract for this one-center atomic family.

The separate multi-center legacy `getmapping(...)` line remains a different
combined inverse-sqrt-density construction. In the repo that path corresponds
to `CombinedInvsqrtMapping`, not to the one-center atomic helper above.

The old one-dimensional COMX-cleaned hybrid constructor is intentionally
omitted from the public reference. It remains only as legacy/internal
experimental code for surrogate comparisons and historical regression checks.
