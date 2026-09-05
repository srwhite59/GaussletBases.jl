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

## Operator conventions and resources

The diagonal-approximation convention and finite Coulomb expansion are
operator-construction inputs, not guarantees of continuum accuracy.

```@docs
AbstractDiagonalApproximation
IntegralDiagonal
position_matrix
coulomb_gaussian_expansion
```

## Bond-aligned QW and Coulomb resources

The QW container is a bond-aligned mapped-product basis built through the
public homonuclear or heteronuclear constructors, not a general geometry type.
Its Coulomb resource is a finite Gaussian expansion, not an exact operator or
a universal-interval accuracy claim.

```@docs
BondAlignedDiatomicQWBasis3D
CoulombGaussianExpansion
```

## Experimental QW geometry diagnostics

These route-specific expert functions return read-only inspection records;
they are neither a general molecular-geometry API nor a stable serialization
schema. For a nested-diatomic source, the diagnostic preserves and inspects
that exact source. The basis overload instead
constructs the normalized nested source through its shared-shell,
endcap/panel, packet, retention, split, and protection controls, so it is
not a lightweight geometry-only query. Its result summarizes geometry and
retention contracts, shell counts and dimensions, actual shared-shell
provenance, fixed dimension, contract audit, and atom-growth anatomy.

The chain and square diagnostics only inspect existing experimental bases.
They report their coordinates, mapped centers, local spacings, mapping kinds,
monotonicity, and symmetry checks; the square also reports x-y center
agreement. They construct no nested source, Hamiltonian, solver input,
arbitrary-orientation model, or heteronuclear model.

```@docs
bond_aligned_diatomic_nested_geometry_diagnostics
bond_aligned_homonuclear_chain_geometry_diagnostics
axis_aligned_homonuclear_square_lattice_geometry_diagnostics
```

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
