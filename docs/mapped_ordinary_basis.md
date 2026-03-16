# Mapped full-line ordinary bases

This note describes the first public package route for the mapped ordinary
Cartesian branch.

The key idea is still the same one established in the Coulomb-expansion note:

- use a short Gaussian expansion of `1/r`
- exploit separability
- build one-dimensional Gaussianized factors
- assemble the 3D hydrogen operator from those separable pieces

The difference is that the one-dimensional ordinary basis is now a public
globally mapped basis object rather than a temporary analysis construction.

## 1. The public 1D mapped basis

The package now provides:

- `MappedUniformBasisSpec`
- `MappedUniformBasis`

This is the full-line analogue of the ordinary uniform basis, but with one
global coordinate map applied to the shared primitive layer.

For example:

```julia
map = fit_asinh_mapping_for_extent(npoints = 9, xmax = 6.0)

basis = build_basis(MappedUniformBasisSpec(:G10;
    count = 9,
    mapping = map,
    reference_spacing = 1.0,
))
```

This is still a **1D** basis object. It is not a general 3D framework by
itself.

## 2. Why this is the right packaging step

The successful Asinh hydrogen study showed that the ordinary Cartesian path
already has a clean scientific story:

- one global 1D map on each axis
- the usual gausslet stencil/contraction machinery
- Coulomb expansion first, not 3D grid first

That is enough to justify making the mapped full-line basis a real package
object before building anything broader.

## 3. The Asinh helper

For the first hydrogen-style studies, the package also provides the narrow
helper

- `fit_asinh_mapping_for_extent(...)`

which builds the symmetric one-parameter family

```text
AsinhMapping(c=s, s=s, tail_spacing=10.0)
```

so that an odd-count mapped full-line basis has its outer centers at
`x = ±xmax`.

This helper is intentionally narrow. It is there to support the first mapped
ordinary Cartesian hydrogen path, not to define every possible mapped basis
workflow.

## 4. What this does not mean

This does **not** mean the package now has:

- a general 3D tensor-product basis framework
- a 3D grid-first ordinary hydrogen path
- a full Cartesian He workflow
- nested contraction machinery on the ordinary branch

The public step here is much smaller:

- a usable mapped 1D ordinary basis
- a clean way to reproduce the successful ordinary Cartesian hydrogen route

## 5. How to read this in the broader package story

The radial path is still the main public entry point.

The mapped ordinary full-line basis should be read as:

- the first public packaging of the ordinary Cartesian hydrogen route
- still on the Coulomb-expansion-first architecture
- useful for hydrogen-focused validation before any He or broader Cartesian
  workflow
