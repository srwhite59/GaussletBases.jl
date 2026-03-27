# Radial Paper-Parity Prototype Note

This note records the named manuscript radial boundary prototype now owned by
`GaussletBases`.

## Named prototype

The settled manuscript prototype is:

- `:paper_parity_g10_k6_x2`

This is a first-class cached prototype/profile above the generic radial family
front door. It is not just an alias for `:G10`, and it is not an atom-specific
mapped basis.

The fixed prototype contract is:

- family: `:G10`
- reference spacing: `1.0`
- odd seed half-width: `L = 24`
- even-tail parameter: `K = 6`
- paper-parity `x`-Gaussians:
  - `0.09358986806`
  - `0.02357750369`

The manuscript build/evaluation provenance recorded with the cached prototype
is:

- `h = 0.001`
- `sigma = 3`
- `s0 = 6.5`
- `rmax_int = 80`

## Cached endpoints

The vendored cache artifact lives at:

- `data/radial/paper_parity_g10_k6_x2.jld2`

It stores both prototype endpoints:

- canonical scientific endpoint:
  - final coefficients in the basis of the seed gausslets plus the two
    `x`-Gaussians
- runtime endpoint:
  - final coefficients expanded into the underlying half-line Gaussian
    primitives plus the two `x`-Gaussians

The public loader is:

- `radial_boundary_prototype(:paper_parity_g10_k6_x2)`

The loaded prototype can then be realized as a runtime `RadialBasis` through:

- `build_basis(prototype; mapping = ...)`

## Strict semantics

This named prototype uses strict manuscript semantics:

- analytic `S` and `X` construction on the undistorted half-line prototype
- high-precision prototype build path
- no silent mode loss

If the expected final mode count is not retained exactly, the build must fail
rather than silently dropping a direction.

That strict policy is specific to the named manuscript prototype. It is not a
global promise about every generic runtime radial build in the package.
