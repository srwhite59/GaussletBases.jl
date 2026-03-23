# Atomic Legacy Named-Basis Supplement Note

This pass replaces the old He-`s`-only supplement entry point as the active
atomic supplement interface.

## What changed

The shared atomic supplement object is now:

- `LegacyAtomicGaussianSupplement`

with the new loader:

- `legacy_atomic_gaussian_supplement(atom, basis_name; lmax = ..., ...)`

The old helper

- `legacy_s_gaussian_data(atom, basis_name; ...)`

still exists, but it is now only a thin compatibility wrapper for
`legacy_atomic_gaussian_supplement(...; lmax = 0)`.

## Legacy model followed

This matches the old `basisaddname` / `lmaxadd` idea in the legacy line:

- load a named atomic Gaussian basis
- keep shells up to `lmax`
- contract shell by shell in the `makeallcontractions(...)` style

That is the right atomic modernization target before comparing atomic anchors.

## Current active split

The shared object now supports two different active consumer styles:

- ordinary-QW and nested-QW use the full named-basis shell metadata through an
  explicit atomic-centered 3D Cartesian shell route for `s, p_x, p_y, p_z`
- the one-dimensional hybrid builder still uses the centered analytic `s`
  projection only

That means:

- `lmax = 0` is the fully active and validated case
- `lmax = 1` now loads the right atomic shell metadata and is a true active
  physical `SP` supplement path for ordinary-QW and nested-QW
- the one-dimensional hybrid builder still rejects `l > 0` honestly, because
  it remains a centered separable 1D route

So the shared atomic interface is now on the intended named-basis-plus-`lmax`
footing, while remaining honest about the narrower one-dimensional hybrid path.

## Scope kept intentionally atomic-only

This pass does not add:

- diatomic or molecular Gaussian placement
- a general molecule-wide basis subsystem
- a new nesting geometry policy

It puts ordinary QW and nested QW onto the same named-basis-plus-`lmax`
atomic supplement interface and makes true atomic `SP` support real there,
without broadening into molecular placement.
