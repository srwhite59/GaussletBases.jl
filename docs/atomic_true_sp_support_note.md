# Atomic True-SP Support Note

This pass closes the atomic `lmax = 1` support gap for the active QW routes.

## What changed

The new active supplement object on the QW side is the explicit internal 3D
shell representation:

- `_AtomicCartesianShellSupplement3D`

It is built from the shared loader object:

- `LegacyAtomicGaussianSupplement`
- `legacy_atomic_gaussian_supplement(atom, basis_name; lmax = 1)`

and expands the named-basis shell metadata into explicit atomic-centered
Cartesian orbitals:

- `s`
- `p_x`
- `p_y`
- `p_z`

## Active paths

`lmax = 1` is now a true active physical path for:

- ordinary QW
- nested QW

Both consume the same shared atomic named-basis supplement object and switch to
the explicit 3D shell route when non-`s` shells are present.

The older centered separable 1D supplement route remains active only for
`lmax = 0`.

## What stayed narrow

This pass is still atomic-only:

- one atom at the origin
- named basis plus `lmax`
- no diatomic placement logic
- no `d` or higher shells

The one-dimensional hybrid builder also remains honestly `s`-only. It still
rejects `l > 0`, because it does not carry a true explicit 3D shell route.

## Legacy model matched

This follows the legacy atomic line in the right structural sense:

- load named basis shells up to `lmax`
- expand them into explicit Cartesian shell content
- carry true 3D shell functions in the active residual-Gaussian route

That is the atomic finishing step needed before diatomics.
