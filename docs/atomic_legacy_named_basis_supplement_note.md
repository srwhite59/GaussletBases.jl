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

## Important current limitation

The present ordinary-QW and nested-QW residual consumers still use a separable
centered 1D supplement route underneath. So in this pass the new shared object
stores:

- the full named-basis shell metadata up to `lmax`
- plus the active centered `s` projection used by the current analytic
  primitive path

That means:

- `lmax = 0` is the fully active and validated case
- `lmax = 1` already loads and preserves the right atomic shell metadata
- but it is not yet a true active `p`-enabled QW/nested supplement path

So this pass fixes the atomic interface footing without pretending the current
separable residual machinery already supports the full angular supplement.

## Scope kept intentionally atomic-only

This pass does not add:

- diatomic or molecular Gaussian placement
- a general molecule-wide basis subsystem
- a new nesting geometry policy

It only puts ordinary QW and nested QW onto the same named-basis-plus-`lmax`
atomic supplement interface so the next atomic anchor comparison can be made on
the right footing.
