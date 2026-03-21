# Ordinary Cartesian GG Pair-Factor Legacy Optimization Note

This note records the next optimization pass after the first term-first `gg`
assembly improvement.

That earlier pass was worthwhile:

- the three-dimensional Coulomb-expanded `gg` assembly now keeps the Coulomb
  term index as the short inner reduction
- the resulting loops look much more like the legacy
  `sum(coulco .* Fx[:,ix,jx] .* Fy[:,iy,jy] .* Fz[:,iz,jz])` pattern

But the timing diagnosis then shifted one level lower.

The main bottleneck was no longer the three-dimensional assembly. It was the
one-dimensional mapped-gausslet `gg` pair-factor builder itself.

## What legacy pattern this pass follows

This pass moves the shared mapped-gausslet `gg` pair-factor construction closer
to the legacy `Velel1d` route.

The relevant old pattern in `PureGaussianGausslet.jl` is:

1. build a one-dimensional raw Gaussian layer
2. form term-resolved one-dimensional Coulomb data there
3. contract that data through the one-dimensional gausslet transforms
4. divide by the resulting one-dimensional weights
5. use the Coulomb-expansion sum only as the short inner reduction in the
   three-dimensional assembly

For the current mapped ordinary basis, that means:

- do not build the active `gg` pair factors from a midpoint-grid kernel in the
  final basis
- instead build a legacy-style mapped pure-Gaussian proxy raw layer at the
  mapped centers and local widths
- fold the local Jacobian amplitude into the contraction matrix
- evaluate the one-dimensional pair factors analytically on that proxy raw
  layer
- then contract them into the working one-dimensional gausslet basis

This is still only the mapped-gausslet `gg` path.
The newer term-first three-dimensional assembly remains in place.

## What implementation was replaced

In the shared numerical `gg` bundle, the active pair-factor path was changed
from:

- `_basis_space_pair_factor_terms(...)`

to:

- `_mapped_legacy_proxy_layer(...)`
- `_pair_gaussian_factor_matrices(proxy_layer; ...)`

So the old midpoint-grid basis-space pair-factor builder is now demoted to a
diagnostic/fallback helper. It is no longer the active path in
`_mapped_ordinary_gausslet_1d_bundle(...)` for
`backend = :numerical_reference`.

## Small iteration timing

On a deliberately small iteration case

- `count = 3`
- `s = 0.8`
- `xmax = 6`
- first Coulomb term only

the before/after comparison was:

- old `_basis_space_pair_factor_terms(...)`: `3.3123 s`
- new `_mapped_legacy_proxy_layer(...)`: `0.2202 s`
- new analytic proxy pair-factor build: `0.5792 s`

So the active one-dimensional pair-factor construction became much faster on
the small diagnostic case.

## Updated light reference-case timing

On the light paper-like reference case

- `count = 9`
- `s = 0.8`
- `xmax = 6`
- He `cc-pVTZ` `s`
- `interaction_treatment = :mwg`

the shared gausslet-side bundle subphases are now approximately:

- `basis_representation(...)`: `15.04 s`
- direct scalar `x^2` plus Gaussian-factor build: `38.69 s`
- `_mapped_legacy_proxy_layer(...)`: `0.21 s`
- analytic proxy pair-factor build: `0.42 s`

with the full `ordinary_cartesian_qiu_white_operators(...; timing = true)`
constructor now printing:

- `shared gausslet-side 1D bundle: 51.98 s`

So phase 1 now completes in a real, measurable regime.
It is still not fast, but it is no longer dominated by the old numerical
pair-factor builder.

## What the next bottleneck becomes

After this correction, the next bottleneck is no longer the shared mapped
gausslet `gg` pair-factor path.

The next phase that failed to clear in the same light reference case was:

- `split gausslet-Gaussian raw-block assembly`

So the next expensive stage now appears to be the Qiu-White-specific
`ga` split-block construction, not the shared `gg` pair-factor builder.

That is the right direction:

- the generic mapped-gausslet `gg` bottleneck was real and is now much smaller
- the remaining slowdown is now closer to genuinely Qiu-White-specific work
- only now does it make sense to inspect that next stage instead of continuing
  to optimize the old `gg` numerical kernel path
