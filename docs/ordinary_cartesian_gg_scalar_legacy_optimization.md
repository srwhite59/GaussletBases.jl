# Ordinary Cartesian GG Scalar/Factor Legacy Optimization Note

This note records the next optimization pass after the legacy-style `gg`
pair-factor correction.

That earlier pair-factor pass was real:

- the active mapped-gausslet `gg` pair-factor path no longer used the old
  midpoint-grid basis-space kernel route
- it moved onto a legacy-style mapped Gaussian proxy raw layer
- the one-dimensional pair terms were then built analytically and contracted
  into the working gausslet basis

But that did not finish the shared `gg` bundle.

On the same light paper-like reference case, the shared bundle was still about
`52 s`, with the remaining cost dominated by the scalar/factor side:

- `basis_representation(...)`: about `15 s`
- direct scalar `x^2` plus Gaussian-factor build: about `39 s`

So the next step in this pass is to follow the legacy code more closely there
too.

## Legacy pattern copied here

The old code did not build `x^2` and one-body Coulomb factor blocks by
midpoint-grid sampling in the final mapped basis.

Instead it built them by contracting raw Gaussian matrices:

- `x1d = basv' * [olapx(...)] * basv`
- `x2d = basv' * [olapx2(...)] * basv`
- `Vpot1d[g, ...]` from raw Gaussian `intVpot(...)` blocks

then transformed those contracted one-dimensional blocks into the later
working basis.

This pass mirrors that structure in the shared mapped-gausslet `gg` bundle:

- keep the mapped Gaussian proxy raw layer already introduced for the
  pair-factor optimization
- use that same proxy layer for `x^2`
- use that same proxy layer for Gaussian-factor blocks
- keep the results in the shared bundle so both ordinary Cartesian IDA and the
  Qiu-White reference path benefit

## What implementation was replaced

For `backend = :numerical_reference`, the shared bundle in
`src/ordinary_mapped_backends.jl` no longer uses:

- `_basis_space_scalar_factor_data(...)`

as the active path for:

- `x^2`
- Gaussian-factor terms

Instead it now uses:

- `_mapped_legacy_proxy_layer(...)`
- `_mapped_legacy_proxy_scalar_data(...)`

which evaluates those matrices analytically on the mapped Gaussian proxy raw
layer and contracts them into the one-dimensional gausslet basis.

The term-first three-dimensional Coulomb-expanded assembly from the previous
pass is unchanged.

## Small iteration timing

On a small development case

- `count = 3`
- `s = 0.8`
- `xmax = 6`
- first two Coulomb terms

the scalar-side before/after comparison was:

- old `_basis_space_scalar_factor_data(...)`: `6.8826 s`
- new `_mapped_legacy_proxy_layer(...)`: `0.5479 s`
- new `_mapped_legacy_proxy_scalar_data(...)`: `2.0202 s`

So the scalar/factor side is materially faster on the same kind of small case.

## Updated light reference-case timing

On the light paper-like reference case

- `count = 9`
- `s = 0.8`
- `xmax = 6`
- He `cc-pVTZ` `s`
- `interaction_treatment = :mwg`

the shared bundle subphases are now approximately:

- `basis_representation(...)`: `25.76 s`
- `_mapped_legacy_proxy_layer(...)`: `0.51 s`
- `_mapped_legacy_proxy_scalar_data(...)`: `1.43 s`
- analytic proxy pair-factor build: `0.63 s`

and the Qiu-White timing hook now reports:

- `shared gausslet-side 1D bundle: 20.99 s`

The subphase and constructor totals are not identical because they were timed
in separate runs, but the practical conclusion is clear:

- the shared `gg` bundle is now materially faster than the earlier `~52 s`
  regime
- the scalar/factor side no longer dominates the shared bundle the way it did
  before

## What the next bottleneck becomes

After this correction, the shared mapped-gausslet `gg` scalar/factor path is
no longer the best next optimization target.

On the same light reference case, after the constructor printed

- `shared gausslet-side 1D bundle: 20.99 s`

the next phase still did not clear within the next observation window.

So the next exposed bottleneck is now the Qiu-White-specific split
gausslet-Gaussian raw-block assembly, not the shared `gg` scalar/factor path.
