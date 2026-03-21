# Ordinary Cartesian Qiu-White GA Cross Legacy Optimization Note

This note records the next optimization pass after the shared `gg` scalar and
pair-factor corrections.

Those shared-bundle fixes were real:

- the mapped-gausslet `gg` pair-factor path moved onto a legacy-style mapped
  Gaussian proxy raw layer
- the shared `gg` scalar/factor path moved there too
- the Qiu-White constructor no longer spent most of its time on the old mixed
  raw-layer or shared `gg` mistakes

That exposed the next bottleneck more clearly: the Qiu-White-specific
gausslet-Gaussian `ga` cross-block route.

The old active `_qwrg_cross_1d_blocks(...)` path was still too numerical:

- midpoint-grid sampling on the final mapped basis
- basis and Gaussian value matrices on that grid
- dense `distance^2` kernels
- repeated `basis' * (kernel * gaussian)` work for each expansion term

That is not the legacy pattern.

## Legacy pattern copied here

The old `PureGaussianGausslet.jl` `doGTO1d(...)` logic forms the off-diagonal
gausslet-Gaussian blocks by:

- keeping the gausslet side in a raw Gaussian proxy basis
- evaluating the needed Gaussian-Gaussian primitive matrices analytically
- contracting those primitive cross blocks into the working gausslet basis

This pass follows that same structure for the Qiu-White `ga` blocks:

- keep `gg` from the shared gausslet bundle
- keep `aa` analytic
- build the `ga` blocks from a dedicated mapped Gaussian proxy layer on the
  gausslet side
- contract only on the gausslet side

The newer term-first three-dimensional Coulomb-expanded assembly is kept
unchanged.

## What implementation was replaced

The active Qiu-White `ga` route no longer uses the midpoint-grid builder as
its constructor path.

It is replaced by:

- `_qwrg_cross_1d_blocks(...)`
  - active path
  - legacy-style mapped proxy contraction
- `_qwrg_cross_1d_blocks_midpoint(...)`
  - retained only as a diagnostic baseline

The new active route uses:

- the shared `_mapped_legacy_proxy_layer(...)`
- analytic Gaussian-Gaussian primitive formulas for
  - overlap
  - kinetic
  - position
  - `x^2`
  - Gaussian-factor blocks
  - pair-factor blocks
- left-side contraction back into the working gausslet basis

## Small warmed timing

On a small development case

- `count = 3`
- `s = 0.8`
- `xmax = 6`
- He `cc-pVTZ` `s`
- first Coulomb term only

the warmed comparison was:

- old midpoint path `_qwrg_cross_1d_blocks_midpoint(...)`: did not finish
  within a `20 s` observation window
- new proxy path `_qwrg_cross_1d_blocks(...)`: about `1.88885e-4 s`

So the active `ga` cross-block path is now effectively free at that scale,
while the demoted midpoint route remains unusably expensive even on a tiny
case.

## Updated light reference-case timing

On the light paper-like reference case

- `count = 9`
- `s = 0.8`
- `xmax = 6`
- He `cc-pVTZ` `s`
- `interaction_treatment = :mwg`

the warmed constructor timing is now:

- shared gausslet-side 1D bundle: `14.668781604 s`
- split gausslet-Gaussian raw-block assembly: `0.421365972 s`
- 3D gausslet overlap assembly: `0.037123932 s`
- residual-space construction: `0.123772256 s`
- 3D gausslet one-body assembly: `0.671861544 s`
- raw one-body transform: `0.034547631 s`
- raw moment-matrix assembly: `0.164400885 s`
- 3D gausslet interaction assembly: `0.020893667 s`
- RG interaction assembly: `0.402140442 s`
- total: `16.544887933 s`

So the `ga` cross-block route is no longer the main problem.

## What the next bottleneck becomes

After this correction, the next bottleneck is no longer Qiu-White-specific.

The dominant remaining cost is again:

- the shared mapped-gausslet one-dimensional bundle

and within the remaining Qiu-White-specific stages the largest single phase on
this light case is now:

- `3D gausslet one-body assembly` at about `0.67 s`

That means the active `ga` midpoint-grid route has been removed as a serious
performance problem, and any next optimization should be chosen from the now
much smaller set of genuinely exposed costs.
