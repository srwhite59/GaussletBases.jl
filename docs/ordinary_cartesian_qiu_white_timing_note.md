# Qiu-White Reference Timing Note

This note records a narrow timing/debugging pass on the separate Qiu-White
residual-Gaussian reference constructor.

The earlier cross-block correction was real:

- the gausslet-Gaussian cross blocks should come from the contracted raw-space
  1D matrix route
- not from midpoint-grid cross-quadrature in the final basis

After that correction, the remaining slowdown is no longer best explained by
the old cross-block mistake. But the next timing pass changed the diagnosis
again: the dominant cost is still earlier than the dense 3D Qiu-White stages.

The goal of this pass is therefore narrow:

- fix any stale correctness issue before trusting timings
- add one opt-in timing path on the Qiu-White reference constructor only
- identify the dominant remaining phase on a single light paper-like He case

This is not a broad timing framework for the library. It is only a debugging
aid for the Qiu-White reference path while that path is still being brought
into algorithmic alignment with the paper and the legacy implementation style.

## What the first timing pass actually found

The first timing/debugging pass changed the diagnosis.

On the light paper-like He case

- `count = 9`
- `s = 0.8`
- `cc-pVTZ`
- `interaction_treatment = :mwg`

the constructor did not even finish the first timed phase,
`exact 1D raw-block build`, within the initial observation window.

That points upstream of the dense 3D assembly.

A backend probe showed why: the mixed raw one-dimensional layer used by
`_qwrg_block_matrices(...)` is still selecting
`_NumericalPrimitiveMatrixBackend`, because the gausslet part of that raw layer
is built from `Distorted{Gaussian, AsinhMapping}` primitives.

So the present performance bottleneck is not mainly the later 3D `kron`-style
assembly that motivated this pass. It is the fact that the supposed
"exact contracted 1D raw-block build" is still numerical at the primitive
operator level for the mixed raw layer.

That led to the split-block correction:

- remove the mixed raw one-dimensional layer
- keep the `gg / ga / aa` split explicitly
- move the Qiu-White cross blocks onto the contracted 1D raw-space route

## What the timing pass found after the split-block refactor

After the split-block correction, the Qiu-White constructor was refactored
again so it no longer privately rebuilds the gausslet-side one-dimensional
data. It now consumes a shared mapped ordinary 1D gausslet bundle.

That changed the ownership of phase 1, but not yet its cost.

On the same light paper-like He case, the refactored constructor still did not
finish phase 1 within a two-minute observation window. Direct timings of the
shared gausslet-side bundle internals showed:

- `basis_representation(...)`: about `15.12 s`
- `_x2_matrix(basis)`: about `18.97 s`
- `gaussian_factor_matrices(...)`: about `5.08 s`
- `_pair_gaussian_factor_matrices(...)`: still running after `60 s`

So the dominant remaining cost is now much clearer:

- it is not the later dense 3D Qiu-White reference stages yet
- it is not mainly the dedicated `ga` cross-block route
- it is the shared mapped-gausslet one-dimensional bundle, especially the
  gausslet-side pair-factor builder

That means the next structural optimization should focus there first. Only
after that shared gausslet-side path is improved will it be meaningful to
judge how much of the remaining cost is intrinsic dense 3D reference work.

## What changed after the first shared-gg optimization pass

The next optimization pass moved the shared `gg` scalar/factor path closer to
the legacy term-first Coulomb-expansion style:

- direct basis-space `x^2` and Gaussian-factor construction
- term-first `F[g,i,j]` storage
- legacy-style short inner reduction in the `gg` 3D Coulomb-expanded assembly

That produced a real improvement on the same light case:

- `basis_representation(...)`: about `15.43 s`
- direct basis-space `x^2` plus Gaussian-factor build: about `4.90 s`

compared with the earlier separate costs of about:

- `basis_representation(...)`: `15.12 s`
- `_x2_matrix(basis)`: `18.97 s`
- `gaussian_factor_matrices(...)`: `5.08 s`

So the shared `gg` scalar/factor side is no longer the main concern.

But the dominant bottleneck still did not move far enough. The new direct
basis-space pair-factor builder still did not finish within about `70 s` on
that light case, and the full Qiu-White constructor still did not get past
phase 1 in a reasonable time.

So the next bottleneck remains:

- the shared mapped-gausslet `gg` pair-factor builder

not yet the later dense 3D reference stages.

## What changed after the legacy-style 1D pair-factor correction

The next optimization pass followed the legacy ordinary code more closely at
the one-dimensional `gg` pair-factor level.

Instead of building the active pair-factor data by midpoint-grid kernel
application in the final mapped basis, the shared numerical `gg` bundle now:

- builds a legacy-style mapped pure-Gaussian proxy raw layer at the mapped
  centers and local widths
- folds the local Jacobian amplitude into the contraction matrix
- evaluates the one-dimensional pair factors analytically on that proxy raw
  layer
- contracts those raw matrices into the working one-dimensional gausslet basis

That removed the old numerical basis-space pair-factor builder from the active
shared `gg` path.

The small before/after check on

- `count = 3`
- `s = 0.8`
- first Coulomb term only

was:

- old `_basis_space_pair_factor_terms(...)`: `3.3123 s`
- new `_mapped_legacy_proxy_layer(...)`: `0.2202 s`
- new analytic proxy pair-factor build: `0.5792 s`

On the light paper-like He case

- `count = 9`
- `s = 0.8`
- `xmax = 6`
- He `cc-pVTZ`
- `interaction_treatment = :mwg`

the shared bundle subphases are now approximately:

- `basis_representation(...)`: `15.04 s`
- direct scalar `x^2` plus Gaussian-factor build: `38.69 s`
- `_mapped_legacy_proxy_layer(...)`: `0.21 s`
- analytic proxy pair-factor build: `0.42 s`

and the constructor timing hook now reports:

- `shared gausslet-side 1D bundle: 51.98 s`

So the earlier diagnosis changed again:

- the shared mapped-gausslet `gg` pair-factor builder is no longer the main
  bottleneck
- the next phase that failed to clear in the same observation window is now
  the Qiu-White-specific split gausslet-Gaussian raw-block assembly

That means the next structural optimization should no longer target the old
shared `gg` numerical pair-factor path first. The remaining slowdown has moved
to the `ga` side of the Qiu-White reference constructor.

## What changed after the legacy-style 1D scalar/factor correction

The next optimization pass followed the legacy ordinary code more closely for
the shared mapped-gausslet `gg` scalar/factor path as well.

Instead of building active `x^2` and Gaussian-factor data by midpoint-grid
basis-space sampling, the shared numerical bundle now:

- reuses the mapped Gaussian proxy raw layer introduced for the pair-factor
  correction
- evaluates `x^2` analytically on that proxy raw layer
- evaluates Gaussian-factor blocks analytically on that proxy raw layer
- contracts those raw matrices into the working one-dimensional gausslet basis

The old `_basis_space_scalar_factor_data(...)` route is therefore no longer
the active shared `gg` scalar/factor path.

On a small development case

- `count = 3`
- `s = 0.8`
- first two Coulomb terms

the scalar-side timing changed from:

- old `_basis_space_scalar_factor_data(...)`: `6.8826 s`

to:

- `_mapped_legacy_proxy_layer(...)`: `0.5479 s`
- `_mapped_legacy_proxy_scalar_data(...)`: `2.0202 s`

On the light paper-like He case

- `count = 9`
- `s = 0.8`
- `xmax = 6`
- He `cc-pVTZ`
- `interaction_treatment = :mwg`

the shared bundle subphases were approximately:

- `basis_representation(...)`: `25.76 s`
- `_mapped_legacy_proxy_layer(...)`: `0.51 s`
- `_mapped_legacy_proxy_scalar_data(...)`: `1.43 s`
- analytic proxy pair-factor build: `0.63 s`

and the constructor timing hook now reports:

- `shared gausslet-side 1D bundle: 20.99 s`

So the shared `gg` bundle is now materially below the earlier `~52 s` regime.

After that timing line printed, the next phase still did not clear within the
next observation window. So the next exposed bottleneck is now the
Qiu-White-specific split gausslet-Gaussian raw-block assembly, not the shared
`gg` scalar/factor side.

## What changed after the legacy-style GA cross-block correction

The next optimization pass then moved the active Qiu-White `ga` cross-block
route onto the legacy-style contracted proxy path as well.

Instead of building active `ga` blocks by midpoint-grid sampling in the final
mapped basis, the constructor now:

- reuses the mapped Gaussian proxy layer already introduced on the gausslet
  side
- evaluates the primitive gausslet-proxy / added-Gaussian cross matrices
  analytically for
  - overlap
  - kinetic
  - position
  - `x^2`
  - Gaussian-factor blocks
  - pair-factor blocks
- contracts those primitive cross blocks only on the gausslet side

The old midpoint-grid implementation was kept only as a diagnostic helper.

On a warmed small development case

- `count = 3`
- `s = 0.8`
- `xmax = 6`
- He `cc-pVTZ`
- first Coulomb term only

the comparison was:

- midpoint diagnostic path: did not finish within a `20 s` observation window
- active proxy cross path: about `1.88885e-4 s`

So the active `ga` cross-block route is no longer a serious cost center.

On the same light paper-like reference case

- `count = 9`
- `s = 0.8`
- `xmax = 6`
- He `cc-pVTZ`
- `interaction_treatment = :mwg`

the warmed constructor timing is now:

- `shared gausslet-side 1D bundle`: `14.668781604 s`
- `split gausslet-Gaussian raw-block assembly`: `0.421365972 s`
- `3D gausslet overlap assembly`: `0.037123932 s`
- `residual-space construction`: `0.123772256 s`
- `3D gausslet one-body assembly`: `0.671861544 s`
- `raw one-body transform`: `0.034547631 s`
- `raw moment-matrix assembly`: `0.164400885 s`
- `3D gausslet interaction assembly`: `0.020893667 s`
- `RG interaction assembly`: `0.402140442 s`
- total: `16.544887933 s`

So the diagnosis changed again:

- the Qiu-White-specific `ga` split raw-block assembly is no longer the next
  bottleneck
- the dominant remaining cost is back on the shared mapped-gausslet 1D bundle
- within the remaining Qiu-White-specific stages, the largest single phase on
  this light case is now the `3D gausslet one-body assembly`
