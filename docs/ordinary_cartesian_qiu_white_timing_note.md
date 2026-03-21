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
