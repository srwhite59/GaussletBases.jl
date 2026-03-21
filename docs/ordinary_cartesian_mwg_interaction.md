# Ordinary Cartesian IDA: Matched-Width-Gaussian (MWG) Residual Interaction

This pass asks a narrower question than the earlier residual-Gaussian note.

The latest legacy-He `s` supplement checks showed that the remaining hybrid
ordinary issue is no longer "we only used a toy Gaussian add-on." A real
Gaussian supplement can make the one-body/core description essentially
hydrogenic while still leaving the raw `1s^2` density-density /
integral-diagonal approximation (IDA) scalar mixed and non-monotone.

That is exactly why the next comparison should focus on interaction modeling,
not on whether the Gaussian supplement itself is still toy-sized.

## Why MWG is the right next step

The earlier residual-Gaussian nearest-center path is useful as a baseline, but
it is not the paper's preferred route. It is a GGT-like approximation:

- build residual Gaussians by orthogonalizing the added Gaussian channel
  against the ordinary backbone
- assign each residual direction to the nearest ordinary center
- reuse ordinary-center interactions there

That can be acceptable when the residual channel occupancy is tiny, but it is
cruder than the paper's recommended matched-width-Gaussian (MWG)
approximation.

The current repo is now ready for the MWG pass because it already has:

- a stable hybrid ordinary basis route
- a legacy-informed He `s` Gaussian supplement adapter
- exact one-dimensional overlap, `x`, and `x^2` matrices in the seed space
- the existing residual-Gaussian seed construction and transfer path

## What MWG means here

For each residual Gaussian direction:

1. orthogonalize the added Gaussian channel against the ordinary backbone
2. compute the exact first and second one-dimensional moments of that
   residual direction
3. derive an effective Gaussian center and width from those moments
4. use those effective Gaussian orbitals to build the residual interaction
   seed
5. transfer that seed interaction back into the final localized hybrid basis

So this pass changes the interaction model, not the one-body hybrid basis
construction.

## Why the scalar target is still useful

The same simple scalar remains the right first validation target:

- take the noninteracting one-body ground orbital
- doubly occupy it
- evaluate the resulting `1s^2` density-density / IDA interaction
  expectation
- compare against the hydrogenic reference `⟨Vee⟩ = (5/8) Z`

For `Z = 2`, that scalar target is `1.25`.

The point is not that a more complete Gaussian supplement should monotonically
improve this scalar by itself. A stronger Gaussian supplement mainly improves
the one-body/core description. The scalar can move in mixed ways unless the
interaction treatment is also reasonable.

## Scope of this pass

This pass stays narrow:

- current hybrid basis API unchanged
- current legacy He `s` adapter reused as-is
- no He solver yet
- no ESOI/EGOI/cusp-correction layer yet
- no broad Gaussian-basis framework

The comparison is:

- `:combined_basis`
- `:residual_gaussian_nearest`
- `:residual_gaussian_mwg`

in the same paper-like backbone regime with real He `s` supplements.

## Current interpretation

If MWG behaves better than nearest-center on the same legacy-He supplement, it
becomes the right residual-Gaussian interaction path for the hybrid ordinary
branch.

That still does **not** make the branch solver-ready by itself. It only
clarifies the next scientifically faithful interaction approximation to use
before any first real ordinary He solve.
