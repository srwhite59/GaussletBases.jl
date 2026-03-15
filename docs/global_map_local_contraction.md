> **Note for new users:** this is a narrow research note about the current 1D contraction / hierarchy direction.  
> It is **not** the best place to start if you are new to GaussletBases.  
> Start instead with `README.md`, `docs/first_radial_workflow.md`, and `docs/recommended_atomic_setup.md`.

# Global map and local contraction note

This note corrects the likely next direction after the current 1D hierarchy and
PGDG-flavored prototype work.

The important correction is:

**do not treat the hierarchy as a license to assign different coordinate maps to
different leaves**

That is not the clean historical nested construction, and it raises interface
problems too early.

## Why per-leaf maps are risky

True per-leaf local maps create several hard questions at once:

- neighboring leaves would have mismatched primitive families
- orthogonality across leaf boundaries becomes nontrivial
- localization near box boundaries may degrade
- one may need explicit buffer or halo regions
- one may need a later global repair step

Those are real design problems, not just implementation details. They should
not be introduced casually.

## How the historical nested idea actually worked

The cleaner historical picture was:

1. define one global coordinate distortion on the axis
2. build one common primitive layer in that globally distorted coordinate
3. define boxes or shells inside that common basis
4. do local contraction or orthogonalization inside boxes or shells
5. organize the retained local functions hierarchically

So the geometry was local, but the primitive layer was still global.

That difference matters.

## Correct next abstraction

The right next abstraction is:

**global map + local shell or box contraction**

not

**different map in each leaf**

In particular:

- the globally defined object should be one common primitive layer over a
  globally mapped 1D region
- the locally defined object should be a contraction or orthogonalization
  defined inside one box or shell of that common layer

## What object is global

The next implementation target should have one globally defined object, for
example in the spirit of:

- `GlobalMappedPrimitiveLayer1D`

Its job would be to store:

- the global coordinate map
- the common `PrimitiveSet1D`
- any global primitive metadata needed for local restriction
- the global primitive overlap or simple one-body matrices, if useful

This object should be shared by every box and shell.

It should also be treated as a meaningful simple basis or intermediate object
in its own right.

That matters conceptually:

- the globally mapped uncontracted primitive or product basis is not only
  scaffolding for later contraction
- it is itself a useful simple object
- and in some cases it may already be a reasonable endpoint to use directly

So local shell or box contraction should be viewed as a refinement layer built
on top of that global basis, not as the only reason the global layer exists.

## What object is local

The next local object should be box- or shell-specific, for example in the
spirit of:

- `LocalBoxContraction1D`
- `LocalShellContraction1D`

Its job would be to store:

- which box or shell it belongs to
- which primitive indices from the global layer are used locally
- the local contraction or orthogonalization matrix
- the retained local functions produced from that contraction

The local object should *not* own a separate coordinate map.

## Where orthogonalization and contraction happen

Orthogonalization or contraction should happen locally, but on submatrices cut
from the globally defined primitive layer.

That means the sequence should be:

1. restrict the global primitive layer to a box or shell window
2. build the local overlap or simple one-body matrices on that restricted set
3. do the local contraction or orthogonalization there
4. retain the resulting local functions in the hierarchy bookkeeping

So the contraction is local, but the primitive family is common.

## Leaves first, parent shells later

For the next narrow implementation step, contraction should happen on leaves
only.

That is the cleanest first proof because:

- the current hierarchy layer already exposes leaves cleanly
- leaf-local bookkeeping is simpler
- one can delay overlap between parent and child retained spaces

Parent-shell contraction should remain part of the design direction, but it is
better treated as the next layer after leaf-local contraction works.

So the short answer is:

- next implementation step: leaves only
- later extension: parent shells too

## How this differs from `LeafLocalPGDG1D`

`LeafLocalPGDG1D` is still useful, but conceptually it does something
different.

It:

- generates a separate local Gaussian family in each leaf
- uses an identity contraction map
- treats the leaf-local generation itself as the basis-building step

The corrected direction should instead:

- start from one common globally mapped primitive layer
- use the hierarchy only to define local contraction regions
- produce retained functions by local contraction of that shared layer

So `LeafLocalPGDG1D` should be viewed as:

- a useful toy or prototype for hierarchy-driven generation infrastructure
- a good test of provenance, grouping, and representation plumbing
- but probably not the right conceptual path to the true historical nested
  construction
- and less faithful to the historical nested idea than a shared global mapped
  layer plus local box or shell contraction

## Recommended next implementation target

The next narrow coding step should be:

1. build one common globally mapped 1D primitive layer
2. restrict that layer to one leaf box at a time
3. do one local contraction inside each leaf
4. expose the retained local functions and their contraction data through the
   existing hierarchy and representation scaffolding

This would validate the corrected conceptual path without yet committing to:

- parent-shell contraction
- historical nested driver logic
- geometry-aware grouping
- named basis-set support

That is the right next step if the goal is to get back onto the historical
nested trajectory cleanly.
