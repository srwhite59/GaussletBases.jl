# Cartesian Nested Face Construction

## Pseudocode

1. Start from one finalized Cartesian fixed block.
   In the intended first practical route, this is the finalized QW-PGDG fixed
   block: the Cartesian product gausslet side after the PGDG-mediated one-body
   and two-electron auxiliary construction, including the final COMX cleanup of
   the base PGDG basis-realization layer.
   Code: current fixed-block inputs live mainly in `src/ordinary_mapped_backends.jl` and `src/ordinary_qiu_white_rg.jl`

2. Partition space into rectangular shells around the atom or atoms.
   The first practical atomic case is a centered rectangular shell sequence.
   Later generalizations may use deeper recursion or more general partitions,
   but the basic primitive is already present at the level of one shell and its
   faces.
   Code: framing for now; related infrastructure exists in `src/hierarchical_partitions.jl`

3. On one one-dimensional interval, define the local side contraction
   (`doside`) from `n` distorted functions to `m < n` retained functions.
   The side basis is built by keeping the low-order local directions of the
   interval and then localizing them with a local COMX/COMU step.
   The intended retained space is not an arbitrary SVD truncation. It is the
   low-order local polynomial or coordinate content of the interval.
   Legacy model: `GaussletModules/PureGaussianGausslet.jl` `getsideu(...)` and `getside(...)`

4. Build the one-dimensional retained side space by a local moment/coordinate
   process.
   In the legacy pattern:
   - start from local weights or local density information
   - build the retained span by repeated multiplication by a local coordinate
     (`u` or `x`)
   - orthogonalize each added direction against the previous ones
   - diagonalize the projected local coordinate operator
   - sign-fix the final side functions with the local weights
   This produces a small, ordered, localized side space for one interval.
   Legacy model: `GaussletModules/PureGaussianGausslet.jl` `getsideu(...)`, `getside(...)`

5. For one face of a rectangular shell, apply the one-dimensional `doside`
   construction in the two tangential directions.
   For an `x-y` face, contract the local `x` interval and the local `y`
   interval; for an `x-z` face, contract `x` and `z`; for a `y-z` face,
   contract `y` and `z`.
   Legacy model: face branches in `GaussletModules/PureGaussianGausslet.jl`

6. Form the two-dimensional face basis as the product of the two tangential
   side spaces.
   This is the first true nested Cartesian object: a product-space basis
   attached to one face interior.
   Legacy model: face assembly in `GaussletModules/PureGaussianGausslet.jl`

7. Use only the interiors of faces so different faces do not overlap.
   Distinct face spaces should remain disjoint by construction.
   This is part of what makes the later nested assembly simple and local.
   Code: framing for now; intended first implementation target

8. Assemble the nested fixed space from the face spaces and any needed retained
   core/interior space.
   The first implementation does not need the full recursive tree. A single
   shell with a small set of face spaces is enough to establish the machinery.
   Code: framing for now

9. Transform the carried operator packet through the same local contractions.
   The fixed block is not only a basis. The carried one-dimensional or
   separable operator data must be transformed consistently through the nesting
   contractions so that the existing Cartesian assembly logic can be reused.
   Code: intended future consumer path on top of `src/ordinary_mapped_backends.jl`

10. Extend later to shell recursion, box hierarchies, and sliced variants.
   Once the one-dimensional `doside` primitive and the first face-product
   construction are in place, the same logic extends naturally to:
   - deeper rectangular shell recursion
   - atomic box hierarchies
   - sliced bases that contract only in transverse directions
   This page stops before those later generalizations.

## References

- Legacy model: `GaussletModules/PureGaussianGausslet.jl`
  - `getsideu(...)`
  - `getside(...)`
  - `getsidexyznew(...)`
  - face branches around the `x-y face`, `x-z face`, and `y-z face` code
- Repo framing page:
  `algorithms/distorted_gausslet_pgdg_refinement_hierarchy.md`

## What This Frames

This page records the intended starting point for the Cartesian nesting line.

The main distinction is that the first real nesting primitive is not a global
three-dimensional contraction. It is the one-dimensional interval contraction
`doside`, followed by a face-product construction in the two tangential
directions.

This keeps the nesting logic local and makes the relationship to the legacy
code explicit.

## Current Repo Status

The repo now has a stabilized base QW-PGDG path that is suitable as the
starting point for this nesting line:

- active PGDG-mediated auxiliary construction is analytic in current use
- the base PGDG basis-realization layer includes its final COMX cleanup
- the current Cartesian fixed block is converging cleanly enough to serve as a
  pre-nesting starting point

What the repo now has is the first narrow primitive described on this page:

- a dedicated one-dimensional `doside` helper on the finalized fixed line
- one first simple `x-y` face-product constructor on that fixed line
- one first opposite-face shell object with shell-level packet propagation
- one first generalized six-face shell-packet interface

What the repo still does not yet have is the broader nested rollout:

- no recursive shell or box nesting
- no operator-packet propagation through a full nested hierarchy
- no recursive shell tree beyond the first generalized single-shell interface

So this page now stabilizes both the algorithm and the intended boundary of the
first implementation step.

## Relation To Other Pages

This page is downstream of, but separate from:

- `algorithms/qiu_white_residual_gaussian_route.md`
- `algorithms/distorted_gausslet_pgdg_refinement_hierarchy.md`

Those pages explain how to build the current Qiu-White / QW-PGDG starting
point. This page explains the next contraction/localization layer that should
sit on top of that starting point.

## Implementation Notes

The first implementation target should stay narrow:

- one finalized Cartesian fixed block
- one local one-dimensional `doside` primitive
- one simple shell or face test

The first implementation should follow the legacy `getsideu/getside` logic
closely rather than inventing a different contraction rule.

Recommended code-comment style:

```julia
# Alg Nested-Face step 3: Build a local 1D doside contraction on one interval.
# See docs/src/algorithms/cartesian_nested_face_construction.md.
```

The same comment rules used elsewhere in the Algorithms section apply here:

- keep step numbers aligned with this page
- keep wording close to the pseudocode
- include the docs path exactly
