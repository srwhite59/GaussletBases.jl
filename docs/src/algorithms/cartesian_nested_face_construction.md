# Cartesian Nested Face Construction

## Pseudocode

1. Start from one finalized Cartesian fixed line.
   In the current repo this means the finalized QW-PGDG fixed line used as the
   parent line for later atomic nesting work.
   Code: `src/ordinary_mapped_backends.jl`, `src/ordinary_qiu_white_rg.jl`

2. Choose one one-dimensional interval on that fixed line.
   The first primitive is local, not global: one interval of finalized parent
   functions that will be compressed into a smaller side space.
   Code: `src/cartesian_nested_faces.jl`

3. Build one local `doside` contraction on that interval.
   The retained span is not an arbitrary truncation. It is the low-order local
   coordinate content of the interval.
   Code: `src/cartesian_nested_faces.jl`
   Legacy model:
   `GaussletModules/PureGaussianGausslet.jl` `getsideu(...)`, `getside(...)`

4. Form the retained side space by the legacy local moment/coordinate rule.
   In the current pattern:
   - start from local weights
   - build directions by repeated multiplication by the local coordinate
   - orthogonalize in the local metric
   - diagonalize the projected local coordinate operator
   - sign-fix with the local weights
   Code: `src/cartesian_nested_faces.jl`

5. Build one tangential face product from those side spaces.
   For one rectangular face:
   - contract the first tangential interval with `doside`
   - contract the second tangential interval with `doside`
   - take the product of the two retained side spaces
   Code: `src/cartesian_nested_faces.jl`

6. Restrict to face interiors so different faces remain disjoint.
   The first primitive shell language already depends on disjoint support
   between opposite or neighboring faces.
   Code: `src/cartesian_nested_faces.jl`

7. Assemble the first shell packet from those face pieces.
   The shell packet carries:
   - one shell coefficient map
   - shell overlap
   - one-body packet pieces
   - Gaussian-factor and pair-factor packet pieces
   Code: `src/cartesian_nested_faces.jl`

8. Stop at the first shell packet.
   This page is only the primitive page:
   - `doside`
   - tangential face products
   - first shell packet
   The full landed atomic nonrecursive route now lives in:
   [Cartesian nested atomic nonrecursive route](cartesian_nested_atomic_nonrecursive_route.md)

## References

- Legacy primitive model: `GaussletModules/PureGaussianGausslet.jl`
  - `getsideu(...)`
  - `getside(...)`
  - `getsidexyznew(...)`
  - face branches around the `x-y face`, `x-z face`, and `y-z face` code
- Upstream fixed-line framing:
  [1D distorted-gausslet PGDG refinement hierarchy](distorted_gausslet_pgdg_refinement_hierarchy.md)
- Landed nonrecursive atomic route:
  [Cartesian nested atomic nonrecursive route](cartesian_nested_atomic_nonrecursive_route.md)

## What This Frames

This page records the local primitive language for Cartesian nesting:

- one-dimensional `doside`
- two-dimensional tangential face products
- one first shell packet with transferred operator data

It is no longer the source-of-truth page for the whole landed atomic nonrecursive
route.

## Current Repo Status

The repo now has the primitive pieces described on this page:

- dedicated one-dimensional `doside` helpers
- tangential face-product constructors
- first shell packet propagation on those face pieces

The repo has moved beyond this primitive stage for the active atomic line. For
the landed nonrecursive fixed-block route, coverage rule, shell-plus-core,
corrected complete-shell source, and fixed-block adapter, use:
[Cartesian nested atomic nonrecursive route](cartesian_nested_atomic_nonrecursive_route.md)

## Relation To Other Pages

- This page is the primitive local-contraction page.
- [Cartesian nested atomic nonrecursive route](cartesian_nested_atomic_nonrecursive_route.md)
  records the landed nonrecursive atomic fixed-block route built from these
  primitives.
- [Qiu-White residual-Gaussian route](qiu_white_residual_gaussian_route.md)
  records the later hybrid completion once a fixed block is already in hand.

## Implementation Notes

Recommended code-comment style:

```julia
# Alg Nested-Face step 3: Build a local 1D doside contraction on one interval.
# See docs/src/algorithms/cartesian_nested_face_construction.md.
```

Guidelines:

- keep this page focused on local primitives
- do not use it as the source-of-truth page for the full atomic nonrecursive
  nesting route
- keep later fixed-block assembly and residual-Gaussian completion on their own
  algorithm pages
