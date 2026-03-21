# Cartesian Nested Face First Primitive

This note records the first concrete implementation step under
`docs/src/algorithms/cartesian_nested_face_construction.md`.

The intended start is narrow:

- one finalized QW-PGDG fixed line
- one local one-dimensional `doside` primitive on one interval
- one simple face-product construction from two tangential `doside` spaces

This is **not** full recursive nesting yet.

## What Was Added

The repo now has an internal one-dimensional side-contraction helper on the
finalized PGDG fixed line:

- `src/cartesian_nested_faces.jl`
- `_nested_doside_1d(...)`
- `_CartesianNestedDoSide1D`

The helper follows the legacy `getside(...)` pattern most closely:

- start from local weights on one interval
- build the retained span by repeated multiplication by the local coordinate
  operator
- metric-orthogonalize each added direction
- diagonalize the projected local coordinate operator
- sign-fix the resulting side functions with the local weights

So the retained side space is meant to keep the low-order local coordinate
content of the interval, not to perform an arbitrary SVD truncation.

## First Face Product

The first simple three-dimensional nested object is:

- `_nested_xy_face_product(...)`
- `_CartesianNestedXYFace3D`

This follows the `x-y face` branch of the legacy code:

- build one local `doside` contraction in `x`
- build one local `doside` contraction in `y`
- take their product space on one fixed `z` index

The current first-pass face support is restricted to the supplied tangential
intervals and the fixed `z` index, so opposite faces remain disjoint by
construction.

## Role Of The Fixed Line

This first primitive sits on top of the stabilized pre-nesting fixed line:

- finalized PGDG/QW-PGDG one-dimensional block
- already COMX-cleaned in the base PGDG basis-realization layer
- already suitable for the current pre-nesting Cartesian convergence family

So this pass does not reopen mapping-family or residual-Gaussian questions. It
only establishes the next contraction/localization layer on top of that fixed
line.

## What This Does Not Do Yet

This first implementation deliberately stops before:

- recursive shell nesting
- deeper box hierarchies
- sliced transverse-only contractions
- operator-packet propagation through a full nested tree

Those later steps should build on this local `doside` and face-product
primitive rather than replacing it.
