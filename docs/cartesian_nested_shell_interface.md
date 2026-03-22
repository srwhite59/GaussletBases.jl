# Cartesian Nested Shell Interface

This note records the next narrow generalization step after the first
shell-packet proof of concept.

The first shell packet established the key numerical point:

- the nesting contractions can push the carried operator packet through cleanly

So the next goal is not recursion yet. It is to move from one special-case
opposite-face shell pair to a shell-packet interface that can carry multiple
face pieces uniformly.

## Goal

The immediate target is a shell-level object that is easier for a later
Cartesian/QW-PGDG consumer to read directly:

- one shell basis built from a uniform collection of face pieces
- one shell-level packet carrying the same transformed data fields as the
  current fixed-line packet
- one explicit map from face pieces to shell columns

That is closer to the shape the existing assembly logic will eventually want.

## What this generalization should cover

The first practical generalized shell interface should now cover:

- all six faces of one rectangular shell
- with the same local ingredients as before:
  - 1D `doside`
  - tangential face-product spaces
  - disjoint face interiors

The main architectural point is that the shell object should now be a uniform
multi-face container rather than a one-off opposite-face special case.

## Packet scope

The shell packet should continue to carry the same transformed operator data:

- overlap
- kinetic
- `x`, `y`, `z`
- `x^2`, `y^2`, `z^2`
- Gaussian-factor terms
- pair-factor terms

That keeps it aligned with the current PGDG/QW-PGDG packet architecture.

## What still remains after this step

Even after the multi-face shell interface lands, the repo still will not yet
have:

- recursive shell nesting
- shell-to-shell hierarchy propagation
- a direct nested consumer inside the existing Cartesian/QW-PGDG constructors

Those are the next layers after the generalized shell packet exists.
