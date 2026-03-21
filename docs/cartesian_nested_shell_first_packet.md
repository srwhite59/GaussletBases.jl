# Cartesian Nested Shell First Packet

This note records the next narrow implementation step after the first local
`doside` and face-product primitive.

The first primitive worked:

- the one-dimensional `doside` contraction sits in the right legacy family
- the first `x-y` face product is structurally clean
- opposite face interiors remain disjoint on the finalized PGDG fixed block

So the next goal is still narrow, but one layer higher:

- assemble one first shell-level nested fixed space from face primitives
- propagate the carried operator packet through the same local contractions

This is still **not** recursive nesting yet.

## First shell scope

The first shell object is deliberately simple:

- one opposite-face pair of `x-y` faces
- one shell-level contraction matrix built by concatenating those face spaces
- one transformed shell-level packet on that same space

That is enough to establish the shell-level architecture without opening a
full rectangular-shell tree.

## Packet propagation scope

The first propagated shell packet carries:

- overlap
- kinetic
- `x`, `y`, and `z` position operators
- `x^2`, `y^2`, and `z^2` operators
- Gaussian-factor term matrices
- pair-factor term matrices

The purpose is to show that the current PGDG/QW-PGDG fixed-line packet can be
pushed through the nesting contractions cleanly, not to settle the final
consumption path yet.

## Structural role

The shell object is meant to be the first bridge between:

- local interval/face contractions
- and a later nested Cartesian fixed block that could be consumed by the
  existing Cartesian assembly logic

What still remains after this step is:

- broader shell coverage
- recursive shell or box nesting
- and the final adaptation of the existing Cartesian/QW-PGDG assembly to a
  nested fixed-space packet
