# Cartesian Nested Representation Completeness

This note records the architectural conclusion from the first nonrecursive
shell-sequence experiments.

The fixed-`nside` pass was useful because it made the real failure explicit:

- the problem is not mainly the final tiny contracted `3 x 3 x 3` core
- even a grow-only multi-shell sequence with an exact raw `5 x 5 x 5` core is
  already physically bad
- so the current face-only shell language is not rich enough to define a viable
  nested QW representation

In other words, the nested fixed block cannot be built from:

- shell faces
- plus one leftover interior/core block

and then stop there.

## Representation Requirement

The stronger requirement is:

every fixed-block basis function in the nested representation should belong to
some contraction object.

If a geometric sector is left as an omitted remainder, or is only represented
indirectly by a different codimension sector, the construction is not complete
enough to be reliable.

For a rectangular shell, that means the contraction language must eventually
cover:

- faces
- edges
- corners
- interior/core

## Complete Shell-layer Decomposition

For one rectangular box pair

- outer box `B_out = I_x^out x I_y^out x I_z^out`
- inner box `B_in  = I_x^in  x I_y^in  x I_z^in`

with `I_x^out`, `I_y^out`, `I_z^out` each one parent site larger on both sides
than the corresponding inner interval, the shell annulus

- `A = B_out \\ B_in`

has the natural disjoint decomposition by codimension.

### Face pieces

These are the open faces, where exactly one coordinate is on the outer
boundary and the other two remain in the inner interval:

- `x`-faces: `x = left/right`, `y in I_y^in`, `z in I_z^in`
- `y`-faces: `y = left/right`, `x in I_x^in`, `z in I_z^in`
- `z`-faces: `z = left/right`, `x in I_x^in`, `y in I_y^in`

Local primitive:

- two-dimensional tangential contraction
- implemented as `doside x doside`

This is what the current code already has.

### Edge pieces

These are the open edges, where exactly two coordinates are on the outer
boundary and the third remains in the inner interval:

- `xy`-edges: `x = left/right`, `y = left/right`, `z in I_z^in`
- `xz`-edges: `x = left/right`, `z = left/right`, `y in I_y^in`
- `yz`-edges: `y = left/right`, `z = left/right`, `x in I_x^in`

Local primitive:

- one-dimensional tangential contraction
- implemented as a single `doside` on the free edge interval

These are the first missing pieces beyond faces.

### Corner pieces

These are the shell corners, where all three coordinates are on the outer
boundary:

- `x = left/right`, `y = left/right`, `z = left/right`

Local primitive:

- zero-dimensional retained object
- the first practical version can simply keep one direct basis function per
  corner, or a tiny local corner block if later symmetry grouping is useful

### Interior/core piece

This is the retained interior block after the shell is removed.

Local primitive:

- direct core block
- or contracted core block
- but in either case it should itself be a declared contraction object, not an
  accidental leftover

## Disjoint-support Rule

The shell annulus should be partitioned by the number of boundary coordinates:

- codimension 1: faces
- codimension 2: edges
- codimension 3: corners

Those strata are disjoint by construction, and together with the retained
interior/core they cover the intended parent basis rows with no leftovers.

That is the right completeness rule for the nested fixed-block language.

## Diagnostic From The Current Face-only Failure

The current failing multi-shell test already shows the incompleteness in simple
counting terms.

For one shell annulus between a `13 x 13 x 13` box and an `11 x 11 x 11` box:

- total shell rows: `13^3 - 11^3 = 866`
- face rows only: `6 * 11^2 = 726`
- missing rows: `866 - 726 = 140`

That missing `140` is not mysterious. It is exactly:

- edges: `12 * 11 = 132`
- corners: `8`

So even one nominal shell already leaves a nontrivial codimension-2 and
codimension-3 remainder if only faces are used.

Likewise:

- `11^3 - 9^3 = 602 = 486` faces `+ 108` edges `+ 8` corners
- `9^3 - 7^3 = 386 = 294` faces `+ 84` edges `+ 8` corners

So the current face-only shell language is incomplete at every shell layer,
not just in the final tiny core.

## Architectural Conclusion

The next viable nested QW design target is not “better face compression.”
It is a complete contraction language for one shell layer:

- face objects
- edge objects
- corner objects
- retained interior/core object

all carried in the same fixed-block packet language the current consumer
already reads.

The consumer model does not need to change again. The missing work is on the
source side: the shell-layer representation itself is not yet complete.

## Smallest Next Implementation Step

The smallest useful next implementation step is:

1. add edge primitives to the existing face language
2. add direct corner pieces
3. assemble one complete nonrecursive shell layer from
   - faces
   - edges
   - corners
   - interior/core
4. propagate the same fixed-block packet through that complete shell layer
5. re-run the same nearest/GGT He check before opening recursion

That is the narrowest step that can test whether completeness, rather than just
compression strength, is the missing ingredient.
