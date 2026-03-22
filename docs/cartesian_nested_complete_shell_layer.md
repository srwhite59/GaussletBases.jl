# Cartesian Nested Complete Shell Layer

This note records the first nonrecursive shell layer that is complete in the
sense described in [cartesian_nested_representation_completeness.md](/Users/srw/Library/CloudStorage/Dropbox/codexhome/repositories/GaussletBases/docs/cartesian_nested_representation_completeness.md).

This first complete-shell pass was run before the later shell-sequence
coverage fix. So its multi-shell complete-sequence numbers are now superseded
by
[cartesian_nested_sequence_coverage_fix.md](/Users/srw/Library/CloudStorage/Dropbox/codexhome/repositories/GaussletBases/docs/cartesian_nested_sequence_coverage_fix.md),
which adds the missing `7^3 - 5^3` annulus and reruns the same He comparison.

The current nested fixed-block consumer is unchanged. The only change in this
pass is on the source side: one shell layer is no longer represented by faces
alone.

## Complete Shell-layer Language

The first complete shell layer now uses the codimension partition:

- faces: open codimension-1 strata, represented by `doside x doside`
- edges: open codimension-2 strata, represented by one free-interval `doside`
- corners: codimension-3 strata, represented directly in the first pass
- interior/core: explicit retained interior block

The support partition is the one from the completeness note:

- faces = exactly one boundary coordinate
- edges = exactly two boundary coordinates
- corners = exactly three boundary coordinates
- interior = no boundary coordinates

So each shell annulus is covered with no leftovers before any compression
question is asked.

## First Complete Shell-layer Diagnostic

For the shell annulus between the `13^3` and `11^3` boxes:

- total shell rows: `866`
- face rows: `726`
- edge rows: `132`
- corner rows: `8`

The new complete shell layer includes all of those rows explicitly.

## First He nearest/GGT Comparison

On the stabilized He fixed-`a` count-17 case:

- baseline unnested:
  - fixed dim `4913`
  - overlap error `4.096500916261903e-12`
  - `E1 = -1.998588440029027`
  - `⟨Vee⟩ = 1.248965557065494`
  - total time `8.117796634 s`

- shell-plus-core:
  - fixed dim `1403`
  - overlap error `2.2354340600827527e-12`
  - `E1 = -1.9985629071304167`
  - `⟨Vee⟩ = 1.2489197930296585`
  - total time `3.381204923 s`

- face-only three-shell sequence:
  - fixed dim `341`
  - overlap error `1.208411649553688e-14`
  - `E1 = -1.6437270930603496`
  - `⟨Vee⟩ = 4.4067758554156065`
  - total time `1.206416259 s`

- complete three-shell sequence:
  - fixed dim `473`
  - overlap error `1.208411649553688e-14`
  - `E1 = -1.8438441482141197`
  - `⟨Vee⟩ = 4.175960247861663`
  - total time `1.374228679 s`

## Interpretation

This is the first shell language complete enough to test the real
representation question fairly.

The result is:

- yes, adding edges and corners improves the face-only shell-sequence physics
- yes, the size/time tradeoff remains promising relative to the unnested and
  shell-plus-core cases
- but no, the complete shell layer is still not physically good enough yet to
  justify opening recursion

So the architectural conclusion is:

- face-only incompleteness was a real problem
- the complete shell layer is the right source language
- but one more improvement is still needed before recursion:
  a better interior replacement policy or a richer low-dimensional retained
  content on the inner boundary hierarchy

The important point is that the consumer model still does not need to change.
The next step remains a source-side refinement in the same fixed-block packet
language.
