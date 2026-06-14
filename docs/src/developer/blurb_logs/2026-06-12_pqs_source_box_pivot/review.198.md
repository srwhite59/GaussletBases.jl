Pass 198 review

Accepted.

This audit caught the key issue before the loop drifted into H1-J/RHF. The
current H2 route is a source-box/boundary diagnostic, not the intended
atom-core-plus-shell physics target.

Accepted inventory finding:

- Parent axis counts are `(x = 9, y = 9, z = 15)`, product `1215`.
- Source boxes are:
  - `pqs_left = 5 x 5 x 5`
  - `pqs_right = 5 x 5 x 5`
  - `product = 5 x 5 x 1`
- Retained counts are:
  - `pqs_left = 98`
  - `pqs_right = 98`
  - `product = 25`
- Dimension arithmetic is `98 + 98 + 25 = 221`.
- The atom-side raw boxes have full `5^3 = 125` support modes, but the retained
  atom-side final basis keeps only boundary modes: `125 - 27 = 98`.
- The 27 interior modes per atom box are not retained as atom-core bound-state
  content.
- There is no multishell atom-centered physics basis in the current H2 target.

Accepted interpretation:

- The positive H1 lowest value from pass 197 is consistent with this basis
  missing near-nuclear atom-core interior modes.
- It should be treated as a route-shape warning, not as a cue to continue into
  H1-J or RHF.
- The current 221-dimensional H2 artifact should be labeled as a route-smoke or
  source-box diagnostic only.

Accepted next-target recommendation:

- Choose option B: replace or extend the H2 driver target to a physical
  gausslet-only basis with retained full atom-core interiors plus shell layers
  before continuing to H1-J/RHF.

Deletion accepted:

- Deleted `test/nested/cartesian_pair_block_one_body_factor_inputs_runtests.jl`.
- The deleted file was private mixed one-body factor-input scaffold pressure,
  not a runner-included endpoint/reference test.

Validation reviewed:

- Deleted-file/testset search had no live `src/test/bin` hits.
- `julia --project=. -e 'using GaussletBases; println("load ok")'` passed.
- `git diff --check` passed.
- Scoped line budget is `0 added / 247 deleted`.

Next manager decision:

Pass 199 should label the current 221-dimensional H2 artifact as
`source_box_diagnostic` / not physics-endpoint-ready, preserving comparison
blocking. It should not yet implement the physical atom-core-plus-shell route.

-- repo-manager@macmini
