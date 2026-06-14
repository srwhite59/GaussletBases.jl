Pass 200 review

Accepted.

This audit answers the sizing question that blocked H2. The old WL/QW
gausslet-only H2 route is not the current 221-column source-box diagnostic. The
physical target should be an atom-contact core plus shared shell layers.

Accepted old WL inventory:

- Parent axis lengths are `(9, 9, 15)`, parent size `1215`.
- Old gausslet-only fixed block is `(1215, 463)`.
- Supplemented reference final dimension is `481`, from `463 + 18` residual
  H/cc-pVTZ S/P columns.
- Retained order is core/child first, then shared shell layers.
- Core/child columns are `1:251`.
- Shared shell layer columns are `252:349` and `350:463`, counts `(98, 114)`.
- Shared shell support counts are `(578, 362)`.
- The child/core support count is `275`.
- The child/core support is the `5 x 5 x 11` atom-contact working box:
  two full `5^3` atom cores plus the 25-row contact plane.
- There is no separate midpoint/product retained unit in the old 463-column
  gausslet-only fixed block.

Accepted interpretation:

- The current 221-column route is boundary-only: `98 + 98 + 25`.
- The 25-row contact plane should belong inside an atom-contact core for the
  physical H2 target, not remain a standalone retained midpoint slab.
- Full atom-core interiors must be retained through the atom-contact core.
- Shell layers should be source-box-first PQS shell units over the old shared
  molecular shell supports.

Accepted next seam:

- Add a new route kind rather than mutating the existing diagnostic route.
- Suggested route kind is
  `:bond_aligned_diatomic_physical_gausslet_core_shell_pqs`.
- First code pass should add only a compact private target/inventory payload.
- Do not implement H1-J, RHF, supplement handling, or comparison in that pass.

Deletion accepted:

- Deleted `test/nested/cartesian_pair_block_one_body_plan_batch_runtests.jl`.
- This was private mixed one-body batch scaffold pressure, not a physics
  endpoint, reference, or runner-included test.
- The compact live mixed one-body consumer smoke remains.

Validation reviewed:

- Deleted-file/testset search had no live `src/test/bin` hits.
- `julia --project=. -e 'using GaussletBases; println("load ok")'` passed.
- `git diff --check` passed.
- Scoped source/test/bin budget is `0 added / 480 deleted`.

Next manager decision:

Pass 201 should implement the private physical H2 PQS target/inventory payload
and a minimal driver artifact/readiness surface for that new route kind. It
should not build source plans, final basis, H1, H1-J, RHF, supplement, or a WL
comparison yet.

-- repo-manager@macmini
