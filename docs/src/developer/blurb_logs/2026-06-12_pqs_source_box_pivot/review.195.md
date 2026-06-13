Pass 195 review

Accepted.

This pass advanced the H2 route one real dependency without crossing into H1:
the driver now carries real bond-aligned diatomic parent and axis-bundle
objects for the gausslet-only H2 readiness target.

Accepted implementation:

- Reused the existing explicit-core-spacing parent-axis probe path with
  carried objects.
- Updated the H2 input to use constructed parent axes rather than a manual
  parent-axis-count fixture.
- Added a compact `parent/*` artifact group for parent availability/type labels.
- Added a final-basis request guard so readiness mode does not cascade into
  final-basis/H1 materialization when `run_h1`, `run_h1_j`, and
  `run_private_rhf` are all false.
- Deleted one stale mixed one-body dispatch scaffold test.

Accepted H2 readiness state:

- Parent axis counts are now constructed as `(x = 9, y = 9, z = 15)`.
- That count product is `1215`, matching the old H2 WL/QW fixed-block row size,
  which is a useful sanity check for the H2 parent lattice.
- Parent basis, wrapped QW basis, and axis-bundle objects are available.
- Source plan is now available:
  `:available_pqs_diatomic_complete_core_shell_source_plan`.
- The route now blocks at the intended next layer:
  `:missing_diatomic_complete_core_shell_final_basis_consumer`.

Validation reviewed:

- Doer's focused H2 readiness test passed, 38/38, in 1m31s.
- I rechecked `julia --project=. -e 'using GaussletBases; println("load ok")'`.
- I rechecked `git diff --check`.
- Deleted-test search had no live `src/test/bin` hits.
- Scoped line budget is `28 added / 265 deleted`, net `-237`.

Deletion accepted:

- Deleted `test/nested/cartesian_pair_block_one_body_dispatch_runtests.jl`.
- The deleted file was metadata-only mixed one-body dispatch scaffolding, not an
  endpoint/reference test and not runner-included.

Next manager decision:

Pass 196 should materialize the H2 final basis only, using the now-available
source plan. It should not require H1/H1-J/RHF, should keep comparison blocked
against supplemented WL/QW references, and should preserve the explicit/manual
driver endpoint style.

-- repo-manager@macmini
