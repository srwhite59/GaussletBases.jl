Pass 199 review

Accepted.

The current H2 221-dimensional artifact now carries the warning labels it
needed. This prevents the diagnostic route from being mistaken for a physics
endpoint.

Accepted labels:

- `route/artifact_role = :source_box_diagnostic`
- `physics/endpoint_ready = false`
- `physics/endpoint_blocker = :retained_atom_core_interiors_missing`
- `basis/retained_atom_core_interiors = false`
- `basis/source_plan_role = :boundary_source_box_diagnostic`

Accepted behavior:

- The current 221-dimensional construction is unchanged.
- H1-J and private RHF remain false.
- Comparison to supplemented WL/QW H2 remains blocked.
- The H2 readiness test asserts the new diagnostic-only labels.

Validation reviewed:

- Doer's focused H2 readiness test passed, 52/52, about 1m35s on the warm
  command and compilation dominated.
- I rechecked `julia --project=. -e 'using GaussletBases; println("load ok")'`.
- I rechecked `git diff --check`.
- Deleted-test search had no live `src/test/bin` hits.
- Scoped line budget is `41 added / 302 deleted`, net `-261`.

Deletion accepted:

- Deleted `test/nested/cartesian_pair_block_one_body_placement_plan_runtests.jl`.
- The deleted file was private one-body placement-plan scaffold pressure, not a
  physics endpoint or reference test.

Next manager decision:

Pass 200 should not continue this 221-dimensional route into H1-J. It should
define the physical H2 gausslet-only driver target: full retained atom-core
interiors plus shell layers, with explicit source-plan/final-basis semantics.

-- repo-manager@macmini
