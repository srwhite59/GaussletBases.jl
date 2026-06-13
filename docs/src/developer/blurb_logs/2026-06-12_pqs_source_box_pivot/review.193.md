Pass 193 review

Accepted.

This pass did the right first H2 move: it added a driver-visible gausslet-only
H2 readiness artifact without pretending that the old supplemented WL/QW H2
HF/ED references are comparable.

Accepted implementation:

- Added `test/driver_inputs/h2_pqs_q5_gausslet_only_r4.jl`.
- Added an explicit/manual readiness test:
  `test/nested/cartesian_ham_builder_h2_pqs_q5_gausslet_only_r4_readiness_runtests.jl`.
- Kept the test out of default/integration runners.
- Extended the visible driver inputs for bond metadata, extents, supplement
  policy, comparison readiness, and H1/H1-J toggles without replacing the
  staged driver script with an opaque wrapper.
- Added a blocker-aware diatomic readiness artifact writer for two-center PQS
  with `supplement_policy = :none`.
- Deleted two stale mixed one-body metadata scaffold tests.

Accepted artifact behavior:

- Records H2 geometry:
  `(0,0,-2)` and `(0,0,2)`, charges `(1,1)`, bond axis `:z`, bond length `4.0`.
- Records `q = 5`, `n_s = 5`, `supplement_policy = :none`.
- Records `comparison_ready = false` and blocker
  `:supplemented_reference_not_comparable_to_gausslet_only`.
- Does not store `comparison/wl_rhf_total` or `comparison/delta_rhf`.
- Does not claim private RHF, H1-J, public API, export, or route-global behavior.
- Records the live H2 materialization blocker:
  `:missing_diatomic_complete_core_shell_source_plan_producer`.

Validation reviewed:

- Focused H2 readiness test passed, 29/29, in 43.3s on my rerun.
- `julia --project=. -e 'using GaussletBases; println("load ok")'` passed.
- `git diff --check` passed.
- Deleted-test search had no live `src/test/bin` hits.
- Scoped budget reported by tracked diff is `237 added / 486 deleted`.
- Including new untracked files, the real scoped budget is `341 added / 486
  deleted`, net `-145`.

Deletion accepted:

- Deleted `test/nested/cartesian_pair_block_one_body_block_set_preflight_runtests.jl`.
- Deleted `test/nested/cartesian_pair_block_one_body_batch_summary_runtests.jl`.
- Both were metadata-only mixed one-body scaffold pressure, not physics
  endpoints and not runner-included coverage.

Watchpoint for the next pass:

The pass adds a fairly broad set of diatomic readiness report aliases. Do not
build another physics layer on top of that surface by inertia. The next pass
should audit the live blocker
`:missing_diatomic_complete_core_shell_source_plan_producer` and identify the
smallest route-owned producer path, preferably consuming compact payloads rather
than adding another report-field cloud.

-- repo-manager@macmini
