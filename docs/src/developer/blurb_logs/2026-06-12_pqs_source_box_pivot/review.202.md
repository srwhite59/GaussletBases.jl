Pass 202 review

Accepted.

The private global-overlap driver hook was stale enough to delete rather than
extract. This is the correct outcome for the cleanup pass: the human-facing
driver helper no longer owns that optional, default-off global-overlap
placement subtree.

Accepted classification:

- Live current target: no.
- Production/public/default behavior: no.
- Driver option default before deletion: off.
- Endpoint/reference tests: none for the private driver hook.
- Private scaffold tests:
  - `test/nested/cartesian_cpb_local_overlap_fingerprint_runtests.jl`
  - `test/nested/cartesian_pair_block_driver_global_overlap_runtests.jl`
- Exact source callers before deletion were driver option plumbing and the
  private materialization stage.

Accepted deletions:

- Removed the private global-overlap helper subtree from
  `src/pqs_source_box_route_driver_helpers.jl`.
- Removed private global-overlap materialization/result/summary plumbing.
- Removed default-off knobs from `bin/cartesian_ham_builder.jl`.
- Deleted `examples/private_global_overlap_option.jl`.
- Deleted the two private scaffold tests above.
- Updated docs to mark the private driver bridge as retired.

Validation reviewed and rerun:

- `rg -n "private_global_overlap" src test bin examples` had no matches.
- Deleted-helper/test/example search had no matches.
- `julia --project=. -e 'using GaussletBases; println("load ok")'` passed.
- `julia --project=. bin/cartesian_ham_builder.jl save_artifact=false save_tsv=false`
  passed and printed `driver complete`.
- `git diff --check` passed.
- Scoped source/test/bin line budget is `0 added / 4510 deleted`.
- Overall diff stat is `23 insertions / 4564 deletions`.

Deletion accepted:

- Deleted private driver-global-overlap code and its preservation tests.
- Preserved lower route-global one-body/oracle paths because they are separate
  lower-layer materialization/reference pressure, not the deleted driver hook.

Remaining blocker:

- No remaining `src/test/bin/examples` caller of the deleted private driver hook.
- Physical H2 remains blocked only on its separate source-plan/final-basis/H1
  work for `:atom_contact_core_plus_pqs_shared_shells`.

-- repo-manager@macmini
