Pass 197 review

Accepted as an H1 diagnostic milestone, not as physics acceptance.

This pass materialized H2 H1 through the existing diatomic complete-core/shell
path and kept the requested trust boundary: no H1-J, no private RHF, no
supplemented WL/QW comparison, no exports, and no public solver behavior.

Accepted implementation:

- H2 input now requests `run_h1 = true` while leaving `run_h1_j = false` and
  `run_private_rhf = false`.
- H1 uses `_pqs_source_box_route_driver_diatomic_complete_core_shell_h1_payload`.
- Artifact records compact H1 diagnostics only:
  `physics/h1_lowest`, H1 matrix finite flag, and H1 symmetry error.
- No broad H1 matrix is saved.
- Deleted one stale mixed one-body record-dispatch scaffold test.

Accepted H2 H1 diagnostic state:

- Final dimension remains `221`.
- Final overlap identity error remains about `2.51e-13`.
- H1 status is `:materialized_pqs_complete_core_shell_final_h1_solve`.
- `physics/h1_lowest = 0.14582426982296057`.
- H1 matrix finite flag is true.
- H1 symmetry error is about `3.66e-15`.
- H1-J and private RHF remain false.
- The next recorded blocker is
  `:missing_diatomic_complete_core_shell_hamiltonian_handoff_payload`.

Important watchpoint:

The positive H1 lowest value is a physics warning, not a success criterion.
For a real H2 one-body Hamiltonian one would expect bound negative one-electron
states. The current 221-dimensional H2 source plan is therefore likely still a
route/source-box diagnostic basis, not the intended atom-core-plus-multishell
physics endpoint. Do not proceed to H1-J or RHF before auditing the H2 source
plan/final-basis inventory against the intended physics target.

Validation reviewed:

- Doer's focused H2 readiness test passed, 47/47, about 1m33s.
- Doer reported the runtime remains compilation dominated.
- I rechecked `julia --project=. -e 'using GaussletBases; println("load ok")'`.
- Deleted-test search had no live `src/test/bin` hits.
- Scoped line budget is `29 added / 394 deleted`, net `-365`.

Deletion accepted:

- Deleted `test/nested/cartesian_pair_block_one_body_record_dispatch_runtests.jl`.
- The deleted file was record-level mixed one-body dispatch scaffolding, not an
  endpoint/reference test.

Next manager decision:

Pass 198 should be an H2 source-plan/final-basis inventory audit, not H1-J.
It should answer whether the current 221-dimensional basis lacks full atom-core
interior modes and/or multishell content, and it should decide what driver
target is needed for the intended H2 physics endpoint.

-- repo-manager@macmini
