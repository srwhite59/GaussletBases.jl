Pass 196 review

Accepted.

This pass materialized the H2 final basis without using H1 as a backdoor. That
is the right layer boundary.

Accepted implementation:

- Added explicit visible driver input `run_final_basis`.
- Default behavior resolves `run_final_basis = nothing` from existing downstream
  requests, so existing H1/H1-J/RHF paths still request final basis as before.
- H2 input now sets `run_final_basis = true` while keeping `run_h1 = false`,
  `run_h1_j = false`, and `run_private_rhf = false`.
- Reused the existing private diatomic final-basis payload path.
- Added a separate H1 request gate so final-basis-only H2 does not solve H1.
- Deleted one stale plan-level mixed one-body metadata scaffold test.

Accepted H2 state:

- Source plan remains available.
- Final basis is available.
- `basis/final_dimension = 221`.
- `basis/final_overlap_identity_error = 2.510171831978689e-13`.
- H1 remains not materialized.
- H1-J and private RHF remain false.
- Comparison to supplemented WL/QW references remains disabled.
- The next blocker is now
  `:missing_diatomic_complete_core_shell_h1_consumer`.

Validation reviewed:

- Doer's focused H2 readiness test passed, 41/41, in 2m30.9s. The reported
  runtime was compilation dominated.
- Doer also ran a focused readback confirming final dimension, overlap identity
  error, and the H1-consumer blocker.
- I rechecked `julia --project=. -e 'using GaussletBases; println("load ok")'`.
- I inspected the diff and deletion.
- Deleted-test search had no live `src/test/bin` hits.
- Scoped line budget is `35 added / 352 deleted`, net `-317`.

Deletion accepted:

- Deleted `test/nested/cartesian_pair_block_one_body_plan_dispatch_runtests.jl`.
- The deleted file was plan-level metadata dispatch scaffolding, not an endpoint
  or reference test.

Next manager decision:

Pass 197 should materialize H2 H1 only through the driver artifact path. It
should assert finite/symmetric H1 and the lowest H1 eigenvalue, but still no
H1-J, density interaction, RHF, supplemented comparison, exports, or public
solver behavior.

-- repo-manager@macmini
