Status: complete.

Performed the targeted PQS final-basis H1 cleanup/retirement audit.

Files changed:

- `docs/src/developer/pqs_source_box_operator_framework.md`
- `docs/src/developer/raw_product_source_retained_transform_policy.md`
- `docs/src/developer/cartesian_coordinate_product_box_contract.md`
- `docs/src/developer/cartesian_route_dictionary.md`
- `docs/src/developer/cartesian_route_retirement_ledger.md`

What was simplified/quarantined:

- Updated source-box/operator docs so `:shell_realization_not_materialized` is explicitly scoped to bridge-level source-space readiness summaries, not the current status of the final-basis H1 route.
- Updated route dictionary/contract wording to distinguish:
  - metadata-only bridge readiness, which still correctly blocks on `:shell_realization_not_materialized`;
  - the explicit PQS final-basis H1 path, which now consumes retained boundary overlap/kinetic, separated by-center nuclear matrices, shell realization, Hamiltonian-stage charge application, and an ordinary final-basis eigensolve.
- Added a short retirement-ledger note recording that the source-box/final-basis H1 probe now succeeds and that `_pqs_current_route_safe_term_matrices(...)` remains oracle/debug vocabulary, not the production PQS final-basis route.

Audit result:

- Active CPBM tests that still assert `:shell_realization_not_materialized` were left intact. They protect the bridge-readiness contract and are not preserving the old H1 scratch route as a live path.
- Raw source-box block tests, retained rule tests, retained source one-body matrix tests, and centered nuclear source tests were left intact because they protect live kernels/seams now consumed by the final-basis route.
- `_pqs_current_route_safe_term_matrices(...)` source and its focused integration usage were not deleted in this pass. The docs now quarantine it as private diagnostic/oracle scaffolding. Deleting it safely needs a caller-driven pass over the current-route authority comparison tests.

Remaining stale/duplicate surfaces:

- `src/CartesianContractedParentMetrics.jl` still contains `_pqs_current_route_safe_term_matrices(...)` and related authority-comparison helpers.
- `test/nested/bond_aligned_diatomic_high_order_recipe_opt_in_source_construction_integration_runtests.jl` still calls `_pqs_current_route_safe_term_matrices(...)` as part of current-route oracle coverage.
- Historical tmp/work retained-source H1 probes remain ignored artifacts only.
- The bridge-level `pqs_source_pair_final_block_readiness_summary(...)` vocabulary remains valid for source-space bridge summaries, but should not be used to describe final-basis H1 readiness.

Validation:

- `julia --project=. -e 'using GaussletBases; println("load ok")'` passed.
- `git diff --check` passed.
- No Julia test was run because only docs/ledger files were edited.

Deletion/shrinkage report:

- No production code or tests were deleted in this pass.
- Active doc wording was simplified to remove the conceptual drift that bridge readiness equals final-basis H1 readiness.
- The new retirement-ledger note earned its cost by recording the post-H1-probe contract boundary and the next retirement target.
- No new tests or metadata fields were added.
- Next deletion-oriented pass should audit `_pqs_current_route_safe_term_matrices(...)` callers and shrink/delete current-route safe-term oracle tests that no longer protect live behavior.

-- repo-doer@macmini
