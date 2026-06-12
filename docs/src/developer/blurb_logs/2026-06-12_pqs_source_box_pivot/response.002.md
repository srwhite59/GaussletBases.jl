Implemented the first compact PQS source-mode retained boundary rule and attached it to PQS retained-unit transform-contract metadata.

Files changed:

- `src/cartesian_raw_product_sources/CartesianRawProductSources.jl`
- `src/cartesian_raw_product_sources/records.jl`
- `src/cartesian_raw_product_sources/summaries.jl`
- `src/cartesian_retained_unit_transform_contracts/unit_contracts.jl`
- `test/nested/cartesian_raw_product_sources_contract_runtests.jl`
- `test/nested/cartesian_retained_unit_transform_contracts_runtests.jl`

What changed:

- Added `PQSBoundaryProductModeRetainedRule` under `CartesianRawProductSources`.
- Added `pqs_boundary_product_mode_retained_rule(...)`.
- Added retained-rule accessors:
  - `retained_mode_indices(rule)`
  - `retained_column_indices(rule)`
- Added compact `summary(rule::PQSBoundaryProductModeRetainedRule)`.
- Attached `raw_product_source_retained_rule` and `raw_product_source_retained_rule_summary` in `_raw_product_source_contract_metadata(unit)` for `:pqs_shell_retained_unit`.
- Missing raw-source metadata now carries `nothing` for the retained-rule fields, rather than implying a rule exists.

Retained-rule shape:

- `source_key`
- `source_mode_dims`
- `source_mode_ordering`
- `retained_rule_kind = :boundary_comx_product_mode_selection`
- `retained_mode_indices`
- `retained_column_indices`
- `retained_count`
- `transform_kind = :source_mode_column_selector`
- `shell_realization_materialized = false`
- `lowdin_cleanup_used = false`
- `metadata`

Evidence:

- `(5, 5, 5)` source modes retain `98` boundary product modes.
- The selector keeps deterministic source-mode order and retained source columns.
- The retained modes are exactly those with at least one local axis index on a boundary (`1` or max axis index).
- PQS retained-unit transform-contract metadata now carries the rule summary; for the existing `(3, 3, 3)` PQS fixture, retained count is `26`.
- The explicit `(5, 4, 3)` PQS raw-source fixture retains `54` boundary modes.

Validation run:

- `julia --project=. test/nested/cartesian_raw_product_sources_contract_runtests.jl` passed.
- `julia --project=. test/nested/cartesian_retained_unit_transform_contracts_runtests.jl` passed.
- `julia --project=. -e 'using GaussletBases; println("load ok")'` passed with `load ok`.
- `git diff --check` passed.

Deletion/shrinkage report:

- No old code was deleted in this pass.
- No old compatibility path became obsolete yet because this is the first source-mode retained-rule object; no pair-block contraction or retained-source operator path consumes it yet.
- The new test coverage is genuinely new live-contract coverage for the source-box-first PQS pivot: it protects the retained source-mode boundary rule and its transform-contract attachment.
- No older test was replaced or shrunk.
- Remaining stale/duplicate surfaces to retire next:
  - old CCPM/source-box reference helpers that still preserve shell/support-row vocabulary should remain oracle-only until retained source-mode contraction exists, then shrink/delete where unused;
  - old raw-box selector metadata paths can be retired or narrowed once the new rule becomes the active retained source-mode authority for pair-block contraction.

-- repo-doer@macmini
