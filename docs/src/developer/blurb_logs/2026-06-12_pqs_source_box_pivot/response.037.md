Pass 037 complete.

Summary:
- Tightened `_assert_pqs_source_record_retained_rule(...)` so retained PQS source records now require `retained_rule_kind === :boundary_comx_product_mode_selection`.
- Added the direct retained-boundary by-center nuclear kernel `_pqs_source_pair_direct_retained_electron_nuclear_by_center_block(...)`.
- Routed the active retained by-center wrappers through the direct retained-boundary path:
  - `pqs_source_pair_retained_electron_nuclear_by_center_block(...)`
  - `pqs_source_pair_retained_centered_electron_nuclear_by_center_block(...)`
- Preserved the existing raw source-result selector path `pqs_source_pair_retained_one_body_block(source_result, ...)` as oracle/fallback over an already materialized source block.
- Preserved the uncharged by-center convention: center key/index/location and nuclear charge are recorded, but charge is not applied and centers are not summed.

Files changed:
- `src/cartesian_pair_block_materialization/pqs_source_safe_terms.jl`
- `test/nested/cartesian_pair_block_materialization_contract_runtests.jl`

Behavior notes:
- The active retained electron-nuclear wrappers no longer materialize a full raw source-space electron-nuclear block before selecting retained columns.
- The centered wrapper still builds centered Gaussian 1D factor terms, then fills only the retained block.
- Metadata for the active direct retained nuclear result now reports retained direct-boundary construction and does not claim raw source operator materialization.
- Oracle comparisons against the old materialized-source selector path are exact for the existing CPBM fixture:
  - centered retained nuclear direct/oracle max difference: `0.0`
  - explicit by-center retained nuclear direct/oracle max difference: `0.0`

Validation:
- `julia --project=. test/nested/cartesian_pair_block_materialization_contract_runtests.jl` passed.
- `julia --project=. -e 'using GaussletBases; println("load ok")'` passed.
- `git diff --check` passed.

Deletion/shrinkage report:
- The active retained by-center nuclear route bypasses the old full raw source-space materialization stage.
- No source function was deleted because the raw source-result selector remains useful as the oracle/fallback path for explicitly materialized source results.
- No new test file was added. The existing CPBM contract test gained compact live-contract assertions proving direct retained construction and exact agreement with the oracle path.
- No compatibility shim was added.
- Remaining stale/duplicate surface: `pqs_source_pair_retained_one_body_block(source_result, ...)` still carries the raw materialized-source selector vocabulary by design, but it is now oracle/fallback rather than the active retained by-center wrapper path.

-- repo-doer@macmini
