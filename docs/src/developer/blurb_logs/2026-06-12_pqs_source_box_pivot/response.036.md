Pass 036 completed.

Files edited:

- `src/cartesian_pair_block_materialization/pqs_source_safe_terms.jl`
- `src/cartesian_pair_block_materialization/pqs_source_one_body.jl`
- `test/nested/cartesian_pair_block_materialization_contract_runtests.jl`
- `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.036.md`

Direct retained-boundary helpers added:

- `_pqs_source_pair_direct_retained_overlap_block`
- `_pqs_source_pair_direct_retained_kinetic_block`
- `_pqs_source_pair_direct_retained_result`
- `_pqs_source_record_retained_rule`
- `_assert_pqs_source_record_retained_rule`

Active retained overlap/kinetic path:

- `pqs_source_pair_retained_overlap_block(record; overlap_1d)` now builds the
  retained block directly from `CRPS.retained_mode_indices(rule)` and 1D
  overlap factors.
- `pqs_source_pair_retained_kinetic_block(record; overlap_1d, kinetic_1d)` now
  builds the retained block directly from retained mode tuples and 1D
  overlap/kinetic factors.
- `pqs_source_pair_retained_one_body_block(record, term; ...)` now dispatches
  retained overlap/kinetic to those direct retained-boundary wrappers.
- The existing `pqs_source_pair_retained_one_body_block(source_result, ...)`
  source-result selector path remains as oracle/fallback.

Metadata/comment update:

- Retained overlap/kinetic wrapper docstrings now describe direct retained
  boundary source-mode construction.
- Direct retained results record:
  - `retained_direct_boundary_product_used = true`;
  - `source_space_input_used = false`;
  - `raw_source_operator_block_materialized = false`;
  - no shell realization, Lowdin, final pair, Hamiltonian, IDA, RHF, driver,
    export, or artifact materialization.

Oracle comparison:

- CPBM contract test compares the new direct retained overlap against the old
  raw-source overlap then selector oracle.
  - `overlap_direct_oracle_maxdiff == 0.0`
- CPBM contract test compares the new direct retained kinetic against the old
  raw-source kinetic then selector oracle.
  - `kinetic_direct_oracle_maxdiff == 0.0`
- Retained shape/count remains `(26, 54)` in the existing `5 x 5 x 5` boundary
  fixture.

Validation run:

- `julia --project=. test/nested/cartesian_pair_block_materialization_contract_runtests.jl`
  passed.
- `julia --project=. -e 'using GaussletBases; println("load ok")'` passed.
- `git diff --check` passed.

Deletion/shrinkage report:

- The raw-source-block plus selector path did not become dead; it remains the
  explicit oracle/fallback path for already materialized source results and for
  non-overlap/kinetic retained source terms.
- Active retained overlap/kinetic no longer need to materialize full raw
  source-space blocks before boundary selection.
- Comments/docstrings were simplified so retained overlap/kinetic wrappers no
  longer claim raw source-space materialization.
- The test update adds compact live kernel coverage rather than broad metadata
  coverage: direct-vs-oracle equality, retained shape, and the key direct-path
  nonclaim flags.
- Remaining next optimization target: direct retained-boundary by-center
  nuclear construction from retained mode tuples and centered Gaussian source
  factors.

-- repo-doer@macmini
