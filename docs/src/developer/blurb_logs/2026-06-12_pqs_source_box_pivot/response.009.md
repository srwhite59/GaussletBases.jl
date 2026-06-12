Pass 009 response: PQS source-axis transform fact plumbing

Status: implemented.

Files changed:

- `src/cartesian_raw_product_sources/CartesianRawProductSources.jl`
- `src/cartesian_raw_product_sources/records.jl`
- `src/cartesian_raw_product_sources/axis_transform_facts.jl`
- `src/cartesian_retained_unit_transform_contracts/unit_contracts.jl`
- `src/cartesian_pair_block_materialization/preflight.jl`
- `test/nested/cartesian_raw_product_sources_contract_runtests.jl`
- `test/nested/cartesian_retained_unit_transform_contracts_runtests.jl`
- `test/nested/cartesian_pair_block_materialization_contract_runtests.jl`

What changed:

- Added `CartesianRawProductSources.axis_source_transform_fact(...)`.
- Added `raw_product_box_plan(...; axis_transform_matrices = ..., axis_transform_facts = ...)`.
- Preserved the default path: plans without supplied transforms still carry
  three `:not_materialized` `AxisSourceTransformFact`s with
  `coefficient_matrix = nothing`.
- Validated materialized coefficient matrices against:
  - `length(source_interval)` rows;
  - `source_mode_dim` columns.
- Allowed supplied transforms as either:
  - `axis_transform_matrices = (x = ..., y = ..., z = ...)` or a 3-tuple;
  - `axis_transform_facts = (...)` after validating axes, intervals,
    dimensions, and materialization status.
- Threaded PQS retained-unit metadata
  `raw_product_source_axis_transform_matrices` or
  `raw_product_source_axis_transform_facts` through
  `CartesianRetainedUnitTransformContracts` into `RawProductBoxPlan`.
- Added `raw_product_source_axis_transform_facts` to the retained-transform
  contract metadata.
- Exposed compact left/right
  `*_raw_product_source_axis_transform_facts` on CPBM PQS source-pair preflight
  records. This carries the existing fact tuple rather than copying a flat
  metadata field cloud.

Axis-transform plumbing evidence:

- Default raw product source plans still report axis transform statuses
  `(:not_materialized, :not_materialized, :not_materialized)`.
- Explicit matrix input now reports
  `(:materialized, :materialized, :materialized)` and stores validated matrix
  shapes such as `(2, 3)`, `(3, 4)`, `(4, 5)`.
- PQS transform contracts now carry materialized axis transform facts when the
  retained unit metadata supplies them.
- CPBM source-pair preflight records expose left/right axis transform facts with
  the expected source support-row by source-mode shapes.
- No Gaussian factor generation, nuclear block generation, shell realization,
  Lowdin cleanup, IDA, Hamiltonian assembly, driver adoption, exports, or
  artifacts were added.

Validation:

- `julia --project=. test/nested/cartesian_raw_product_sources_contract_runtests.jl` passed.
- `julia --project=. test/nested/cartesian_retained_unit_transform_contracts_runtests.jl` passed.
- `julia --project=. test/nested/cartesian_pair_block_materialization_contract_runtests.jl` passed.
- `julia --project=. -e 'using GaussletBases; println("load ok")'` passed.
- `git diff --check` passed.

Deletion/shrinkage report:

- No old code, test, metadata, or compatibility path became unnecessary yet.
  This pass made the existing `AxisSourceTransformFact` placeholder usable
  rather than replacing a live older implementation.
- Nothing was deleted or simplified.
- New tests are genuinely new live-contract coverage for the newly active
  source-axis transform seam:
  - default raw-source facts stay not materialized;
  - explicit source-axis matrices validate and become materialized facts;
  - retained-transform contracts and CPBM preflight records carry those facts.
- Remaining stale/duplicate surfaces to retire next:
  - CCPM centered local-Gaussian source-box wrappers remain oracle/reference
    only until CPBM consumes these facts to generate PQS source-mode Gaussian
    factor arrays directly.
  - The next pass can now target the source-axis Gaussian factor helper without
    using CCPM physical/nuclear wrappers.

-- repo-doer@macmini
