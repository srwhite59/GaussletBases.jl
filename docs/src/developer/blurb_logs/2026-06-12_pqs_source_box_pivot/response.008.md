Pass 008 response: PQS source-axis Gaussian factor source audit

Status: blocked, no implementation added.

Exact blocker:

```julia
:missing_pqs_source_axis_gaussian_factor_source
```

Inspected surfaces:

- `src/cartesian_raw_product_sources/CartesianRawProductSources.jl`
- `src/cartesian_raw_product_sources/records.jl`
- `src/cartesian_raw_product_sources/axis_transform_facts.jl`
- `src/cartesian_retained_unit_transform_contracts/unit_contracts.jl`
- `src/cartesian_retained_unit_transform_contracts/records.jl`
- `src/cartesian_pair_block_materialization/preflight.jl`
- `src/cartesian_pair_block_materialization/pqs_source_safe_terms.jl`
- `src/CartesianParentAxisFactors.jl`
- `src/CartesianContractedParentMetrics.jl`
- `src/CartesianCPBBlockProviders.jl`
- `src/cartesian_retained_units/lower_contract_units.jl`
- `src/cartesian_terminal_lowering/region_contracts.jl`
- `docs/src/developer/pqs_source_box_operator_framework.md`

Audit result:

- The live path for the blurb's `src/CartesianRawProductSources.jl` target is
  `src/cartesian_raw_product_sources/CartesianRawProductSources.jl`.
- `RawProductBoxPlan` currently owns the filled/source CPB, source intervals,
  source-mode dimensions, source-mode ordering, source-mode indices, and
  metadata-only `AxisSourceTransformFact`s.
- The default axis transform facts explicitly set
  `coefficient_status = :not_materialized` and `coefficient_matrix = nothing`.
  They do not carry source-axis basis/backend objects or source-mode coefficient
  maps.
- `CartesianRetainedUnitTransformContracts` is explicitly metadata-only. For
  PQS units it reconstructs a raw-product source plan from the single filled
  source CPB and source-mode dimensions, then attaches only the boundary
  retained rule and summaries.
- `preflight.jl` records source-mode dimensions, counts, ordering, retained
  rules, and transform/source contract keys. It does not attach source-axis
  layers, primitive/backend data, or coefficient matrices to the
  `PairBlockMaterializationRecord`.
- Pass 007's `pqs_source_pair_electron_nuclear_by_center_block(...)` is
  correctly shaped for caller-supplied term-first arrays, but it has no current
  route-owned way to build those arrays from the record.
- `CartesianParentAxisFactors` and `CartesianCPBBlockProviders` have
  parent-axis PGDG/electron-nuclear factor conventions, but those are indexed
  by parent axis data. The PQS source-pair record does not carry the source-axis
  transform data needed to project or reinterpret them as source-mode arrays of
  shape `(nterms, left_axis_count, right_axis_count)`.
- `CartesianContractedParentMetrics` still contains private centered
  local-Gaussian source-box wrappers using `gaussian_factor_matrices(...)`, but
  those wrappers require explicit `axis_layers` and are CCPM oracle/reference
  infrastructure. The blurb correctly forbids routing the new CPBM helper
  through `_pqs_pqs_source_box_*nuclear_attraction*` wrappers.

Decision:

I did not implement a Gaussian factor source helper because the active CPBM
PQS source-pair route currently carries dimensions/order/retained rules only.
Producing factor arrays would require inventing or silently recovering missing
source-axis representation data. That would add late binding at exactly the
boundary this pass is trying to clarify.

Recommended next seam:

- Add a narrow route-owned PQS source-axis factor source object at the raw
  product source or retained-transform boundary.
- It should carry, per axis, the source-mode coefficient/axis representation
  needed to project analytic `gaussian_factor_matrices(...)` output into
  source-mode coordinates.
- The object should be constructed from explicit source-axis/parent-axis
  bundle facts, not from CCPM physical wrappers, and should then feed the pass
  007 CPBM electron-nuclear source block.

Validation:

- `julia --project=. -e 'using GaussletBases; println("load ok")'` passed.
- `git diff --check` passed.
- No production code or tests changed, so the CPBM contract test was not run.

Deletion/shrinkage report:

- No old code, test, metadata, or compatibility path became unnecessary in this
  pass because no replacement source-axis factor object exists yet.
- Nothing was deleted or simplified.
- No new test or artifact was added; this was a blocker audit and the blurb says
  not to add tests when blocked.
- Remaining stale/duplicate surfaces to retire later:
  - CCPM centered local-Gaussian source-box wrappers remain oracle/reference
    only until CPBM owns a real source-axis Gaussian factor source.
  - `RawProductBoxPlan` documentation says future adapters may attach
    externally built axis transforms, while the live plan still has
    `:not_materialized` axis facts. That is the next contract gap to close.

-- repo-doer@macmini
