Purpose:

Plumb externally materialized PQS source-axis transforms to the pair-block
preflight boundary.

Why now:

Pass 008 showed that PQS source-mode Gaussian factors cannot be generated from
the current pair-block record because the route carries source dimensions and
retained rules, but not source-axis transform data. `RawProductBoxPlan` already
has `AxisSourceTransformFact`, and its docs say future adapters may attach
externally built axis transforms. Make that placeholder usable in the narrowest
way.

Exact task:

Add a compact path for externally supplied source-axis coefficient transforms:

```text
axis source support rows x source modes
```

Recommended scope:

1. In `CartesianRawProductSources`, add a small constructor or keyword path for
   materialized `AxisSourceTransformFact`s.
2. Validate each coefficient matrix shape against:
   - `length(source_interval)` rows;
   - `source_mode_dim` columns.
3. Preserve the default `:not_materialized` facts when no matrices are supplied.
4. In `CartesianRetainedUnitTransformContracts`, if a PQS retained unit carries
   explicit source-axis transform facts or matrices in metadata, pass them into
   the raw product source plan.
5. In `CartesianPairBlockMaterialization` preflight metadata, expose enough of
   the left/right raw product source axis transform facts for the next Gaussian
   factor helper to consume them. Prefer carrying the compact object/facts over
   copying a flat field cloud.

If the retained-unit metadata seam for explicit axis transforms does not exist
and adding it would require broad retained-unit redesign, do not fake it.
Report:

```julia
:missing_pqs_retained_unit_axis_transform_metadata_seam
```

Trust boundary:

No Gaussian factor generation, electron-nuclear H1 diagnostic, shell
realization, Lowdin cleanup, density-density, IDA, Hamiltonian assembly, driver
adoption, exports, artifacts, full-parent fallback, direct Cartesian fallback,
or old CCPM wrapper adoption.

Test policy:

Add only compact module-contract checks:

- default raw product plans still have three `:not_materialized` axis facts;
- a plan with explicit axis transform matrices records materialized facts and
  validates shapes;
- a PQS transform contract/preflight record can carry or expose those facts
  without materializing shell realization or Lowdin.

Do not add a physics test. Do not add broad metadata checks.

Deletion/shrinkage report required:

- what old code, test, metadata, or compatibility path became unnecessary;
- what was deleted or simplified;
- if nothing was deleted, why no existing surface was made obsolete yet;
- whether any new test/artifact was added and why it earned its cost;
- any remaining stale or duplicate surfaces to retire next.

Validation:

- `julia --project=. test/nested/cartesian_raw_product_sources_contract_runtests.jl`
- `julia --project=. test/nested/cartesian_retained_unit_transform_contracts_runtests.jl`
- `julia --project=. test/nested/cartesian_pair_block_materialization_contract_runtests.jl`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git diff --check`

If the pass blocks without production changes, only load check and diff check
are required.

Report back:

- write `.agent_handoffs/response.009.md.tmp`, then atomically rename to
  `.agent_handoffs/response.009.md`;
- also write the curated copy to
  `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.009.md`;
- include files changed or exact blocker;
- include axis-transform plumbing evidence;
- include validation run;
- include deletion/shrinkage report;
- sign `-- repo-doer@macmini`.

-- repo-manager@macmini
