Ran the 5 x 5 x 5 PQS retained-source H1 readiness probe using the real
raw-product source-box vocabulary.

Artifacts:

- `tmp/work/pqs_source_box_5x5_h1_readiness_probe.jl`
- `tmp/work/pqs_source_box_5x5_h1_readiness_probe_output.txt`

Probe shape:

- Built a `CartesianRawProductSources.raw_product_box_plan(...)` with
  `source_mode_dims = (5, 5, 5)`.
- Built the `pqs_boundary_product_mode_retained_rule(...)` and confirmed the
  retained boundary count is `98`.
- Attached materialized identity source-axis transform facts, explicitly
  labeled as `:diagnostic_identity_source_mode_transform`.
- Built retained-source `S`, `T`, and uncharged by-center `V_unit` through the
  existing CPBM helpers.
- Formed `H = T + Z * V_unit` only inside the probe.
- No production code, tests, shell realization, Lowdin cleanup, IDA,
  density-density, RHF, global driver, exports, artifacts, full-parent
  fallback, direct Cartesian fallback, or CCPM wrapper path was added.

Key diagnostics:

- raw source mode dims/count: `(5, 5, 5)` / `125`
- retained boundary count/dimension: `98`
- transform label: `:diagnostic_identity_source_mode_transform`
- transform statuses/shapes:
  `(:materialized, :materialized, :materialized)` /
  `((5, 5), (5, 5), (5, 5))`
- `S`, `T`, `V_unit`, and `H` shapes: `(98, 98)`
- finite entries: true for all four matrices
- symmetry error: `0.0` for all four matrices
- overlap eigenvalue range:
  `0.9999999999999988` to `1.000000000000004`
- overlap rank: `98`
- overlap condition estimate: `1.0000000000000053`
- overlap identity error: `1.9984014443252818e-15`
- solve kind: `ordinary_retained_source_diagnostic`
- lowest diagnostic H1 energy: `1.2657907439492142`
- nuclear charge recorded: `true`
- nuclear charge applied in source block: `false`
- centers summed: `false`
- shell realization materialized: `false`
- Lowdin cleanup used: `false`
- IDA data materialized: `false`
- Hamiltonian data materialized in source route: `false`

Readiness result and blocker:

- The real raw-product source-box retained-rule vocabulary works for the
  retained-source H1 wiring probe.
- This run is diagnostic only because the source-axis transforms are identity
  matrices supplied by the probe.
- Exact blocker for physical/final PQS H1:
  `:missing_pqs_source_axis_transform_builder`.
- The ordinary solve is not a physical acceptance result; it is an
  identity-source-mode retained-source diagnostic before real PQS source-axis
  transforms, shell realization, and Lowdin/final-basis construction.

Validation:

- `julia --project=. tmp/work/pqs_source_box_5x5_h1_readiness_probe.jl`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git diff --check`

Deletion/shrinkage report:

- No old code, test, metadata, or compatibility path became unnecessary in this
  probe-only pass.
- Nothing was deleted or simplified because the pass intentionally avoided
  production changes and used ignored `tmp/work` artifacts.
- No new long-term test was added; the probe artifact earns its cost by
  documenting the first 5 x 5 x 5 raw-product retained-source H1 readiness
  check and the precise source-axis transform blocker.
- Remaining stale/duplicate surfaces to retire next: old CCPM source-box H1 or
  nuclear wrappers remain quarantine candidates once a real PQS source-axis
  transform builder feeds this retained-source path.

-- repo-doer@macmini
