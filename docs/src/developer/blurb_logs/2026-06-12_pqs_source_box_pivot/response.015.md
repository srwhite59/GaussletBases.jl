Ran the 5 x 5 x 5 PQS retained-source H1 probe with non-identity source-axis
transforms supplied by the old nested helper as an oracle.

Artifacts:

- `tmp/work/pqs_source_box_old_transform_h1_readiness_probe.jl`
- `tmp/work/pqs_source_box_old_transform_h1_readiness_probe_output.txt`

Probe shape:

- Built a 13-point mapped basis using the pattern from
  `test/nested/cartesian_nested_face_first_primitive_runtests.jl`.
- Used the old private nested helper only to produce source-axis coefficient
  matrices:
  `GaussletBases._cartesian_raw_product_box_plan(...)`.
- Extracted 11 x 5 per-axis local coefficient matrices over source interval
  `2:12`.
- Fed those matrices into a new
  `CartesianRawProductSources.raw_product_box_plan(...)` as
  `axis_transform_matrices`.
- Built the PQS boundary retained rule and confirmed retained count `98`.
- Built retained-source `S`, `T`, and uncharged by-center `V_unit` through the
  current CPBM helpers.
- Formed `H = T + Z * V_unit` only inside the probe.
- No production code, tests, shell realization, Lowdin cleanup, IDA,
  density-density, RHF, global driver, exports, artifacts, full-parent
  fallback, direct Cartesian fallback, fixed-block matrix authority, or CCPM
  wrapper path was added.

Key diagnostics:

- transform source: `:old_nested_source_axis_transform_oracle`
- production/physical PQS H1 blocker:
  `:missing_repo_owned_pqs_source_axis_transform_builder`
- basis count: `13`
- source interval: `2:12`
- old transform source mode dims: `(5, 5, 5)`
- old transform max axis overlap error: `5.266241002305216e-16`
- old transform source product modes orthogonal: `true`
- transform shapes: `((11, 5), (11, 5), (11, 5))`
- transform nonidentity: `true`
- raw source mode dims/count: `(5, 5, 5)` / `125`
- retained boundary count/dimension: `98`
- record transform statuses/shapes:
  `(:materialized, :materialized, :materialized)` /
  `((11, 5), (11, 5), (11, 5))`
- `S`, `T`, `V_unit`, and `H` shapes: `(98, 98)`
- finite entries: true for all four matrices
- symmetry error: `0.0` for all four matrices
- overlap eigenvalue range:
  `0.9999999999999973` to `1.0000000000000033`
- overlap rank: `98`
- overlap condition estimate: `1.000000000000006`
- overlap identity error: `6.661338147750939e-16`
- solve kind: `ordinary_retained_source_diagnostic`
- lowest diagnostic H1 energy: `0.0320561000473788`
- nuclear charge recorded: `true`
- nuclear charge applied in source block: `false`
- centers summed: `false`
- shell realization materialized: `false`
- Lowdin cleanup used: `false`
- IDA data materialized: `false`
- Hamiltonian data materialized in source route: `false`

Readiness result and blocker:

- The old nested helper proves the non-identity source-axis coefficient shape
  can feed the current `CartesianRawProductSources` and CPBM retained-source
  path.
- The old helper remains oracle/kernel reference only; it was not made route
  authority.
- Exact production blocker:
  `:missing_repo_owned_pqs_source_axis_transform_builder`.
- The ordinary H1 value is still a retained-source diagnostic only, before
  production-owned PQS transform construction, shell realization, Lowdin/final
  basis construction, and global/driver adoption.

Validation:

- `julia --project=. tmp/work/pqs_source_box_old_transform_h1_readiness_probe.jl`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git diff --check`

Deletion/shrinkage report:

- No old code, test, metadata, or compatibility path became unnecessary in this
  probe-only pass.
- Nothing was deleted or simplified because the pass intentionally used ignored
  `tmp/work` artifacts and kept the old nested helper as oracle only.
- No new long-term test was added; the probe artifact earns its cost by
  proving the non-identity transform shape and sharpening the blocker from
  generic missing transforms to a missing repo-owned transform builder.
- Remaining stale/duplicate surfaces to retire next: after a repo-owned PQS
  source-axis transform builder exists, the old nested source-box transform
  helpers should be quarantined as oracle/reference or deleted where no longer
  called.

-- repo-doer@macmini
