Ran the probe-first one-center PQS retained-source H1 readiness check.

Artifacts:

- `tmp/work/pqs_retained_source_h1_readiness_probe.jl`
- `tmp/work/pqs_retained_source_h1_readiness_probe_output.txt`

Probe shape:

- Built one synthetic PQS retained self-pair with materialized identity
  source-axis transforms over a 3 x 3 x 3 source box.
- Used explicit `UniformBasisSpec` axis layers.
- Built retained-source `S` with
  `pqs_source_pair_retained_overlap_block(...)`.
- Built retained-source `T` with
  `pqs_source_pair_retained_kinetic_block(...)`.
- Built one uncharged by-center retained-source `V_unit` with
  `pqs_source_pair_retained_centered_electron_nuclear_by_center_block(...)`.
- Formed `H = T + Z * V_unit` only inside the probe, as the Hamiltonian-stage
  charge-application diagnostic.
- No production code, tests, route/global driver, shell realization, Lowdin,
  IDA, density-density, RHF, exports, artifacts, full-parent fallback, direct
  Cartesian fallback, or CCPM wrapper path was added.

Key diagnostics:

- retained dimension: `26`
- source mode dims/count: `(3, 3, 3)` / `27`
- `S`, `T`, `V_unit`, and `H` shapes: `(26, 26)`
- finite entries: true for all four matrices
- symmetry error: `0.0` for all four matrices
- overlap eigenvalue range:
  `0.999999999999999` to `1.0000000000000024`
- overlap rank: `26`
- overlap condition estimate: `1.0000000000000036`
- overlap identity error: `1.3322676295501878e-15`
- solve kind: `ordinary_retained_source_diagnostic`
- lowest diagnostic H1 energy: `0.9464341674675832`
- nuclear charge recorded: `true`
- nuclear charge applied in source block: `false`
- centers summed: `false`
- shell realization materialized: `false`
- Lowdin cleanup used: `false`
- IDA data materialized: `false`
- Hamiltonian data materialized in source route: `false`

Readiness result:

- The retained-source one-center pieces compose coherently for this synthetic
  self-pair probe.
- The retained-source overlap is identity to roundoff in this fixture, so the
  probe solve is an ordinary retained-source diagnostic, not a raw generalized
  diagnostic.
- This is still not a final PQS acceptance path. It remains before shell
  realization and Lowdin/final-basis construction, and the reported H1 number
  is only a source-box wiring diagnostic.

Validation:

- `julia --project=. tmp/work/pqs_retained_source_h1_readiness_probe.jl`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git diff --check`

Deletion/shrinkage report:

- No old code, test, metadata, or compatibility path became unnecessary in this
  probe-only pass.
- Nothing was deleted or simplified because the pass intentionally avoided
  production code changes and used ignored `tmp/work` artifacts.
- No new long-term test was added; the probe artifact earns its cost by
  recording the first retained-source H1 wiring/readiness check without adding
  suite carrying cost.
- Remaining stale/duplicate surfaces to retire next: old CCPM source-box H1 or
  nuclear vocabulary should stay quarantined until a downstream PQS shell/
  Lowdin route consumes this retained-source block path or proves the old
  helper is no longer needed as an oracle.

-- repo-doer@macmini
