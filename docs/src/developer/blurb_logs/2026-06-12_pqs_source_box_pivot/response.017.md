Ran the 5 x 5 x 5 PQS retained-source H1 readiness probe with the repo-owned source-axis transform facts.

Artifacts:

- `tmp/work/pqs_source_box_repo_transform_h1_readiness_probe.jl`
- `tmp/work/pqs_source_box_repo_transform_h1_readiness_probe_output.txt`

What the probe did:

- Built the same 13-point mapped-basis / interval `2:12` / source dims `(5, 5, 5)` fixture used in pass 015.
- Called `CartesianPairBlockMaterialization.pqs_source_axis_transform_facts_from_pgdg_axes(...)`.
- Passed the returned `axis_transform_facts` into `CartesianRawProductSources.raw_product_box_plan(...)`.
- Built the PQS boundary retained rule and retained-source `S`, `T`, and uncharged by-center `V_unit` through current CPBM helpers.
- Formed `H = T + Z * V_unit` only inside the probe.
- Did not call the old source-box helper, old fixed-block matrix authority, CCPM wrappers, shell realization, Lowdin cleanup, IDA, density-density, RHF, global driver, exports, artifacts, full-parent CPB, or direct Cartesian fallback.

Key diagnostics:

- `transform_source = :repo_owned_pgdg_doside_source_axis_transform`
- `old_source_box_helper_called = false`
- `old_fixed_block_matrix_authority_used = false`
- `basis_count = 13`
- `source_interval = 2:12`
- `transform_source_mode_dims = (5, 5, 5)`
- `transform_shapes = ((11, 5), (11, 5), (11, 5))`
- `transform_max_axis_overlap_error = 5.12225151339942e-16`
- `raw_source_mode_count = 125`
- `retained_boundary_count = 98`
- `retained_dimension = 98`
- `overlap_rank = 98`
- `overlap_condition_estimate = 1.000000000000006`
- `overlap_identity_error = 6.661338147750939e-16`
- `solve_kind = :ordinary_retained_source_diagnostic`
- `lowest_h1_energy = 0.0320561000473788`
- `h1_delta_from_pass015 = 0.0`
- `overlap_identity_error_delta_from_pass015 = 0.0`

Readiness result:

The repo-owned PGDG doside transform builder reproduces the pass-015 old-transform oracle H1 and overlap diagnostics to roundoff. The source-axis transform blocker is cleared for this probe.

The exact next blocker for physical/final PQS H1 is:

```text
:missing_pqs_shell_realization_lowdin_final_basis_construction
```

The ordinary solve remains a retained-source diagnostic only. It is not a final physical PQS acceptance solve because shell realization, Lowdin/final-basis construction, and the source-box-to-final PQS basis boundary are not yet materialized.

Validation:

- `julia --project=. tmp/work/pqs_source_box_repo_transform_h1_readiness_probe.jl`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git diff --check`

Deletion/shrinkage report:

- No tracked production code or tests changed in this pass; this was intentionally probe-only.
- Nothing was deleted yet because the old helper is still useful as historical/oracle context until the downstream physical PQS final-basis path consumes the repo-owned transform builder.
- The pass-017 probe supersedes the pass-015 old-helper probe as the active H1 readiness artifact for source-axis transforms.
- No new long-term test was added. The ignored `tmp/work` artifact earned its temporary carrying cost by proving the H1 probe no longer depends on the old nested source-box helper.
- Remaining stale/duplicate surfaces to retire next: old nested source-box readiness probes and helper vocabulary should be quarantined or deleted once the shell-realization/Lowdin final-basis boundary is implemented and covered by a live PQS acceptance path.

-- repo-doer@macmini
