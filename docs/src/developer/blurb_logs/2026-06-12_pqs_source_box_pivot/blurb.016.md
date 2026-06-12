Purpose:

Implement the narrow repo-owned PQS source-axis transform fact builder exposed
by the pass-015 probe.

Why now:

The new raw-source and CPBM retained-source one-body path works with
non-identity source-axis transforms, but pass 015 still borrowed those
transforms from the old private nested wrapper:

```text
_cartesian_raw_product_box_plan(...)
```

The next blocker is:

```text
:missing_repo_owned_pqs_source_axis_transform_builder
```

Exact task:

Add a narrow helper in CPBM or another appropriate non-CRPS seam. Do not put
PGDG/doside numerical construction into `CartesianRawProductSources`, because
CRPS is intentionally metadata-only and depends only on `CartesianCPB`.

Recommended API shape:

```text
pqs_source_axis_transform_facts_from_pgdg_axes(
    axis_sources;
    source_intervals,
    source_mode_dims,
    enforce_symmetric_odd = false,
)
```

or a better local name if the codebase suggests one.

Required behavior:

- accept x/y/z axis PGDG intermediates or mapped ordinary gausslet bundles;
- call the low-level doside/PGDG source transform kernel directly, e.g.
  `_nested_doside_1d(...)`;
- produce `CartesianRawProductSources.AxisSourceTransformFact`s with
  materialized coefficient matrices shaped
  `(length(source_interval_axis), source_mode_dim_axis)`;
- record compact diagnostics: requested/resolved source dims, max overlap
  error, source mode adjustment flag, transform source, no shell/Lowdin/IDA;
- not call `_cartesian_raw_product_box_plan(...)` or old fixed-block matrices.

If a complete raw plan helper is clearly smaller and cleaner, it may compose
these facts with `CartesianRawProductSources.raw_product_box_plan(...)`, but
keep the source-axis transform builder visible and testable.

Test policy:

Add one compact test in the existing CPBM contract file or another existing
focused contract file. Compare the new helper against the old nested helper as
an oracle for the 5 x 5 x 5 fixture used in pass 015:

- source intervals `2:12`;
- source dims `(5, 5, 5)`;
- per-axis coefficient matrices match old `axis_local_coefficients`;
- max overlap error is small;
- no shell realization, Lowdin, IDA, Hamiltonian, driver, or artifact claim.

Do not add a new test file unless there is a strong reason.

Trust boundary:

No H1 solve, shell realization, Lowdin cleanup, IDA, density-density, RHF,
global driver, exports, artifacts, full-parent fallback, direct Cartesian
fallback, CCPM wrapper adoption, or old fixed-block matrix authority.

Deletion/shrinkage report required:

- what old code, test, metadata, or compatibility path became unnecessary;
- what was deleted or simplified;
- if nothing was deleted, why no existing surface was made obsolete yet;
- whether any new test/artifact was added and why it earned its cost;
- any remaining stale or duplicate surfaces to retire next.

Validation:

- `julia --project=. test/nested/cartesian_pair_block_materialization_contract_runtests.jl`
  or the focused existing test file touched by the change;
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git diff --check`

Report back:

- write `.agent_handoffs/response.016.md.tmp`, then atomically rename to
  `.agent_handoffs/response.016.md`;
- also write the curated copy to
  `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.016.md`;
- include implementation or exact blocker;
- include validation run;
- include deletion/shrinkage report;
- sign `-- repo-doer@macmini`.

-- repo-manager@macmini
