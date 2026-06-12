Purpose:

Rerun the 5 x 5 x 5 PQS retained-source H1 readiness probe using the new
repo-owned source-axis transform fact builder.

Why now:

Pass 016 implemented:

```text
pqs_source_axis_transform_facts_from_pgdg_axes(...)
```

The next check is to remove the pass-015 probe's direct dependency on the old
nested raw-product-box helper. The probe should now use the new builder to
create materialized CRPS axis transform facts, then feed those facts into the
same raw-source and CPBM retained-source one-body path.

Exact task:

Create or update a `tmp/work` probe. Do not add a new test.

The probe should:

- build the same 13-point mapped-basis / interval `2:12` / source dims
  `(5, 5, 5)` fixture used in pass 015;
- call `pqs_source_axis_transform_facts_from_pgdg_axes(...)`;
- pass the returned `axis_transform_facts` into
  `CartesianRawProductSources.raw_product_box_plan(...)`;
- build the PQS boundary retained rule and confirm retained count `98`;
- build retained-source `S`, `T`, and uncharged by-center `V_unit` through the
  current CPBM helpers;
- form `H = T + Z * V_unit` only inside the probe;
- report overlap rank/condition/identity error, solve kind, lowest diagnostic
  H1 value if solved, and whether any old source-box helper was called.

Expected result:

The H1 diagnostics should match the pass-015 old-transform oracle probe to
roundoff. The blocker should advance from
`:missing_repo_owned_pqs_source_axis_transform_builder` to the next real
physical/final-basis blocker, likely shell realization / Lowdin final-basis
construction or a clearer source-box-to-final PQS basis boundary.

Trust boundary:

No production code unless the probe exposes a tiny wiring issue. No shell
realization, Lowdin cleanup, IDA, density-density, RHF, global driver, exports,
artifacts, full-parent fallback, direct Cartesian fallback, CCPM wrapper
adoption, or old fixed-block matrix authority.

Test policy:

No new test in this pass. Use `tmp/work` artifacts and report paths.

Deletion/shrinkage report required:

- what old code, test, metadata, or compatibility path became unnecessary;
- what was deleted or simplified;
- if nothing was deleted, why no existing surface was made obsolete yet;
- whether any new test/artifact was added and why it earned its cost;
- any remaining stale or duplicate surfaces to retire next.

Validation:

- `julia --project=. <tmp/work probe>`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git diff --check`

Report back:

- write `.agent_handoffs/response.017.md.tmp`, then atomically rename to
  `.agent_handoffs/response.017.md`;
- also write the curated copy to
  `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.017.md`;
- include probe artifact path and key diagnostics;
- include exact next blocker for physical/final PQS H1;
- include validation run;
- include deletion/shrinkage report;
- sign `-- repo-doer@macmini`.

-- repo-manager@macmini
