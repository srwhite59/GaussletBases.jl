Purpose:

Run the one-center PQS retained-source H1 readiness probe on the real
raw-product source-box vocabulary, using the 5 x 5 x 5 boundary retained rule.

Why now:

Pass 013 proved the retained-source one-body pieces compose, but only with a
synthetic 3 x 3 x 3 identity-transform fixture. The next probe should use the
actual raw source-box objects and retained-rule convention already in the repo:

- `CartesianRawProductSources.raw_product_box_plan(...)`
- `CartesianRawProductSources.pqs_boundary_product_mode_retained_rule(...)`
- `source_mode_dims = (5, 5, 5)`
- retained boundary mode count = `98`

Exact task:

Create a `tmp/work` probe. Do not add a long-term test.

Use these source surfaces as starting points:

- `src/cartesian_raw_product_sources/records.jl`
- `src/cartesian_raw_product_sources/axis_transform_facts.jl`
- `test/nested/cartesian_raw_product_sources_contract_runtests.jl`
- `test/nested/cartesian_retained_unit_transform_contracts_runtests.jl`

The probe should:

- build a `RawProductBoxPlan` with `source_mode_dims = (5, 5, 5)`;
- build the PQS boundary retained rule and confirm retained count `98`;
- attach materialized source-axis transform facts if possible;
- if the only available transform is identity, label the run
  `:diagnostic_identity_source_mode_transform`, not physical PQS;
- build retained-source `S`, `T`, and uncharged by-center `V_unit` using the
  existing CPBM helpers from pass 013;
- form `H = T + Z * V_unit` only inside the probe;
- report matrix shapes, symmetry, overlap rank/condition/identity error, solve
  kind if solved, and exact blocker for physical/final PQS H1.

Do not hide the distinction between identity source modes and real PQS
source-axis transforms. If a real source-axis transform builder is missing,
report `:missing_pqs_source_axis_transform_builder` or a more precise blocker.

Trust boundary:

No production route unless the probe exposes a tiny missing seam whose fix is
clearly less than the probe scaffolding. No shell realization, Lowdin cleanup,
IDA, density-density, RHF, global driver, exports, artifacts, full-parent
fallback, direct Cartesian fallback, or old CCPM wrapper adoption.

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

- write `.agent_handoffs/response.014.md.tmp`, then atomically rename to
  `.agent_handoffs/response.014.md`;
- also write the curated copy to
  `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.014.md`;
- include probe artifact path and key diagnostics;
- include exact blocker for physical/final PQS H1;
- include validation run;
- include deletion/shrinkage report;
- sign `-- repo-doer@macmini`.

-- repo-manager@macmini
