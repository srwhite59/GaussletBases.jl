Purpose:

Rerun the 5 x 5 x 5 PQS retained-source H1 readiness probe with real
non-identity source-axis transforms sourced from the old nested helper as an
oracle/kernel reference.

Why now:

Pass 014 proved the `CartesianRawProductSources` raw-plan and 98-mode boundary
rule wiring, but only with identity source-axis transforms. The repo already
has older private source-box transform helpers:

- `src/cartesian_nested_faces.jl`
- `_cartesian_source_box_axis_transform(...)`
- `_cartesian_source_box_axis_transform_plan(...)`
- `_cartesian_raw_product_box_plan(...)`
- `test/nested/cartesian_nested_face_first_primitive_runtests.jl`

Use them only to supply non-identity axis coefficient matrices for a probe.
Do not make the old helper the new route authority.

Exact task:

Create a `tmp/work` probe. Do not add production code or tests in this pass
unless the probe exposes a very small missing seam.

The probe should:

- build an old nested raw product-box plan with source dims `(5, 5, 5)` using
  the pattern in `cartesian_nested_face_first_primitive_runtests.jl`;
- extract its per-axis local coefficient matrices;
- build a new `CartesianRawProductSources.raw_product_box_plan(...)` using
  those matrices as `axis_transform_matrices`;
- build the PQS boundary retained rule and confirm retained count `98`;
- build retained-source `S`, `T`, and uncharged by-center `V_unit` through the
  current CPBM helpers;
- form `H = T + Z * V_unit` only in the probe;
- report overlap rank/condition/identity error, solve kind, lowest diagnostic
  H1 value if solved, and whether the transforms are non-identity;
- label the transform source as `:old_nested_source_axis_transform_oracle`.

Expected outcome:

If the probe works, the next production blocker should become:

```text
:missing_repo_owned_pqs_source_axis_transform_builder
```

That means the old helper proved the numerical transform shape, but the new
route still needs a module-owned builder before physical/final PQS H1.

Trust boundary:

No shell realization, Lowdin cleanup, IDA, density-density, RHF, global driver,
exports, artifacts, full-parent fallback, direct Cartesian fallback, or old
CCPM wrapper adoption. Do not call old fixed-block matrices as authority.

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

- write `.agent_handoffs/response.015.md.tmp`, then atomically rename to
  `.agent_handoffs/response.015.md`;
- also write the curated copy to
  `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.015.md`;
- include probe artifact path and key diagnostics;
- include exact blocker for production/physical PQS H1;
- include validation run;
- include deletion/shrinkage report;
- sign `-- repo-doer@macmini`.

-- repo-manager@macmini
