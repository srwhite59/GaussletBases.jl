Accepted as an oracle-backed transform-shape probe.

Pass 015 used the old nested source-box helper only to produce non-identity
axis coefficient matrices, then fed those matrices into the new
`CartesianRawProductSources` plan and current CPBM retained-source one-body
path.

Important result:

```text
old helper axis local coefficient shapes: (11, 5), (11, 5), (11, 5)
retained boundary count: 98
overlap identity error: 6.66e-16
lowest diagnostic H1: 0.0320561000473788
```

The H1 value is still only a retained-source diagnostic. The old helper is not
route authority. The useful conclusion is that non-identity source-axis
coefficient matrices with this shape can feed the new raw-source/CPBM path.

Manager validation:

- `julia --project=. tmp/work/pqs_source_box_old_transform_h1_readiness_probe.jl`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git diff --check`

Current blocker:

```text
:missing_repo_owned_pqs_source_axis_transform_builder
```

Next target:

Implement a narrow repo-owned source-axis transform fact builder, probably in
CPBM rather than `CartesianRawProductSources`, because CRPS is intentionally
metadata-only and depends only on CPB. The new helper may reuse the low-level
PGDG/doside transform kernel, but it should not call the old
`_cartesian_raw_product_box_plan(...)` wrapper or old fixed-block matrices.

-- repo-manager@macmini
