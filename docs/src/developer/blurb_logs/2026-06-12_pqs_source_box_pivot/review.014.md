Accepted as a raw-source vocabulary readiness probe.

Pass 014 moved from the synthetic 3 x 3 x 3 fixture to the repo-owned
`CartesianRawProductSources` vocabulary:

```text
RawProductBoxPlan(source_mode_dims = (5, 5, 5))
PQS boundary retained rule -> retained count 98
retained-source S/T/V_unit blocks
```

The probe correctly labeled its materialized identity axis transforms as
`:diagnostic_identity_source_mode_transform`. The matrix wiring works, but the
ordinary H1 diagnostic is not physical PQS because the source-axis transforms
are identity matrices supplied by the probe.

Manager validation:

- `julia --project=. tmp/work/pqs_source_box_5x5_h1_readiness_probe.jl`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git diff --check`

Current blocker:

```text
:missing_pqs_source_axis_transform_builder
```

Next target:

Use the existing old nested source-box transform helper only as an
oracle/source of non-identity axis coefficient matrices, then feed those
matrices through the new `CartesianRawProductSources` plan and CPBM retained
source path. This should stay a `tmp/work` probe. It should not adopt the old
helper as route authority.

Known old surfaces:

- `src/cartesian_nested_faces.jl`
- `_cartesian_source_box_axis_transform(...)`
- `_cartesian_source_box_axis_transform_plan(...)`
- `_cartesian_raw_product_box_plan(...)`
- `test/nested/cartesian_nested_face_first_primitive_runtests.jl`

-- repo-manager@macmini
