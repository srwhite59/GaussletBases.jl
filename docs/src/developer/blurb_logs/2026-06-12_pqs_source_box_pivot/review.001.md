Accepted.

The doer stayed inside the read-only audit boundary and identified the correct
module-spine starting point for the PQS pivot.

Key conclusions:

- `CartesianRawProductSources` owns raw product source-box facts through
  `RawProductBoxPlan`.
- `CartesianPairBlockMaterialization` already owns PQS/PQS raw source-space
  safe one-body blocks in `pqs_source_safe_terms.jl`.
- Existing shell-realization/readiness layers are metadata bridge/readiness
  surfaces and should remain oracle/adapter/debug only for now.
- Older CCPM and fixed-block helpers are reference/migration surfaces, not new
  route authority.
- The first implementation target should be a source-mode retained boundary
  rule/transform, not shell projection or Lowdin.

The first retained rule target is concrete and reviewable:

```text
source modes: 5 x 5 x 5 = 125
interior modes: 3 x 3 x 3 = 27
boundary retained modes: 98
```

This is source-mode boundary selection, not support-row or shell-row boundary
selection.

Validation reported by doer:

- `git status --short --branch`
- `git diff --check`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`

Next target:

Implement a compact PQS boundary source-mode retained rule under
`CartesianRawProductSources`, then attach its summary to PQS retained-unit
transform-contract metadata. Do not implement pair-block contraction, shell
realization, Lowdin, electron-nuclear, IDA, or drivers in the next pass.

-- repo-manager@macmini
