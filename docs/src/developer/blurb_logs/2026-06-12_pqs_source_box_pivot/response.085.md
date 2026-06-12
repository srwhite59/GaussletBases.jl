Pass 085 complete: did the PQS multi-layer file-boundary audit and performed one mechanical same-module extraction.

Files changed:

- `src/GaussletBases.jl`
- `src/pqs_multilayer_shell_source_plan.jl`
- `src/pqs_multilayer_shell_region_plan.jl`

Recommended file split:

1. `pqs_multilayer_shell_region_plan.jl`
   - shellification/lowering-backed region-plan records and constructor;
   - compact coverage/duplicate support checks;
   - shared plan utility helpers used by source realization.
2. `pqs_multilayer_shell_source_plan.jl`
   - explicit-box bridge layer specs;
   - region-plan layer specs;
   - common repeated one-cell PQS shell source realization;
   - public `pqs_multilayer_shell_source_plan(...)` entry points.
3. Future split: `pqs_multilayer_support_one_body.jl`
   - dense support product helper;
   - support kinetic helper;
   - support electron-nuclear by-center helper;
   - center/factor-term helpers.
4. Future split: `pqs_multilayer_complete_core_shell_final_basis_payloads.jl` or two smaller files:
   - complete core/shell final-basis adapter;
   - H1 payload.

Recommended include order:

```julia
include("pqs_multilayer_shell_region_plan.jl")
include("pqs_multilayer_shell_source_plan.jl")
```

If the future support/H1 files are extracted, the dependency order should be:

```julia
include("pqs_multilayer_shell_region_plan.jl")
include("pqs_multilayer_shell_source_plan.jl")
include("pqs_multilayer_support_one_body.jl")
include("pqs_multilayer_h1_payload.jl")
```

The support one-body file depends on source-plan shape and shared support ordering. The H1 payload depends on support one-body plus `CartesianFinalBasisRealization`.

Mechanical extraction implemented:

- Moved the region-plan section into `src/pqs_multilayer_shell_region_plan.jl`.
- Added that file before `src/pqs_multilayer_shell_source_plan.jl` in `src/GaussletBases.jl`.
- Preserved function names and return shapes.
- Did not move code across modules.
- Did not add tests or new route behavior.

Effect:

- `src/pqs_multilayer_shell_source_plan.jl` shrank from 1111 lines to 873 lines.
- The new region-plan file is 239 lines.
- The active source-plan file now starts at source/support realization rather than shellification/lowering region ownership.

Validation:

- `julia --project=. test/nested/pqs_direct_retained_final_h1_runtests.jl`
  - passed: 46 tests.
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - passed.
- `git diff --check`
  - passed.
- Whitespace scan on the touched/new source files found no trailing whitespace.

Deletion/shrinkage report:

- Conceptual responsibility separated: shellification/lowering-backed region planning is now in `pqs_multilayer_shell_region_plan.jl`, away from PQS source realization and H1 payload machinery.
- The original source-plan file still owns explicit-box bridge specs, source realization, dense support-space product/kinetic/nuclear helpers, final-basis adapter, and H1 payload.
- No obsolete production path was deleted; this was a same-module mechanical split to reduce route-island pressure without behavior changes.
- No test was added; existing H1 gate covers the live behavior touched by the include split.
- Remaining private route-island pressure:
  - support one-body helpers remain in the source-plan file but should move next if a follow-up split is wanted;
  - final-basis adapter and H1 payload remain in the same file and can be separated after support one-body is split;
  - explicit-box bridge remains in the source-plan file as compatibility/probe machinery.

-- repo-doer@macmini
