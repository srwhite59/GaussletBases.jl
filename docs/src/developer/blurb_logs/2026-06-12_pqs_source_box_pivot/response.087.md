Pass 087 complete: extracted the complete core/shell final-basis and H1 helpers into a same-module file.

Files changed:

- `src/GaussletBases.jl`
- `src/pqs_multilayer_shell_source_plan.jl`
- `src/pqs_multilayer_complete_core_shell_h1.jl`

What moved:

- `pqs_multilayer_complete_core_shell_final_basis`
- `pqs_multilayer_complete_core_shell_h1_payload`
- `_blocked_pqs_multilayer_complete_core_shell_final_basis`

What stayed:

- `src/pqs_multilayer_shell_source_plan.jl` now owns only:
  - explicit-box bridge layer specs;
  - region-plan layer specs;
  - common repeated one-cell PQS shell source realization;
  - public source-plan entry points.

Include order used:

```julia
include("pqs_multilayer_shell_region_plan.jl")
include("pqs_multilayer_support_one_body.jl")
include("pqs_multilayer_shell_source_plan.jl")
include("pqs_multilayer_complete_core_shell_h1.jl")
```

This matches the requested dependency order: region planning first, support one-body helpers second, source realization third, and final-basis/H1 seam last.

Effect:

- `src/pqs_multilayer_shell_source_plan.jl` shrank from 527 lines to 265 lines.
- `src/pqs_multilayer_complete_core_shell_h1.jl` is 263 lines.
- The source-plan file now mostly owns source realization plus explicit-box bridge machinery.

Validation:

- `julia --project=. test/nested/pqs_direct_retained_final_h1_runtests.jl`
  - passed: 46 tests.
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - passed.
- `git diff --check`
  - passed.
- Whitespace scan on touched/new source files found no trailing whitespace.

Deletion/shrinkage report:

- Conceptual responsibility separated: complete core/shell final-basis adaptation and H1 payload assembly are now outside the source realization file.
- No old behavior was deleted; this was a same-module mechanical split preserving function names and return shapes.
- No tests were added; the focused H1 gate remains the live contract for this split.
- `pqs_multilayer_shell_source_plan.jl` is now primarily source realization plus explicit-box bridge logic, so it is substantially less of a private route island.
- Later module-boundary decision remains: the final-basis/H1 seam may eventually belong closer to `CartesianFinalBasisRealization`, but that is a conceptual module move and was intentionally not done in this file-boundary pass.
- Remaining stale or duplicate surfaces to retire next:
  - explicit-box bridge remains compatibility/probe machinery;
  - dense support one-body helpers are still H1-seam/reference helpers, not general PQS operator assembly;
  - future H1/J work should avoid growing these files back into a private density/RHF route island.

-- repo-doer@macmini
