Pass 086 complete: extracted the multi-layer PQS support one-body helpers into a same-module file.

Files changed:

- `src/GaussletBases.jl`
- `src/pqs_multilayer_shell_source_plan.jl`
- `src/pqs_multilayer_support_one_body.jl`

What moved:

- `_pqs_multilayer_support_product_matrix`
- `_pqs_multilayer_axis_tuple`
- `_pqs_multilayer_center_records`
- `_pqs_multilayer_center_property`
- `_pqs_multilayer_center_summary`
- `_pqs_multilayer_explicit_factor_terms`
- `_pqs_multilayer_centered_factor_terms`
- `_pqs_multilayer_term_first_factor_array`
- `_pqs_multilayer_validate_factor_terms`
- `_pqs_multilayer_support_electron_nuclear_matrix`
- `pqs_multilayer_support_kinetic_matrix`
- `pqs_multilayer_support_electron_nuclear_by_center_matrices`

What stayed:

- `src/pqs_multilayer_shell_source_plan.jl` still owns explicit-box bridge specs, region-plan layer specs, common shell source realization, public source-plan entry points, complete core/shell final-basis adapter, H1 payload, and the blocked final-basis result helper.

Include order used:

```julia
include("pqs_multilayer_shell_region_plan.jl")
include("pqs_multilayer_support_one_body.jl")
include("pqs_multilayer_shell_source_plan.jl")
```

This order keeps shellification/lowering region utilities first, then support one-body helpers, then source realization and payload code. No circular include pressure appeared because the support helpers consume the source plan structurally and only require `_pqs_multilayer_property(...)` plus shared source-plan field names.

Effect:

- `src/pqs_multilayer_shell_source_plan.jl` shrank from 873 lines to 527 lines.
- `src/pqs_multilayer_support_one_body.jl` is 347 lines.
- The dense support-space one-body machinery is now visibly scoped to H1 seam/reference machinery rather than mixed into source-plan realization.

Validation:

- `julia --project=. test/nested/pqs_direct_retained_final_h1_runtests.jl`
  - passed: 46 tests.
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - passed.
- `git diff --check`
  - passed.
- Whitespace scan on touched/new source files found no trailing whitespace.

Deletion/shrinkage report:

- Conceptual responsibility separated: dense support-space product/kinetic/electron-nuclear one-body helpers moved out of the source-plan realization file.
- No obsolete production behavior was deleted; this was a mechanical file-boundary split preserving function names and return shapes.
- No tests were added; the existing H1 gate covers the live behavior touched by the include/file split.
- This makes the source-plan file less of a private route island: it now starts with explicit-box/region source realization rather than dense operator helper machinery.
- Remaining stale or duplicate surfaces to retire next:
  - complete core/shell final-basis adapter and H1 payload still live in `pqs_multilayer_shell_source_plan.jl`;
  - explicit-box bridge remains compatibility/probe machinery in the source-plan file;
  - support one-body helpers remain dense H1-seam helpers, not general PQS operator assembly.

-- repo-doer@macmini
