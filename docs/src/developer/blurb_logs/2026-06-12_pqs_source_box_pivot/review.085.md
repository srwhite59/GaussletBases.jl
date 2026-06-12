Review 085:

Accepted. This was the right kind of correction pass: a same-module mechanical
file split that separates shellification/lowering-backed PQS region planning
from the source-realization/H1 seam without changing route behavior.

What I checked:

- `src/pqs_multilayer_shell_region_plan.jl` now owns the region-plan records,
  coverage checks, and `pqs_multilayer_shell_region_plan(...)` constructor.
- `src/pqs_multilayer_shell_source_plan.jl` no longer starts by defining the
  shellification/lowering region-plan concept.
- `src/GaussletBases.jl` includes the new region-plan file before the source
  realization file.
- Function names and return shapes were preserved.
- No new tests, physics behavior, driver wiring, IDA, density, RHF, exports, or
  artifacts were added.

Validation rerun by manager:

- `julia --project=. test/nested/pqs_direct_retained_final_h1_runtests.jl`
  passed, 46 tests.
- `julia --project=. -e 'using GaussletBases; println("load ok")'` passed.
- `git diff --check` passed.

Deletion / shrinkage:

- Conceptual responsibility separated: shellification/lowering-backed region
  planning is no longer embedded in the source-plan/H1 file.
- `src/pqs_multilayer_shell_source_plan.jl` shrank by moving the region-plan
  section out.
- No obsolete production path was deleted because this was deliberately a
  behavior-preserving file-boundary pass.
- No test was added; the existing H1 gate remained the behavior check.

Remaining concern:

- `src/pqs_multilayer_shell_source_plan.jl` still mixes source realization,
  dense support one-body helpers, final-basis adapter, and H1 payload. The next
  cleanup should extract the support one-body helpers as another same-module
  mechanical split before adding more PQS physics.

-- repo-manager@macmini
