Pass 071 response: PQS multilayer support kinetic helper

Summary:
- Added `pqs_multilayer_support_kinetic_matrix(plan)` in `src/pqs_multilayer_shell_source_plan.jl`.
- The helper requires a `:pqs_multilayer_shell_source_plan` with status `:available_pqs_multilayer_shell_source_plan`.
- It builds the support-space kinetic matrix over `vcat(plan.core_support_states, plan.shell_support_states)` using the standard three-term Cartesian product form:
  - `Kx * Sy * Sz`
  - `Sx * Ky * Sz`
  - `Sx * Sy * Kz`
- It returns the support-space matrix directly and makes no final-basis transfer, H1, nuclear, IDA, RHF, driver, export, or artifact claim.
- Updated the tracked H1 gate to consume `GaussletBases.pqs_multilayer_support_kinetic_matrix(fixture.plan)`.
- Deleted the local `_pqs_h1_support_kinetic_matrix` helper.
- Updated `docs/src/developer/pqs_near_term_final_basis_realization_plan.md` to mark support kinetic implemented while leaving the electron-nuclear wait note intact.

Validation:
- `julia --project=. -e 't = @elapsed include("test/nested/pqs_direct_retained_final_h1_runtests.jl"); println("elapsed_s=", t)'`
  - PASS: 31 tests, elapsed_s=62.293904292
  - Note: this included package precompilation.
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - PASS: `load ok`
- `git diff --check`
  - PASS

Diff:
- `src/pqs_multilayer_shell_source_plan.jl`
- `test/nested/pqs_direct_retained_final_h1_runtests.jl`
- `docs/src/developer/pqs_near_term_final_basis_realization_plan.md`
- Total: 55 insertions, 40 deletions

Deletion/shrinkage report:
- Became unnecessary:
  - the test-local `_pqs_h1_support_kinetic_matrix` helper.
- Deleted/simplified:
  - removed `_pqs_h1_support_kinetic_matrix`;
  - replaced its use in the H1 gate with the route-owned multi-layer support kinetic helper.
- No new test file was added. The existing H1 gate now covers the live helper through the current endpoint workflow.
- Remaining stale/duplicate surface:
  - `_pqs_h1_support_nuclear_matrix` remains test-local because support electron-nuclear still needs the separated uncharged by-center convention, sign convention, charge-application stage, and centered/off-origin Gaussian factor source contract before promotion.

-- repo-doer@macmini
