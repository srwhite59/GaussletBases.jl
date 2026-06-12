Pass 073 response: PQS support electron-nuclear by-center helper

Summary:
- Added `pqs_multilayer_support_electron_nuclear_by_center_matrices(plan; ...)` in `src/pqs_multilayer_shell_source_plan.jl`.
- The helper requires an available `:pqs_multilayer_shell_source_plan`.
- It returns one support-space matrix record per center, ordered over core support rows followed by shell support rows.
- Each center matrix follows the documented convention:
  - negative unit-charge electron-nuclear attraction, `-1/r_A`;
  - nuclear charge recorded but not applied;
  - centers not summed;
  - no final-basis transfer, H1, IDA, density-density, RHF, driver wiring, exports, or artifacts.
- Implemented the explicit factor path used by the current H1 gate:
  - one origin `pgdg_intermediate.gaussian_factor_terms` array can be reused for x/y/z for a single center.
- Implemented the axis-layer centered factor path:
  - when explicit factors are not supplied, the helper uses `axis_layers`, Coulomb expansion exponents, and center locations through the low-level `gaussian_factor_matrices(...; center = ...)` convention.
- Updated `test/nested/pqs_direct_retained_final_h1_runtests.jl` to use the route-owned helper and deleted `_pqs_h1_support_nuclear_matrix`.
- Updated `docs/src/developer/pqs_near_term_final_basis_realization_plan.md` to mark the helper implemented and clarify explicit-origin versus off-origin factor sources.

Validation:
- `julia --project=. -e 't = @elapsed include("test/nested/pqs_direct_retained_final_h1_runtests.jl"); println("elapsed_s=", t)'`
  - PASS: 31 tests, elapsed_s=61.486079042
  - Note: this included package precompilation.
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - PASS: `load ok`
- `git diff --check`
  - PASS

Diff:
- `src/pqs_multilayer_shell_source_plan.jl`
- `test/nested/pqs_direct_retained_final_h1_runtests.jl`
- `docs/src/developer/pqs_near_term_final_basis_realization_plan.md`
- Total: 298 insertions, 52 deletions

Deletion/shrinkage report:
- Became unnecessary:
  - the test-local `_pqs_h1_support_nuclear_matrix` helper.
- Deleted/simplified:
  - removed `_pqs_h1_support_nuclear_matrix`;
  - routed the H1 gate's support nuclear matrix through `pqs_multilayer_support_electron_nuclear_by_center_matrices`;
  - kept the existing H1 gate assertion shape rather than adding a new metadata test suite.
- No new test file was added. The existing H1 gate now exercises the live helper through the endpoint workflow.
- Remaining stale/duplicate surfaces:
  - raw-support H1 probes and old fixed-block/WL nuclear matrices remain oracle/reference only;
  - future cleanup can look for probe-local support nuclear assembly that still duplicates the route-owned helper.

-- repo-doer@macmini
