Pass 074 response: PQS support electron-nuclear axis-layer correction

Summary:
- Corrected the pass-073 axis-layer path in `pqs_multilayer_support_electron_nuclear_by_center_matrices`.
- Added `_pqs_multilayer_term_first_factor_array(...)` so the result of `gaussian_factor_matrices(...)` is converted from a collection of 2D matrices into the term-first `Array{Float64,3}` shape expected by the support electron-nuclear contraction.
- Kept the explicit origin-factor path from pass 073.
- Added one focused assertion inside the existing H1 gate:
  - origin axis-layer centered factors via `gaussian_factor_matrices(...; center = 0.0)` produce the same support nuclear matrix as explicit `pgdg_intermediate.gaussian_factor_terms`;
  - comparison tolerance: `atol = 1.0e-12`, `rtol = 0.0`.
- Updated the docs note to keep the axis-layer support claim accurate.

Axis-layer status:
- Implemented and tested for the origin comparison path.
- The helper now validates term-first factor count and matrix shape before contraction.
- It does not use origin PGDG factors as off-origin authority; off-origin use must supply axis layers so centered factors are generated at each center location.

Validation:
- `julia --project=. -e 't = @elapsed include("test/nested/pqs_direct_retained_final_h1_runtests.jl"); println("elapsed_s=", t)'`
  - PASS: 32 tests, elapsed_s=62.70530675
  - Note: this included package precompilation.
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - PASS: `load ok`
- `git diff --check`
  - PASS

Diff:
- `src/pqs_multilayer_shell_source_plan.jl`
- `test/nested/pqs_direct_retained_final_h1_runtests.jl`
- `docs/src/developer/pqs_near_term_final_basis_realization_plan.md`
- Current combined pass-073/pass-074 diff: 343 insertions, 52 deletions

Deletion/shrinkage report:
- Became unnecessary:
  - the test-local `_pqs_h1_support_nuclear_matrix` helper remains deleted from pass 073.
- Deleted/simplified:
  - no additional deletion beyond the pass-073 shrinkage; pass 074 corrected a live helper path before acceptance.
- Test coverage:
  - one focused assertion was added to replace a false implementation claim with live-contract coverage for the axis-layer factor-shape bug.
- Remaining stale/duplicate surfaces:
  - raw-support H1 probes and old fixed-block/WL nuclear matrices remain oracle/reference only;
  - future cleanup can retire any other probe-local support nuclear assembly that duplicates the route-owned helper.

-- repo-doer@macmini
