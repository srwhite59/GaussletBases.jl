Purpose:
  Correct pass 073 before it is accepted. The explicit origin-factor path is
  useful, but the claimed off-origin `axis_layers` path appears to return the
  wrong factor shape.

Issue to check:
  In `src/pqs_multilayer_shell_source_plan.jl`, the axis-layer path currently
  appears to do:

  ```julia
  factors = ntuple(axis -> gaussian_factor_matrices(...), 3)
  return (factors[1], factors[2], factors[3], :centered_axis_layers)
  ```

  But `gaussian_factor_matrices(...)` returns a collection of 2D matrices, while
  `_pqs_multilayer_validate_factor_terms(...)` expects each axis term object to
  be a term-first 3D array with shape `(nterms, n, n)`.

Task:
  Fix the pass-073 helper before acceptance.

  Preferred fix:
  - Convert the axis-layer `gaussian_factor_matrices(...)` result for each axis
    into a term-first `Array{Float64,3}` matching the explicit
    `pgdg_intermediate.gaussian_factor_terms` shape.
  - Validate term count and matrix shape.
  - Add a focused test that exercises the `axis_layers` path, ideally with an
    origin center where it can be compared against the explicit
    `pgdg_intermediate.gaussian_factor_terms` path or against the final H1
    gate's existing nuclear matrix.

  Acceptable fallback if this is not local:
  - Narrow the helper/docs to say only explicit `gaussian_factor_terms_by_center`
    is implemented.
  - Return a precise blocker/error for missing axis-layer centered factor
    support.
  - Remove any claim that off-origin centered axis-layer support is implemented.

Do not:
  - commit the current pass-073 state as accepted without this correction;
  - add H1/RHF/IDA/density-density features;
  - change fixture-rule policy;
  - add exports or driver wiring;
  - use origin PGDG factors as off-origin authority.

Validation:
  - focused H1 test or the smallest test covering the axis-layer helper path;
  - `julia --project=. -e 'using GaussletBases; println("load ok")'`;
  - `git diff --check`.

Report:
  - whether the axis-layer path is now implemented and tested, or explicitly
    blocked;
  - any numerical comparison used for the test;
  - deletion/shrinkage report as usual.

-- repo-manager@macmini
