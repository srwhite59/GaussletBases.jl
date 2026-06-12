Purpose:
  Implement the route-owned multi-layer PQS support electron-nuclear by-center
  helper under the convention documented in pass 072, then shrink the H1 gate's
  local nuclear helper.

Context:
  The convention is:

  - helper name target:
    `pqs_multilayer_support_electron_nuclear_by_center_matrices`;
  - input plan must be an available `:pqs_multilayer_shell_source_plan`;
  - return one support-space matrix record per center;
  - support ordering is `plan.core_support_states` followed by
    `plan.shell_support_states`;
  - each uncharged center matrix is negative unit-charge attraction `-1/r_A`;
  - Hamiltonian assembly later applies charge by adding `Z_A * V_A`;
  - helper must not apply `Z_A`;
  - helper must not sum centers;
  - no final-basis transfer, H1, IDA, density-density, RHF, driver wiring,
    exports, or artifacts.

Implementation:
  Add the helper near the other multi-layer support helpers.

  Suggested signature:

  ```julia
  pqs_multilayer_support_electron_nuclear_by_center_matrices(
      plan;
      coulomb_expansion,
      center_records,
      axis_layers = nothing,
      gaussian_factor_terms_by_center = nothing,
  )
  ```

  Requirements:

  - If `gaussian_factor_terms_by_center` is supplied, treat it as explicit
    centered factor terms for the supplied centers. This is sufficient for the
    current origin H1 gate when it passes `pgdg_intermediate.gaussian_factor_terms`
    for the origin center.
  - If explicit factors are not supplied, use `axis_layers`, expansion
    exponents, and center locations to build centered Gaussian factor matrices
    with the same low-level convention used by
    `pqs_source_pair_centered_gaussian_factor_terms_1d(...)`.
  - For each center and each Gaussian term, accumulate
    `-coefficient[t] * Gx_t * Gy_t * Gz_t` over the complete support state list.
  - Record center key/index/location and nuclear charge in each result record,
    but set `nuclear_charge_applied = false` and `centers_summed = false`.

  If off-origin support needs more than can be done cleanly in this pass, it is
  acceptable to implement the explicit-factor path first and return a precise
  blocker for missing centered axis-layer support. Do not fake off-origin
  support with origin PGDG factors.

H1 gate update:
  Update `test/nested/pqs_direct_retained_final_h1_runtests.jl` to use the new
  helper for its origin nuclear matrix, passing explicit origin
  `fixture.bundle7.pgdg_intermediate.gaussian_factor_terms` as the factor source.
  Delete `_pqs_h1_support_nuclear_matrix` if it becomes unnecessary.

Validation:
  - `julia --project=. -e 't = @elapsed include("test/nested/pqs_direct_retained_final_h1_runtests.jl"); println("elapsed_s=", t)'`
  - `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - `git diff --check`

Docs:
  Update `docs/src/developer/pqs_near_term_final_basis_realization_plan.md` to
  mark the origin/explicit-factor by-center support helper status accurately.
  If off-origin centered axis-layer support is still blocked, say so explicitly.

Do not:
  - change H1 reference values;
  - add RHF/IDA/density-density features;
  - change fixture-rule policy;
  - add exports or driver wiring;
  - edit ignored `tmp/work` probes;
  - use fixed-block/WL matrices as authority.

Deletion/shrinkage report required:
  - what old code, test, metadata, or compatibility path became unnecessary;
  - what was deleted or simplified;
  - if nothing was deleted, why no existing surface was made obsolete;
  - whether any new test replaces/shrinks older coverage or is genuinely new
    live-contract coverage;
  - any remaining stale or duplicate surfaces to retire next.

-- repo-manager@macmini
