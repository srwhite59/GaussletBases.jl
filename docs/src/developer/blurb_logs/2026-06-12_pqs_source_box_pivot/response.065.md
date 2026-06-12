Implemented pass 065.

Changes:
- Extended the generic retained PQS one-body selector in
  `src/cartesian_pair_block_materialization/pqs_source_one_body.jl` to accept:
  - `:position_x`, `:position_y`, `:position_z`
  - `:x2_x`, `:x2_y`, `:x2_z`
- Routed those generic retained moment requests through the pass-064 direct
  retained wrappers. The path does not materialize full raw source matrices.
- Extended `pqs_source_pair_retained_one_body_blocks(...)` and
  `pqs_retained_source_one_body_matrix(plan, term; ...)` with `position_1d`
  and `x2_1d` keyword forwarding.
- Allowed `pqs_retained_source_one_body_matrix(batch_result)` to accept square
  retained self-pair position/x2 terms.
- Did not add exports; the named retained moment wrappers remain
  qualified/internal.
- Updated
  `docs/src/developer/pqs_near_term_final_basis_realization_plan.md` to note
  that the generic retained selector now covers overlap, kinetic, position, and
  x2 source terms.

Validation:
- `julia --project=. -e 't = @elapsed include("test/nested/cartesian_pair_block_materialization_contract_runtests.jl"); println("elapsed_s=", t)'`
  - Passed.
  - `CartesianPairBlockMaterialization PQS source-pair batch selectors` now
    reports 304 passing checks.
  - Elapsed: 71.93s, including precompilation/startup.
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - Passed.
- `git diff --check`
  - Passed.

Deletion/shrinkage report:
- No old production code became unnecessary. The generic selector is now a thin
  dispatch surface over the direct retained wrappers instead of introducing a
  parallel route.
- No compatibility shim, metadata layer, fixture rule, driver wiring, export,
  artifact, H1/RHF/IDA/density-density, or GTO path was added.
- No standalone test file was added. The new checks extend existing compact
  CPBM retained/PQS contract coverage and protect the live selector/matrix
  contract.
- Remaining stale/duplicate surface to retire next: none identified in this
  seam. The raw-source selector path remains intentionally available as the
  oracle/reference path.

-- repo-doer@macmini
