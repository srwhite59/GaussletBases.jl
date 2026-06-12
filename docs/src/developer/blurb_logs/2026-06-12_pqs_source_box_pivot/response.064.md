Implemented pass 064.

Changes:
- Added exported direct retained PQS moment wrappers in `CartesianPairBlockMaterialization`:
  - `pqs_source_pair_retained_position_x_block`
  - `pqs_source_pair_retained_position_y_block`
  - `pqs_source_pair_retained_position_z_block`
  - `pqs_source_pair_retained_x2_x_block`
  - `pqs_source_pair_retained_x2_y_block`
  - `pqs_source_pair_retained_x2_z_block`
- Added one shared internal retained axis-product helper in
  `src/cartesian_pair_block_materialization/pqs_source_safe_terms.jl`.
- Extended the existing compact PQS source-pair contract section in
  `test/nested/cartesian_pair_block_materialization_contract_runtests.jl` to
  compare each new direct retained position/x2 wrapper against the raw
  source-space block plus retained-selector oracle path.
- Updated `docs/src/developer/pqs_near_term_final_basis_realization_plan.md`
  to mark direct retained position/x2 one-body source blocks as present.

Validation:
- `julia --project=. -e 't = @elapsed include("test/nested/cartesian_pair_block_materialization_contract_runtests.jl"); println("elapsed_s=", t)'`
  - Passed: 1140 tests total across the CPBM contract file.
  - The edited `CartesianPairBlockMaterialization PQS source-pair preflight`
    section passed with 422 tests.
  - Elapsed: 71.03s, including precompilation/startup.
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - Passed.
- `git diff --check`
  - Passed.

Deletion/shrinkage report:
- No old production path became unnecessary in this pass. The raw source-space
  block plus retained-selector path remains the intended small-fixture oracle
  for the new direct retained wrappers.
- No compatibility shim was added.
- No broad metadata assertions or new standalone test file were added. The test
  addition is genuinely new live-contract coverage for the six new direct
  retained moment wrappers and is localized to the existing CPBM contract file.
- Remaining stale/duplicate surface to consider later: the generic retained
  term selector still only exposes overlap/kinetic. I did not widen it because
  this blurb asked for named direct retained moment wrappers, and widening the
  generic selector would be a separate API decision.

-- repo-doer@macmini
