Implemented the retained PQS source-mode safe one-body selector/batch layer for overlap and kinetic.

Files changed:

- `src/cartesian_pair_block_materialization/CartesianPairBlockMaterialization.jl`
- `src/cartesian_pair_block_materialization/pqs_source_one_body.jl`
- `test/nested/cartesian_pair_block_materialization_contract_runtests.jl`

What changed:

- Added `pqs_source_pair_retained_one_body_block(record, term; overlap_1d, kinetic_1d = nothing)`.
- Added `pqs_source_pair_retained_one_body_blocks(plan, term; overlap_1d, kinetic_1d = nothing)`.
- Exported `pqs_source_pair_retained_one_body_blocks`.
- The retained selector delegates through the existing raw source-space one-body selector, then applies the retained source-mode contraction from pass 003.
- Supported retained selector terms in this pass:
  - `:overlap`
  - `:kinetic`
- Position and x2 were not added; adding them would broaden this pass beyond the requested overlap/kinetic selector layer.

Selector/batch evidence:

- Retained overlap selector result matches `pqs_source_pair_retained_overlap_block(...)`.
- Retained kinetic selector result matches `pqs_source_pair_retained_kinetic_block(...)`.
- Retained overlap batch materialized the same ready pair count and skipped count as the raw overlap batch.
- Retained kinetic batch materialized the same ready pair count as the raw kinetic batch.
- Retained batch records report `block_space = :retained_pqs_source_modes`.
- Retained selector/batch records keep shell realization and Lowdin cleanup false.
- The retained path still does not claim final shell-realized pair blocks, Hamiltonian data, IDA data, exports, or artifacts.

Validation run:

- `julia --project=. test/nested/cartesian_pair_block_materialization_contract_runtests.jl` passed.
- `julia --project=. -e 'using GaussletBases; println("load ok")'` passed with `load ok`.
- `git diff --check` passed.

Deletion/shrinkage report:

- No old code was deleted in this pass.
- No existing compatibility path became obsolete yet because this pass wires the retained source-mode contraction into the selector/batch surface but does not replace downstream shell-realization/readiness or global route consumers.
- The new test coverage is genuinely new live-contract coverage for the retained source-mode selector/batch surface. It extends the existing PQS source-pair batch selector test rather than adding a new file.
- Remaining stale/duplicate surfaces to retire next:
  - raw source-space selector tests can be reduced later once retained source-mode blocks become the active route consumer;
  - shell/support-row source-box reference helpers remain candidates for oracle-only quarantine after a retained source-mode route layer consumes these retained blocks directly.

-- repo-doer@macmini
