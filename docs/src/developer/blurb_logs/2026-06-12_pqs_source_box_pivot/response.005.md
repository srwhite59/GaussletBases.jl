Implemented the first one-unit retained PQS source-mode global one-body matrix helper.

Files changed:

- `src/cartesian_pair_block_materialization/CartesianPairBlockMaterialization.jl`
- `src/cartesian_pair_block_materialization/pqs_source_one_body.jl`
- `test/nested/cartesian_pair_block_materialization_contract_runtests.jl`

What changed:

- Added `pqs_retained_source_one_body_matrix(batch_result)`.
- Added `pqs_retained_source_one_body_matrix(plan, term; overlap_1d, kinetic_1d = nothing)`.
- Exported `pqs_retained_source_one_body_matrix`.
- The helper accepts retained PQS source-mode overlap/kinetic batch results only.
- It requires exactly one materialized retained self-pair block.
- It rejects skipped records that are ready PQS retained self-pairs.
- It uses the retained self-pair block directly as the dense retained source-mode matrix.
- It reports `matrix_space = :retained_pqs_source_modes`.
- It keeps the nonclaims: no shell realization, no Lowdin, no electron-nuclear, no density-density, no IDA, no Hamiltonian, no driver route, no exports, no artifacts.

Retained global matrix evidence:

- Symmetric one-unit retained overlap batch materializes one retained self-pair block.
- `pqs_retained_source_one_body_matrix(...)` returns that retained overlap block as the matrix by identity.
- Symmetric one-unit retained kinetic batch materializes one retained self-pair block.
- `pqs_retained_source_one_body_matrix(...)` returns that retained kinetic block as the matrix by identity.
- Matrix retained dimension equals the retained block row count.
- Matrix entries are finite.
- The overlap and kinetic matrix fixtures are symmetric.
- Shell realization and Lowdin cleanup remain false in the matrix result.

Validation run:

- `julia --project=. test/nested/cartesian_pair_block_materialization_contract_runtests.jl` passed.
- `julia --project=. -e 'using GaussletBases; println("load ok")'` passed with `load ok`.
- `git diff --check` passed.

Deletion/shrinkage report:

- No old code was deleted in this pass.
- No existing compatibility path became obsolete yet because this is the first one-unit retained source-mode matrix layer and does not replace the existing shell-realization bridge, final-readiness summaries, or future multi-unit placement.
- The new test coverage is genuinely new live-contract coverage for one-unit retained PQS source-mode matrix construction. It extends the existing PQS source-pair selector test rather than adding a new file.
- Remaining stale/duplicate surfaces to retire next:
  - source-space-only bridge/readiness paths should be narrowed once a retained source-mode route consumer exists above this one-unit matrix helper;
  - old shell/support-row source-box oracle helpers remain candidates for quarantine/deletion after retained source-mode route assembly replaces their test pressure.

-- repo-doer@macmini
