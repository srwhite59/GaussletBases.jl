Purpose:

Add the smallest retained PQS source-mode global one-body matrix layer.

Why now:

Passes 003 and 004 made retained PQS source-mode overlap/kinetic pair blocks
available through the normal selector/batch surface. The next source-box-first
step is to assemble a one-unit self-pair retained matrix from those retained
blocks.

This is the intended progression:

```text
raw source-space block
-> retained source-mode block
-> one-unit retained source-mode global matrix
```

This is still not shell realization and not final shell-realized PQS.

Exact task:

Add a narrow helper in `CartesianPairBlockMaterialization` for a one-unit
self-pair retained source-mode matrix.

Recommended API shape:

```text
pqs_retained_source_one_body_matrix(batch_result)
```

or, if cleaner:

```text
pqs_retained_source_one_body_matrix(plan, term; overlap_1d, kinetic_1d = nothing)
```

Required behavior:

- accept only retained PQS source-mode overlap/kinetic batch results;
- require exactly one materialized self-pair result for this pass;
- require zero or more skipped records only if they are not PQS-ready retained
  source-mode self-pairs;
- use the retained block directly as the dense global retained matrix;
- report retained dimension equal to the block size;
- report `matrix_space = :retained_pqs_source_modes` or equivalent;
- preserve nonclaims: no shell realization, no Lowdin, no electron-nuclear,
  no IDA, no Hamiltonian, no driver, no exports/artifacts.

If a natural compact object exists locally, use it. Do not create a broad route
object or a large report field cloud.

Trust boundary:

No multi-unit placement unless the needed retained column ranges already exist
and can be used without new framework. No shell projection, Lowdin cleanup,
support-row contraction, electron-nuclear, density-density, IDA, Hamiltonian
assembly, RHF, driver adoption, exports, artifacts, full-parent fallback, or
direct Cartesian fallback.

Test policy:

Extend the existing compact PQS source-pair contract test. Check only:

- one-unit self-pair retained overlap global matrix equals the retained overlap
  block;
- one-unit self-pair retained kinetic global matrix equals the retained kinetic
  block;
- dimensions match retained counts;
- matrix is finite and symmetric for the self-pair fixture;
- shell realization and Lowdin remain false.

Do not add broad metadata checks. Do not add a new test file unless the existing
file becomes genuinely less clear.

Deletion/shrinkage report required:

- what old code, test, metadata, or compatibility path became unnecessary;
- what was deleted or simplified;
- if nothing was deleted, why no existing surface was made obsolete yet;
- whether any new test replaces/shrinks older coverage or is genuinely new
  live-contract coverage;
- any remaining stale or duplicate surfaces to retire next.

Validation:

- `julia --project=. test/nested/cartesian_pair_block_materialization_contract_runtests.jl`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git diff --check`

Report back:

- write `.agent_handoffs/response.005.md.tmp`, then atomically rename to
  `.agent_handoffs/response.005.md`;
- also write the curated copy to
  `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.005.md`;
- include files changed;
- include retained global matrix evidence;
- include validation run;
- include deletion/shrinkage report;
- sign `-- repo-doer@macmini`.

-- repo-manager@macmini
