Purpose:

Implement the minimal PQS source-space electron-nuclear by-center block from
caller-supplied Gaussian factor matrices.

Why now:

Pass 006 established the convention: the new CPBM PQS path should build a
negative unit-charge by-center attraction matrix using `-c_t`, record center
charge, and leave charge application/summing for later diagnostic or
Hamiltonian assembly.

Do not implement source-axis Gaussian factor generation yet. This pass should
only implement the 3D product contraction when the caller already supplies
term-first axis factor matrices.

Exact task:

Add a narrow helper in `CartesianPairBlockMaterialization`:

```text
pqs_source_pair_electron_nuclear_by_center_block(
    record;
    coulomb_expansion,
    center_record,
    gaussian_factor_terms_1d,
)
```

Equivalent argument names are fine if they match local style.

Required behavior:

- require a ready PQS/PQS source-pair record;
- require term-first axis factor arrays with shape
  `(nterms, left_axis_count, right_axis_count)` for x/y/z;
- require `length(coulomb_expansion.coefficients) == nterms`;
- build the raw source-space block as:

```text
sum_t (-coefficients[t]) * Gx[t] * Gy[t] * Gz[t]
```

- use the same source-mode ordering/product fill conventions as the existing
  PQS safe one-body helpers;
- return a `PairBlockMaterializationResult` in raw product source mode;
- metadata must say:
  - by-center data;
  - nuclear charge recorded;
  - nuclear charge not applied;
  - centers not summed;
  - uncharged by-center convention;
  - no shell realization;
  - no Lowdin;
  - no IDA;
  - no Hamiltonian/driver/export/artifact.

Retained contraction:

Add the smallest retained wrapper only if it is just:

```text
pqs_source_pair_retained_one_body_block(source_result)
```

Do not add a global matrix, H1 diagnostic, or source-axis Gaussian factor
generator in this pass.

Code surfaces:

```text
src/cartesian_pair_block_materialization/pqs_source_safe_terms.jl
src/cartesian_pair_block_materialization/pqs_source_one_body.jl
src/cartesian_pair_block_materialization/CartesianPairBlockMaterialization.jl
test/nested/cartesian_pair_block_materialization_contract_runtests.jl
```

If a small new CPBM file is cleaner than expanding `pqs_source_safe_terms.jl`,
that is acceptable. Do not create a broad route module.

Test policy:

Extend the existing compact PQS source-pair contract test. Use synthetic
term-first Gaussian factor arrays and coefficients. Check only:

- source nuclear block equals a manual two-term product contraction;
- source term/block space are correct;
- charge is recorded but not applied;
- centers are not summed;
- retained nuclear block equals source block selected by retained columns;
- no shell realization and no Lowdin.

Do not use old CCPM as production authority. Do not add broad metadata checks
or a new test file unless the existing test becomes unclear.

Trust boundary:

No source-axis Gaussian factor generator, shell projection, Lowdin cleanup,
support-row contraction, density-density, IDA, Hamiltonian assembly, H1 solve,
RHF, driver adoption, exports, artifacts, full-parent fallback, direct
Cartesian fallback, or old CCPM wrapper adoption.

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

- write `.agent_handoffs/response.007.md.tmp`, then atomically rename to
  `.agent_handoffs/response.007.md`;
- also write the curated copy to
  `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.007.md`;
- include files changed;
- include source/retained nuclear block evidence;
- include validation run;
- include deletion/shrinkage report;
- sign `-- repo-doer@macmini`.

-- repo-manager@macmini
