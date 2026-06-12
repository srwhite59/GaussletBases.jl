Purpose:

Project supplied parent/source-axis Gaussian term matrices into PQS source-mode
Gaussian factor arrays using the materialized axis transform facts.

Why now:

Pass 009 plumbed materialized source-axis coefficient transforms to the CPBM
preflight record. Pass 007 already consumes term-first source-mode Gaussian
factor arrays for electron-nuclear blocks. The missing narrow layer is:

```text
parent/source-axis Gaussian term matrices
+ left/right source-axis transforms
-> source-mode Gaussian factor arrays
```

This is projection only. Do not implement analytic Gaussian matrix generation
yet.

Exact task:

Add a CPBM helper that, for one PQS/PQS source-pair record, takes supplied
term-first axis Gaussian matrices in source-support coordinates and projects
them to source-mode coordinates:

```text
G_source[t] = C_left' * G_parent[t] * C_right
```

Recommended API shape:

```text
pqs_source_pair_gaussian_factor_terms_1d(
    record;
    gaussian_factor_terms_axis,
)
```

where `gaussian_factor_terms_axis` is an x/y/z tuple or named tuple of
term-first arrays with shape:

```text
(nterms, left_source_support_count, right_source_support_count)
```

The helper should return x/y/z term-first arrays with shape:

```text
(nterms, left_source_mode_dim, right_source_mode_dim)
```

Required behavior:

- require materialized left/right `AxisSourceTransformFact`s on the record;
- validate coefficient matrices and term matrices by axis;
- project each term with `C_left' * G * C_right`;
- preserve no-shell/no-Lowdin/no-IDA/no-Hamiltonian nonclaims in any summary;
- do not call CCPM `_pqs_pqs_source_box_*nuclear_attraction*` wrappers;
- do not apply nuclear charge or build electron-nuclear blocks in this helper.

It is acceptable to add a thin convenience wrapper that feeds the projected
arrays into `pqs_source_pair_electron_nuclear_by_center_block(...)`, but only
if it is trivial and does not add new route claims.

Test policy:

Extend the existing compact PQS source-pair contract test. Use synthetic
axis-transform matrices and synthetic parent/source-axis term matrices. Check:

- projected x/y/z arrays equal manual `C_left' * G * C_right` for one or two
  terms;
- shape is `(nterms, left_source_mode_dim, right_source_mode_dim)`;
- missing/not-materialized axis facts block or throw clearly;
- if a nuclear convenience wrapper is added, its source block matches the
  pass-007 supplied-factor block fed with the projected arrays.

Do not add a new test file, physics test, H1 diagnostic, broad metadata checks,
or analytic Gaussian generation.

Trust boundary:

No analytic Gaussian factor generation from centers/exponents, shell
projection, Lowdin cleanup, density-density, IDA, Hamiltonian assembly, H1
solve, RHF, driver adoption, exports, artifacts, full-parent fallback, direct
Cartesian fallback, or old CCPM wrapper adoption.

Deletion/shrinkage report required:

- what old code, test, metadata, or compatibility path became unnecessary;
- what was deleted or simplified;
- if nothing was deleted, why no existing surface was made obsolete yet;
- whether any new test/artifact was added and why it earned its cost;
- any remaining stale or duplicate surfaces to retire next.

Validation:

- `julia --project=. test/nested/cartesian_pair_block_materialization_contract_runtests.jl`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git diff --check`

Report back:

- write `.agent_handoffs/response.010.md.tmp`, then atomically rename to
  `.agent_handoffs/response.010.md`;
- also write the curated copy to
  `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.010.md`;
- include files changed;
- include projection evidence;
- include validation run;
- include deletion/shrinkage report;
- sign `-- repo-doer@macmini`.

-- repo-manager@macmini
