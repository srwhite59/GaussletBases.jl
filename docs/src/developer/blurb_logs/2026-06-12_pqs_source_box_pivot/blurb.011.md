Purpose:

Generate PQS source-support Gaussian term matrices from explicit axis layers,
then project them to source-mode Gaussian factor arrays.

Why now:

Pass 010 can project supplied support-space Gaussian term matrices through
materialized source-axis transforms. The remaining pre-H1 source-factor step is
to produce those support-space matrices from explicit axis-layer objects using
the repo's low-level analytic Gaussian factor API:

```text
gaussian_factor_matrices(layer; exponents, center)
```

Do not call the old CCPM positive-gaussian/nuclear wrapper helpers.

Exact task:

Add a narrow CPBM helper, or block precisely if the low-level API is not
available from this module without broad dependency changes.

Recommended API shape:

```text
pqs_source_pair_centered_gaussian_factor_terms_1d(
    record;
    axis_layers,
    coulomb_expansion,
    center_record,
)
```

or equivalent local naming.

Required behavior if implemented:

- require explicit axis layers for x/y/z;
- require `coulomb_expansion.exponents`;
- extract center coordinate from `center_record`;
- for each axis, call the low-level analytic Gaussian factor API to obtain one
  matrix per exponent;
- slice those matrices to the left/right source intervals from the pair record
  or axis transform facts;
- build term-first support arrays with shape
  `(nterms, left_support_count, right_support_count)`;
- call `pqs_source_pair_gaussian_factor_terms_1d(...)` from pass 010 to project
  to source-mode arrays;
- preserve no-shell/no-Lowdin/no-IDA/no-Hamiltonian nonclaims in any summary.

If the explicit axis-layer object shape is not available in the existing test
fixtures, use a `tmp/work` probe or return a precise blocker. Do not invent a
fake analytic backend.

Trust boundary:

No electron-nuclear block assembly unless it is a trivial optional wrapper
around pass 007 and clearly uses the projected arrays. No H1 solve, RHF,
density-density, IDA, shell realization, Lowdin cleanup, Hamiltonian assembly,
driver adoption, exports, artifacts, full-parent fallback, direct Cartesian
fallback, or old CCPM wrapper adoption.

Test policy:

Prefer one compact test only if there is an existing small analytic axis-layer
fixture. Check:

- generated/projected arrays match direct projection of
  `gaussian_factor_matrices(layer; exponents, center)` sliced to the source
  intervals;
- term count and shapes are correct;
- no CCPM wrapper path is claimed.

If no suitable fixture exists, do not add a synthetic fake-backend test. Report
the blocker and leave code unchanged.

Deletion/shrinkage report required:

- what old code, test, metadata, or compatibility path became unnecessary;
- what was deleted or simplified;
- if nothing was deleted, why no existing surface was made obsolete yet;
- whether any new test/artifact was added and why it earned its cost;
- any remaining stale or duplicate surfaces to retire next.

Validation:

- if production code/tests changed:
  `julia --project=. test/nested/cartesian_pair_block_materialization_contract_runtests.jl`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git diff --check`

Report back:

- write `.agent_handoffs/response.011.md.tmp`, then atomically rename to
  `.agent_handoffs/response.011.md`;
- also write the curated copy to
  `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.011.md`;
- include implementation or exact blocker;
- include validation run;
- include deletion/shrinkage report;
- sign `-- repo-doer@macmini`.

-- repo-manager@macmini
