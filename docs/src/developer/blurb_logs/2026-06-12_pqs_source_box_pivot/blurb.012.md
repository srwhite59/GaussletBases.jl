Purpose:

Connect the generated PQS centered Gaussian source factors to the existing
uncharged by-center electron-nuclear source block.

Why now:

Pass 007 already materialized an uncharged by-center PQS source
electron-nuclear block from caller-supplied source-mode Gaussian factors. Pass
010 projected support-space Gaussian factors into source-mode coordinates. Pass
011 generated those support-space Gaussian factors from explicit axis layers.
The next useful source-box step is the direct composition:

```text
axis layers + Coulomb expansion + center
-> centered Gaussian source-mode factors
-> uncharged by-center electron-nuclear source block
```

Exact task:

Add a narrow CPBM convenience helper, or a better local equivalent:

```text
pqs_source_pair_centered_electron_nuclear_by_center_block(
    record;
    axis_layers,
    coulomb_expansion,
    center_record,
)
```

It should:

- call `pqs_source_pair_centered_gaussian_factor_terms_1d(...)`;
- pass the result to `pqs_source_pair_electron_nuclear_by_center_block(...)`;
- preserve the current convention: nuclear charge recorded but not applied,
  centers not summed, by-center block uncharged;
- not call CCPM source-box wrappers;
- not build shell realization, Lowdin cleanup, IDA data, Hamiltonians, global
  routes, drivers, exports, or artifacts.

If it is useful and does not broaden the pass, add the retained contraction
wrapper too:

```text
pqs_source_pair_retained_centered_electron_nuclear_by_center_block(...)
```

That wrapper should only compose the centered source block with the existing
retained source-mode contraction. Do not add a global H1 route in this pass.

Test policy:

Use the existing CPBM source-pair contract file. Add only compact assertions
that the centered helper matches the already-tested supplied-factor path using
the centered factors from pass 011. If the retained wrapper is added, check it
matches `pqs_source_pair_retained_one_body_block(source_result)`.

Do not add a new test file. Do not add broad metadata checks beyond the live
convention flags: by-center, charge recorded/not applied, centers not summed,
source/retained block shape, and no shell/Lowdin/IDA/Hamiltonian claims.

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

- write `.agent_handoffs/response.012.md.tmp`, then atomically rename to
  `.agent_handoffs/response.012.md`;
- also write the curated copy to
  `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.012.md`;
- include implementation or exact blocker;
- include validation run;
- include deletion/shrinkage report;
- sign `-- repo-doer@macmini`.

-- repo-manager@macmini
