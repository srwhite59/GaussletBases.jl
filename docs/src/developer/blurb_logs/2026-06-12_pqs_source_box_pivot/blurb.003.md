Purpose:

Implement the first retained PQS safe one-body block contraction from existing
raw source-space PQS/PQS blocks.

Why now:

Pass 002 added the source-mode boundary retained rule. The next source-box-first
operator step is:

```text
raw source-space one-body block O_source
retained source-mode selector columns
-> O_retained = T_left' * O_source * T_right
```

For the current selector transform, this is equivalent to:

```text
O_retained = O_source[left_columns, right_columns]
```

This is still source-mode retained block construction. It is not shell
realization, not shell-row contraction, and not Lowdin.

Exact task:

Add a narrow retained-block helper in `CartesianPairBlockMaterialization` for
PQS/PQS safe one-body terms.

Primary code surfaces:

```text
src/cartesian_pair_block_materialization/pqs_source_safe_terms.jl
test/nested/cartesian_pair_block_materialization_contract_runtests.jl
```

Use existing raw source-space block helpers as input:

```text
pqs_source_pair_overlap_block(...)
pqs_source_pair_kinetic_block(...)
```

Recommended helper shape:

```text
pqs_source_pair_retained_one_body_block(source_result, left_rule, right_rule)
```

or, if cleaner for existing records:

```text
pqs_source_pair_retained_overlap_block(record; overlap_1d)
pqs_source_pair_retained_kinetic_block(record; overlap_1d, kinetic_1d)
```

The retained helper should:

- require raw source-space source result materialized;
- require left/right PQS retained rules from transform-contract metadata;
- use `retained_column_indices(rule)` for left/right selectors;
- return a block with shape `(left_retained_count, right_retained_count)`;
- report `block_space = :retained_pqs_source_modes` or similarly clear label;
- report source-space input was used;
- report shell realization not materialized;
- report Lowdin not used;
- report no final Hamiltonian, IDA, export, or artifact data.

Start with overlap and kinetic. Do not add position/x2 unless it is essentially
free and does not broaden the test burden.

Trust boundary:

No shell projection, Lowdin, support-row contraction, electron-nuclear,
density-density, IDA, Hamiltonian, RHF, global route assembly, driver adoption,
exports, artifacts, full-parent fallback, or direct Cartesian fallback.

Test policy:

One compact module-contract test is allowed. It should check:

- raw overlap or kinetic source block shape;
- retained overlap or kinetic block shape;
- retained block equals `source_block[left_columns, right_columns]`;
- no shell realization and no Lowdin;
- source-space raw result remains available as input, not route authority for
  final shell realization.

Avoid broad metadata checks and driver tests.

Validation:

- run `julia --project=. test/nested/cartesian_pair_block_materialization_contract_runtests.jl`;
- run any directly affected raw-source/retained-transform test if needed;
- run `julia --project=. -e 'using GaussletBases; println("load ok")'`;
- run `git diff --check`.

Deletion/shrinkage report required:

- what old code, test, metadata, or compatibility path became unnecessary;
- what was deleted or simplified;
- if nothing was deleted, why no existing surface was made obsolete yet;
- whether any new test replaces/shrinks older coverage or is genuinely new
  live-contract coverage;
- any remaining stale or duplicate surfaces to retire next.

Report back:

- write `.agent_handoffs/response.003.md.tmp`, then atomically rename to
  `.agent_handoffs/response.003.md`;
- also write the curated copy to
  `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.003.md`;
- include files changed;
- include retained block shape/equality evidence;
- include validation run;
- include deletion/shrinkage report;
- sign `-- repo-doer@macmini`.

-- repo-manager@macmini
