Purpose:

Wire retained PQS source-mode safe one-body blocks through the existing
one-body selector/batch surface.

Why now:

Pass 003 added single-record retained source-mode contraction for overlap and
kinetic. The next small source-box-first step is to let callers request
retained PQS source-mode safe-term blocks at the same selector/batch level that
already exists for raw source-space blocks.

This should make the intended consumer path clearer:

```text
PQS/PQS source-pair preflight
-> raw source-space safe one-body block
-> retained source-mode safe one-body block
-> later global retained source-mode route
```

It should not jump to shell realization or final PQS blocks.

Exact task:

Add a narrow retained one-body selector/batch layer in
`CartesianPairBlockMaterialization`.

Primary code surfaces:

```text
src/cartesian_pair_block_materialization/pqs_source_one_body.jl
src/cartesian_pair_block_materialization/pqs_source_safe_terms.jl
test/nested/cartesian_pair_block_materialization_contract_runtests.jl
```

Recommended API shape:

```text
pqs_source_pair_retained_one_body_block(record, term; overlap_1d, ...)
pqs_source_pair_retained_one_body_blocks(plan, term; overlap_1d, ...)
```

The selector should delegate to existing raw source-space one-body helpers and
then call the retained source-mode contraction added in pass 003.

Required supported terms for this pass:

```text
:overlap
:kinetic
```

Optional only if very small and directly reuses existing raw helpers:

```text
:position_x, :position_y, :position_z
:x2_x, :x2_y, :x2_z
```

If adding position/x2 would broaden the test or implementation materially, do
not add them in this pass.

Trust boundary:

No shell projection, Lowdin cleanup, support-row contraction, electron-nuclear,
density-density, IDA, Hamiltonian assembly, RHF, global route assembly, driver
adoption, exports, artifacts, full-parent fallback, or direct Cartesian
fallback.

Test policy:

Keep the test compact. Prefer extending the existing PQS source-pair batch
selector test rather than adding a new file.

Check only:

- retained selector returns the same block as the term-specific retained helper;
- retained batch materializes the same ready pair count as the raw source batch;
- retained batch record has `block_space = :retained_pqs_source_modes`;
- shell realization and Lowdin remain false.

Do not add broad report-field coverage.

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

- write `.agent_handoffs/response.004.md.tmp`, then atomically rename to
  `.agent_handoffs/response.004.md`;
- also write the curated copy to
  `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.004.md`;
- include files changed;
- include selector/batch evidence;
- include validation run;
- include deletion/shrinkage report;
- sign `-- repo-doer@macmini`.

-- repo-manager@macmini
