Purpose:

Perform a caller-driven retirement audit for
`_pqs_current_route_safe_term_matrices(...)` and related current-route
safe-term authority comparison helpers.

Why now:

The final-basis PQS H1 route has been proven by probe. The docs now quarantine
`_pqs_current_route_safe_term_matrices(...)` as oracle/debug vocabulary, but
the source and integration callers still exist. The next cleanup should be
based on actual callers, not broad deletion.

Exact task:

Search for and inspect callers/mentions of:

```text
_pqs_current_route_safe_term_matrices
_pqs_current_route_safe_term_authority_comparison
current_route_safe_term
safe_term_authority
shell_realized_pqs_pairs_use_oracle_not_algorithm
```

Classify each live caller as one of:

```text
delete_now
quarantine_oracle
keep_live_kernel_contract
keep_paper/benchmark/reference
needs_manager_decision
```

Then make only safe cleanup edits:

- delete or simplify obsolete tests/helpers only if they no longer protect a
  live kernel, benchmark, or oracle comparison;
- otherwise, update docs/ledger/comments so the remaining surface is clearly
  oracle/debug and not a final-basis route authority;
- do not add new tests;
- do not touch the final-basis H1 route implementation;
- do not implement IDA, density-density, RHF, drivers, exports, or artifacts.

Likely files to inspect:

- `src/CartesianContractedParentMetrics.jl`
- `test/nested/bond_aligned_diatomic_high_order_recipe_opt_in_source_construction_integration_runtests.jl`
- `docs/src/developer/cartesian_route_retirement_ledger.md`
- `docs/src/developer/pqs_source_box_operator_framework.md`
- `docs/src/developer/raw_product_source_retained_transform_policy.md`
- `docs/src/developer/projected_q_shell_policy.md`

Decision rule:

If a caller is the only coverage for an old paper/benchmark/reference route,
do not delete it. Quarantine it in wording and report the deletion blocker.

If a caller only preserves old helper vocabulary now replaced by the explicit
final-basis H1 path, delete or simplify it.

Validation:

- if a test file is edited, run that focused test only;
- `julia --project=. -e 'using GaussletBases; println("load ok")'`;
- `git diff --check`.

Deletion/shrinkage report required:

- exact caller classification table;
- what old code/test/doc wording became unnecessary;
- what was deleted, simplified, or quarantined;
- if nothing was deleted, why no existing surface was safely obsolete;
- whether any new doc/log artifact was added and why it earned its cost;
- remaining stale or duplicate surfaces to retire next.

Report back:

- write `.agent_handoffs/response.029.md.tmp`, then atomically rename to
  `.agent_handoffs/response.029.md`;
- also write the curated copy to
  `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.029.md`;
- include files changed;
- include validation run;
- include deletion/shrinkage report;
- sign `-- repo-doer@macmini`.

-- repo-manager@macmini
