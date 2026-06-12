Purpose:

Do a targeted PQS cleanup/retirement audit after the final-basis H1 probe.

Why now:

Pass 027 proved the projected-q-shell final-basis one-electron route against a
shell-support oracle to roundoff. That means some older retained-source H1
scratch probes, readiness blockers, and current-route safe-term oracle
vocabulary may now be stale or should be clearly quarantined.

Exact task:

Search active tests/docs/source for stale PQS one-electron readiness surfaces,
especially:

```text
retained-source H1
retained_source_h1
pqs_retained_source
shell_realization_not_materialized
missing_pqs_final_one_electron_solve
_pqs_current_route_safe_term_matrices
current-route safe-term
safe-term oracle
```

Inspect the hits. Then:

- delete or simplify any active test/doc wording that now preserves the old
  retained-source/H1 scratch route as if it were the live path;
- quarantine old current-route safe-term oracle language as oracle/debug only
  if it still has live comparison value;
- do not delete raw source-box block tests, retained rule tests, centered
  nuclear source tests, or compact module-contract tests that protect live
  kernels;
- do not add new tests in this pass;
- do not turn the tmp H1 probe into a permanent gate yet.

If no deletion is safe, write a precise retirement ledger/log note explaining
which surfaces remain and why.

Suggested places to inspect first:

- `test/nested/runtests.jl`
- `test/nested/integration_runtests.jl`
- `test/nested/cartesian_pair_block_materialization_contract_runtests.jl`
- `docs/src/developer/pqs_source_box_operator_framework.md`
- `docs/src/developer/raw_product_source_retained_transform_policy.md`
- `docs/src/developer/projected_q_shell_policy.md`
- `docs/src/developer/cartesian_route_retirement_ledger.md`
- `docs/src/developer/numerical_contracts.md`

Trust boundary:

Cleanup/retirement only. Do not implement a new physics route, new H1 helper,
IDA, density-density, RHF, drivers, exports, or artifacts. Do not add broad
metadata assertions.

Test policy:

No new tests. Run only the smallest validation needed for edited files. If only
docs/logs change, no Julia test is required beyond a load check.

Deletion/shrinkage report required:

- what old code, test, metadata, or compatibility path became unnecessary;
- what was deleted or simplified;
- if nothing was deleted, why no existing surface was made obsolete yet;
- whether any new doc/log artifact was added and why it earned its cost;
- any remaining stale or duplicate surfaces to retire next.

Validation:

- smallest touched test if a test file is edited;
- `julia --project=. -e 'using GaussletBases; println("load ok")'`;
- `git diff --check`.

Report back:

- write `.agent_handoffs/response.028.md.tmp`, then atomically rename to
  `.agent_handoffs/response.028.md`;
- also write the curated copy to
  `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.028.md`;
- include files changed;
- include what was deleted/simplified/quarantined;
- include remaining stale surfaces;
- include validation run;
- include deletion/shrinkage report;
- sign `-- repo-doer@macmini`.

-- repo-manager@macmini
