Purpose:

Run the slow nested validation now that the route-driver stack overflow has
been fixed, and classify any remaining blocker.

Why now:

Pass 029 edited an integration test by deleting helper-vocabulary assertions,
but the proper slow nested harness failed before reaching that file. Pass 030
fixed the earlier route-driver stack overflow. We now need to know whether the
pass 029 cleanup is validated or whether the slow harness has another live
blocker.

Exact task:

Run:

```text
TMPDIR=/private/tmp GAUSSLETBASES_TEST_GROUPS=nested GAUSSLETBASES_SLOW_TESTS=1 julia --project=. test/runtests.jl
```

Use a concrete Julia binary if plain `julia` hits the Juliaup lockfile issue in
this sandbox.

Decision rule:

- If the harness passes, report it.
- If it fails before reaching
  `bond_aligned_diatomic_high_order_recipe_opt_in_source_construction_integration_runtests.jl`,
  classify the failure as an unrelated validation blocker and do not edit code
  unless the fix is as narrow and obvious as pass 030.
- If it fails in the edited file, inspect whether the pass 029 assertion
  deletion caused the failure and fix only that.
- Do not broaden this into source-route work.

Do not:

- add tests;
- change source code unless a narrow test-harness blocker prevents validation;
- edit docs except the response/review log;
- touch PQS route behavior, final-basis H1, IDA, density-density, RHF, drivers,
  exports, or artifacts.

Validation:

- The slow nested harness command above is the validation.
- Also run `julia --project=. -e 'using GaussletBases; println("load ok")'`
  and `git diff --check` if any file is edited.

Deletion/shrinkage report required:

- if no files are changed, say so explicitly;
- if a narrow harness fix is made, say what old assertion/helper behavior was
  simplified or why no deletion was possible;
- report whether any new test was avoided;
- report remaining validation blockers, if any.

Report back:

- write `.agent_handoffs/response.031.md.tmp`, then atomically rename to
  `.agent_handoffs/response.031.md`;
- also write the curated copy to
  `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.031.md`;
- include validation result;
- include files changed, if any;
- include deletion/shrinkage report;
- sign `-- repo-doer@macmini`.

-- repo-manager@macmini
