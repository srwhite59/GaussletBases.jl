Purpose:

Validate the new CPBM PQS final-basis object on a real projected-q-shell
fixture before adding one-body operators.

Why now:

Pass 019 implemented the final-basis object, but its first module-contract test
uses an identity shell realization. Before implementing
`O_final = R' * O_shell_support * R`, validate the object against the existing
old projected-q-shell shell-realization oracle.

Exact task:

Create a `tmp/work` probe. Do not change production code or tests unless a
small missing seam is exposed.

Use the existing focused surfaces:

- `src/CartesianContractedParentMetrics.jl`
  - `_pqs_shell_realization_plan(...)`
  - `_pqs_product_box_support_overlap_matrix(...)`
  - `_pqs_raw_product_box_plan(...)` if needed for source metadata only
- `test/nested/pqs_projected_q_shell_local_layer_integration_runtests.jl`
  for fixture setup patterns.

The probe should:

- build a real projected-q-shell descriptor and metrics fixture;
- obtain the old shell plan from `_pqs_shell_realization_plan(...)`;
- obtain the shell support overlap matrix, using the old helper only as
  oracle/input source for this probe;
- build a CRPS raw source plan and PQS boundary retained rule with matching
  source dims/key/order;
- call `pqs_source_shell_realization_final_basis(...)` with:
  - `shell_support_indices`;
  - `shell_overlap`;
  - `shell_projection = shell_plan.shell_projection_matrix`;
  - `lowdin_cleanup = shell_plan.lowdin_cleanup`;
- compare:
  - `final_shell_coefficients` to `shell_plan.shell_isometry_matrix`;
  - `final_overlap` to `shell_plan.realized_overlap`;
  - final identity error to `shell_plan.isometry_error`;
- confirm the one-body blocker remains
  `:missing_pqs_shell_projected_one_body_operator_materialization`.

Trust boundary:

The old shell plan is oracle/input source only. Do not call
`_pqs_current_route_safe_term_matrices(...)`. Do not materialize one-body
operators, H1, IDA, density-density, RHF, drivers, exports, or artifacts.

Test policy:

No new test in this pass. Use `tmp/work` artifacts and report paths.

Deletion/shrinkage report required:

- what old code, test, metadata, or compatibility path became unnecessary;
- what was deleted or simplified;
- if nothing was deleted, why no existing surface was made obsolete yet;
- whether any new test/artifact was added and why it earned its cost;
- any remaining stale or duplicate surfaces to retire next.

Validation:

- `julia --project=. <tmp/work probe>`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git diff --check`

Report back:

- write `.agent_handoffs/response.020.md.tmp`, then atomically rename to
  `.agent_handoffs/response.020.md`;
- also write the curated copy to
  `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.020.md`;
- include probe artifact path and key diagnostics;
- include exact next blocker;
- include validation run;
- include deletion/shrinkage report;
- sign `-- repo-doer@macmini`.

-- repo-manager@macmini
