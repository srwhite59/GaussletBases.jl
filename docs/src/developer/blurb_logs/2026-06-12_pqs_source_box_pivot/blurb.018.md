Purpose:

Audit the PQS shell-realization / Lowdin final-basis boundary before
implementation.

Why now:

The retained-source H1 path now has repo-owned source-axis transforms and
source-box one-body matrices. The current blocker is:

```text
:missing_pqs_shell_realization_lowdin_final_basis_construction
```

This is a design boundary. Do not start coding until the route is stated
precisely enough to avoid turning shell-row support contraction into the PQS
algorithm.

Exact task:

Read the relevant docs and code, then write an audit/target-card response. No
source code changes and no tests.

Required surfaces to inspect:

- `docs/src/developer/pqs_source_box_operator_framework.md`
- `docs/src/developer/raw_product_source_retained_transform_policy.md`
- `docs/src/developer/projected_q_shell_policy.md`
- `src/CartesianContractedParentMetrics.jl`
  - `_pqs_shell_realization_plan(...)`
  - `_pqs_product_box_realization_plan(...)`
  - `_pqs_current_route_shell_realization_transform_fact(...)`
  - `_pqs_current_route_safe_term_matrices(...)`
- `src/cartesian_nested_faces.jl`
  - projected q-shell / symmetric Lowdin helpers around the current shell
    realization code
- focused tests that exercise these surfaces, especially
  `test/nested/pqs_projected_q_shell_local_layer_integration_runtests.jl` and
  `test/nested/bond_aligned_diatomic_high_order_recipe_opt_in_source_construction_integration_runtests.jl`

Questions to answer:

1. What is the minimal final-basis object needed for a one-center PQS H1 probe?
2. What is the exact transform direction from retained source modes to final
   shell-realized/Lowdin-cleaned columns?
3. Which matrices should transform one-body operators, and by what formula?
4. What overlap identity/isometry checks are required before an ordinary final
   solve?
5. Which old helper can be used as oracle/kernel reference, and which old
   helper must not be adopted as route authority?
6. What is the smallest next implementation pass?
7. What old code/test surface would become less necessary after that pass?

Important policy:

- Shell projection plus Lowdin belongs to final realization, not raw-box
  operator construction.
- Do not use Lowdin cleanup alone as the full raw-to-final transform.
- Do not treat support-local shell-row contraction as the PQS algorithm.
- Old current-route safe-term matrices are oracle/debug, not production
  authority.
- Do not add metadata fields, helpers, tests, or docs in this audit pass.

Deliverable:

Write the audit in:

```text
.agent_handoffs/response.018.md
docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.018.md
```

No tracked source/docs/test changes except the curated response file.

Deletion/shrinkage report required:

- what old code, test, metadata, or compatibility path would become
  unnecessary after the proposed next implementation;
- if nothing can be retired yet, explain why;
- identify any existing tests likely to shrink once final-basis PQS H1 exists.

Validation:

- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git diff --check`

Report back:

- include exact next implementation blurb recommendation;
- include validation run;
- include deletion/shrinkage report;
- sign `-- repo-doer@macmini`.

-- repo-manager@macmini
