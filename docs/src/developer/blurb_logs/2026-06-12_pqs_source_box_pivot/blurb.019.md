Purpose:

Implement the first CPBM-owned PQS shell-realization final-basis object.

Why now:

Pass 018 identified the next blocker:

```text
:missing_pqs_shell_realization_final_basis_object
```

The goal of this pass is only to materialize and validate the final basis
itself:

```text
R = shell_projection * lowdin_cleanup
S_final = R' * S_shell_support * R
```

Do not transform one-body operators in this pass.

Exact task:

Add a narrow helper in CPBM, or block precisely if the existing data cannot be
provided without adopting old route authority.

Recommended API shape:

```text
pqs_source_shell_realization_final_basis(
    raw_source_plan,
    retained_rule;
    shell_support_indices,
    shell_overlap,
    shell_projection,
    lowdin_cleanup,
)
```

or a better local equivalent.

Required behavior:

- validate that `retained_rule` is the PQS boundary source-mode rule;
- validate shape consistency:
  - `shell_projection`: shell support rows x boundary source modes;
  - `lowdin_cleanup`: boundary source modes x final retained columns;
  - `shell_overlap`: shell support rows x shell support rows;
- compute `final_shell_coefficients = shell_projection * lowdin_cleanup`;
- compute `projected_boundary_overlap = shell_projection' * shell_overlap * shell_projection`;
- compute `final_overlap = final_shell_coefficients' * shell_overlap * final_shell_coefficients`;
- report rank/eigenvalue/identity diagnostics;
- require full rank for this first one-center fixture; if rank drops, block
  precisely instead of silently accepting a reduced basis;
- record nonclaims: no one-body operators, no H1, no IDA, no RHF, no driver,
  no exports/artifacts.

Oracle/reference:

Use `_pqs_shell_realization_plan(...)` or focused projected-q-shell fixtures
only as an oracle/reference for shape and isometry. Do not call
`_pqs_current_route_safe_term_matrices(...)` and do not treat support-local
operator contraction as the algorithm.

Test policy:

Add one compact module-contract test only if it uses an existing focused
projected-q-shell fixture or a minimal old-helper oracle object. Avoid a broad
metadata-field test. The test should check:

- `R = P * L` shape;
- final overlap identity error is small;
- no one-body/H1/IDA/RHF claims;
- operator materialization remains blocked by
  `:missing_pqs_shell_projected_one_body_operator_materialization`.

If the fixture setup is too entangled, do not invent a large synthetic test;
return the exact blocker and leave code unchanged.

Trust boundary:

No H1 solve, no one-body operator transform, no support-local safe-term matrix
authority, no IDA, no density-density, no RHF, no global driver, no exports,
no artifacts, no direct Cartesian fallback, no full-parent fallback, no old
fixed-block matrix authority.

Deletion/shrinkage report required:

- what old code, test, metadata, or compatibility path became unnecessary;
- what was deleted or simplified;
- if nothing was deleted, why no existing surface was made obsolete yet;
- whether any new test/artifact was added and why it earned its cost;
- any remaining stale or duplicate surfaces to retire next.

Validation:

- if production code/tests changed:
  `julia --project=. test/nested/cartesian_pair_block_materialization_contract_runtests.jl`
  or the focused existing test file touched by the change;
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git diff --check`

Report back:

- write `.agent_handoffs/response.019.md.tmp`, then atomically rename to
  `.agent_handoffs/response.019.md`;
- also write the curated copy to
  `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.019.md`;
- include implementation or exact blocker;
- include validation run;
- include deletion/shrinkage report;
- sign `-- repo-doer@macmini`.

-- repo-manager@macmini
