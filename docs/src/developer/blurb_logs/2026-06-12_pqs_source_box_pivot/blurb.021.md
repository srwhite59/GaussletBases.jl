Purpose:

Add the first CPBM-owned shell-projected final-basis one-body operator helper
for PQS. This should clear the algebraic part of:

```text
:missing_pqs_shell_projected_one_body_operator_materialization
```

without yet generating kinetic, nuclear, H1, IDA, density-density, RHF, driver,
export, or artifact behavior.

Why now:

Pass 020 proved that `pqs_source_shell_realization_final_basis(...)` matches the
old projected-q-shell shell-realization oracle on a real fixture. The next live
seam is no longer the final-basis object itself; it is projecting a supplied
shell-support operator into that final orthonormal basis.

Exact task:

Implement a narrow helper in CPBM, likely near
`src/cartesian_pair_block_materialization/pqs_source_shell_final_basis.jl`:

```julia
pqs_source_shell_projected_one_body_matrix(final_basis, shell_operator; term)
```

or an equally small API if the existing file shape suggests a better name.

The helper should:

- require `final_basis.status == :available_pqs_shell_realization_final_basis`;
- accept a caller-supplied real `shell_operator` matrix in shell-support rows;
- validate shape against `final_basis.shell_support_count`;
- validate finite entries;
- for symmetric real one-body terms, validate symmetry and report the symmetry
  error;
- compute:

```text
O_boundary = P' * O_shell_support * P
O_final = R' * O_shell_support * R
```

where `P = final_basis.shell_projection` and
`R = final_basis.final_shell_coefficients`;

- optionally cross-check `O_final` against `L' * O_boundary * L`, where
  `L = final_basis.lowdin_cleanup`, and report the difference;
- return a compact result with:
  - status;
  - blocker;
  - term;
  - boundary operator;
  - final operator;
  - shape/rank/symmetry/finite diagnostics;
  - explicit nonclaims that no H1 solve, charge summing, IDA, density-density,
    RHF, driver route, export, or artifact was materialized.

Trust boundary:

This is an operator projection seam only. Do not call
`_pqs_current_route_safe_term_matrices(...)`. Do not generate analytic
kinetic/nuclear shell-support operators in this pass. Do not use old fixed-block
matrices as route authority. The old shell plan may be used only in a `tmp/work`
probe as oracle input for the final-basis object, as in pass 020.

Test policy:

Add at most one compact module-contract test if needed. It should validate the
algebraic projection seam using a small synthetic shell operator and should not
be a broad metadata test.

Also run a `tmp/work` real-fixture probe if useful:

- use the pass-020 projected-q-shell fixture;
- use `shell_overlap` as a supplied operator;
- confirm the helper returns a final operator equal to the final overlap already
  validated by `pqs_source_shell_realization_final_basis(...)`.

Do not add this real-fixture probe as a permanent test unless a non-obvious bug
requires it.

Expected blocker after this pass:

If the algebraic projection helper is available, the remaining blocker should
advance from:

```text
:missing_pqs_shell_projected_one_body_operator_materialization
```

to the next constructive source blocker, likely:

```text
:missing_pqs_shell_support_one_body_operator_source
```

or the exact name the implementation justifies. Do not pretend H1 is ready.

Deletion/shrinkage report required:

- what old code, test, metadata, or compatibility path became unnecessary;
- what was deleted or simplified;
- if nothing was deleted, why no existing surface was made obsolete yet;
- whether any new test/artifact was added and why it earned its cost;
- any remaining stale or duplicate surfaces to retire next.

Validation:

- focused module-contract test if one is added or edited;
- `julia --project=. -e 'using GaussletBases; println("load ok")'`;
- `git diff --check`.

Report back:

- write `.agent_handoffs/response.021.md.tmp`, then atomically rename to
  `.agent_handoffs/response.021.md`;
- also write the curated copy to
  `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.021.md`;
- include files changed;
- include any probe artifact path;
- include exact blocker after this pass;
- include validation run;
- include deletion/shrinkage report;
- sign `-- repo-doer@macmini`.

-- repo-manager@macmini
