Accepted.

Pass 020 validated the new CPBM final-basis object against the real
projected-q-shell shell-realization oracle:

```text
final_shell_coefficients == shell_plan.shell_isometry_matrix
final_overlap == shell_plan.realized_overlap
final identity error == shell_plan.isometry_error
```

The probe used `_pqs_shell_realization_plan(...)` and
`_pqs_product_box_support_overlap_matrix(...)` only as oracle/input sources. It
did not call `_pqs_current_route_safe_term_matrices(...)` and did not
materialize one-body operators, H1, IDA, density-density, RHF, drivers, exports,
or old fixed-block matrix authority.

Manager validation:

- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git diff --check`

Next blocker remains:

```text
:missing_pqs_shell_projected_one_body_operator_materialization
```

The next pass should add the algebraic shell-projected operator seam for a
caller-supplied shell-support operator matrix:

```text
O_boundary = P' * O_shell_support * P
O_final = R' * O_shell_support * R
```

It should not yet generate kinetic/nuclear shell operators or run H1.

-- repo-manager@macmini
