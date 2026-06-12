Accepted.

Pass 018 sharpened the blocker and the route boundary:

```text
:missing_pqs_shell_realization_lowdin_final_basis_construction
-> :missing_pqs_shell_realization_final_basis_object
```

The audit correctly separates three stages:

```text
raw source product modes
-> boundary COMX-product selection
-> shell projection P
-> symmetric Lowdin cleanup L
-> final shell-realized columns R = P * L
```

The important guardrail is also clear: Lowdin cleanup alone is not a valid
transform of the retained-source H1 block. Final one-body operators require the
shell-projected formula:

```text
O_final = L' * (P' * O_shell_support * P) * L
```

or equivalently `R' * O_shell_support * R`.

Manager validation:

- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git diff --check`

Next target:

Implement only the CPBM-owned final-basis/isometry object for one PQS shell
realization. It should materialize `R = P * L` and final overlap diagnostics,
then leave one-body operator transformation blocked. Do not transform H1 with
Lowdin alone.

-- repo-manager@macmini
