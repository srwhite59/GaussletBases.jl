Accepted.

Pass 022 answered the key boundary question with a real projected-q-shell
fixture:

```text
retained_source_overlap == P' * S_shell_support * P
retained_source_kinetic == P' * T_shell_support * P
```

to roundoff. The probe also confirmed:

```text
L' * (P' * S_shell_support * P) * L == final_basis.final_overlap
```

and validated the kinetic shell-support projection through the pass-021 helper.

Manager validation:

- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git diff --check`
- reviewed `tmp/work/pqs_shell_support_operator_boundary_audit_probe_output.txt`

Important interpretation:

Lowdin-transforming retained-source overlap/kinetic is valid for this fixture
because the retained-source blocks have now been shown to equal the
shell-projected boundary operators `P' O_shell P`. It is not a generic license
to apply Lowdin cleanup to arbitrary raw source operators.

Next implementation target:

Add a CPBM final-from-boundary one-body helper for overlap/kinetic:

```text
O_final = L' * O_boundary * L
```

where `O_boundary` is a proven retained boundary operator, preferably supplied
as the existing `pqs_retained_source_one_body_matrix(...)` result rather than a
loose matrix. H1 should remain blocked until the electron-nuclear boundary
operator exists.

-- repo-manager@macmini
