Accepted.

Pass 021 added the intended CPBM algebraic seam:

```text
pqs_source_shell_projected_one_body_matrix(final_basis, shell_operator; term)
```

It validates a caller-supplied shell-support operator, computes:

```text
O_boundary = P' * O_shell_support * P
O_final = R' * O_shell_support * R
```

and cross-checks the equivalent Lowdin form:

```text
O_final = L' * O_boundary * L
```

The implementation did not generate kinetic/nuclear operators, did not run H1,
and did not call `_pqs_current_route_safe_term_matrices(...)`. The added test is
compact live-contract coverage for the algebraic projection seam.

Manager validation:

- `julia --project=. test/nested/cartesian_pair_block_materialization_contract_runtests.jl`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git diff --check`

Minor note:

The result still carries several explicit nonclaim booleans. That is acceptable
for this short-lived transition seam because those flags guard against the
specific drift risk: treating a projected one-body matrix as H1, IDA, RHF, or a
driver route. Once a real final PQS H1 probe exists, these report fields should
be candidates for slimming.

Next blocker:

```text
:missing_pqs_shell_support_one_body_operator_source
```

Before implementing a source, the next pass should audit whether the existing
PQS retained-source overlap/kinetic blocks are already equal to the
shell-projected boundary operators `P' O_shell P` on the real projected-q-shell
fixture. If they are, the next implementation can transform proven boundary
operators with `L' O_boundary L`. If they are not, a true shell-support operator
source is needed.

-- repo-manager@macmini
