Accepted with one follow-up requirement.

Pass 019 added the intended CPBM-owned final-basis object:

```text
pqs_source_shell_realization_final_basis(...)
```

It computes:

```text
R = shell_projection * lowdin_cleanup
projected_boundary_overlap = P' * S_shell * P
final_overlap = R' * S_shell * R
```

and it leaves one-body materialization blocked with:

```text
:missing_pqs_shell_projected_one_body_operator_materialization
```

The implementation does not transform H1, does not use Lowdin alone on a
retained-source operator, and does not call current-route safe-term matrices.

Manager validation:

- `julia --project=. test/nested/cartesian_pair_block_materialization_contract_runtests.jl`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git diff --check`

Follow-up before one-body operators:

The first test uses an identity shell realization. That is acceptable as a
small algebraic contract check, but it is not enough to validate the real
projected-q-shell data path. The next pass should be a `tmp/work` probe that
feeds this helper from a real old `_pqs_shell_realization_plan(...)` fixture
and compares `R`, final overlap, and isometry error to the old shell plan as
oracle.

-- repo-manager@macmini
