Implemented the first CPBM-owned PQS shell-realization final-basis object.

Implementation:

- Added `src/cartesian_pair_block_materialization/pqs_source_shell_final_basis.jl`.
- Exported and included `pqs_source_shell_realization_final_basis(...)` from `CartesianPairBlockMaterialization`.
- The helper consumes a `CartesianRawProductSources.RawProductBoxPlan` plus a `PQSBoundaryProductModeRetainedRule`, validates that the rule is the raw PQS boundary source-mode selector, then computes:

```text
R = shell_projection * lowdin_cleanup
projected_boundary_overlap = shell_projection' * shell_overlap * shell_projection
final_overlap = R' * shell_overlap * R
```

- It reports compact rank/eigenvalue/symmetry/identity diagnostics for the projected-boundary and final overlaps.
- It returns blocked status instead of accepting a reduced or non-identity final basis:
  - `:shell_realization_rank_deficient`
  - `:final_overlap_not_identity`
- It explicitly keeps one-body/H1/IDA/RHF/driver/export/artifact claims false and records:

```text
one_body_operator_blocker = :missing_pqs_shell_projected_one_body_operator_materialization
```

Test coverage:

- Added one compact CPBM contract check in `test/nested/cartesian_pair_block_materialization_contract_runtests.jl`.
- The test uses a CRPS `(5, 5, 5)` PQS boundary retained rule and a minimal identity shell realization to verify:
  - final-basis shape `98 x 98`;
  - `R = P * L`;
  - final overlap identity;
  - no one-body/H1/IDA/RHF/driver/artifact claims;
  - one-body materialization remains blocked by `:missing_pqs_shell_projected_one_body_operator_materialization`.

No old route authority was adopted:

- No call to `_pqs_current_route_safe_term_matrices(...)`.
- No support-local safe-term operator contraction.
- No H1 solve.
- No Lowdin-only transform of retained-source H1.
- No IDA, density-density, RHF, global driver, exports, artifacts, direct Cartesian fallback, full-parent fallback, or old fixed-block matrix authority.

Validation:

- `julia --project=. test/nested/cartesian_pair_block_materialization_contract_runtests.jl`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git diff --check`

Next blocker:

```text
:missing_pqs_shell_projected_one_body_operator_materialization
```

The next implementation should materialize final shell-projected one-body operators with:

```text
O_final = L' * (P' * O_shell_support * P) * L
```

or the equivalent `R' * O_shell_support * R`, without treating the existing retained-source block as transformable by Lowdin alone.

Deletion/shrinkage report:

- No old code was deleted in this pass. The new final-basis object is a first live production seam, but source-box bridge/readiness and old oracle surfaces are still needed until one-body materialization exists.
- No compatibility adapter was added to preserve old tuple/metadata shapes.
- The new test is genuinely new live-contract coverage for the CPBM final-basis seam; it does not duplicate a broad old integration path.
- Once shell-projected one-body materialization exists, CPBM bridge/readiness tests that only assert `:shell_realization_not_materialized` should shrink to one blocker smoke.
- Remaining stale surfaces to retire next: pass-015/pass-017 retained-source H1 probes should become historical/oracle artifacts after final shell-projected H1 exists; current-route safe-term tests that preserve support-local PQS oracle vocabulary should then be reduced or quarantined.

-- repo-doer@macmini
