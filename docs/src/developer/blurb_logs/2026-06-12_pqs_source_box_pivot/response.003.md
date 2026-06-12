Implemented the first retained PQS safe one-body source-mode contraction from existing raw source-space PQS/PQS blocks.

Files changed:

- `src/cartesian_pair_block_materialization/CartesianPairBlockMaterialization.jl`
- `src/cartesian_pair_block_materialization/preflight.jl`
- `src/cartesian_pair_block_materialization/pqs_source_safe_terms.jl`
- `test/nested/cartesian_pair_block_materialization_contract_runtests.jl`

What changed:

- Added `pqs_source_pair_retained_one_body_block(source_result, left_rule, right_rule)`.
- Added `pqs_source_pair_retained_one_body_block(source_result)` using retained rules carried from preflight metadata.
- Added convenience helpers:
  - `pqs_source_pair_retained_overlap_block(record; overlap_1d)`
  - `pqs_source_pair_retained_kinetic_block(record; overlap_1d, kinetic_1d)`
- PQS/PQS source-pair preflight now requires left/right `PQSBoundaryProductModeRetainedRule` objects from transform-contract metadata.
- PQS/PQS source-pair preflight now carries the retained rules and compact retained-rule summaries into pair-block materialization record metadata.
- Raw PQS source-space block results now preserve left/right retained rules and side-specific source-mode ordering metadata so retained contraction can validate the selector source.
- Fixed an active preflight seam for `CartesianUnitPairs.UnitPairIndexTable`: `pair_block_materialization_plan` now walks pair-operator records and unit pairs together instead of indexing the lightweight pair table as if it were a rich tuple/vector.

Retained block evidence:

- Fixture source-space overlap block shape: `(27, 60)`.
- Left retained selector count: `26`.
- Right retained selector count: `54`.
- Retained overlap block shape: `(26, 54)`.
- Retained overlap block equals `source_overlap[left_columns, right_columns]`.
- Retained kinetic block shape: `(26, 54)`.
- Retained kinetic block equals `source_kinetic[left_columns, right_columns]`.
- Retained results report `block_space = :retained_pqs_source_modes`.
- Retained results report source-space input was used.
- Retained results keep `shell_realization_materialized = false` and `lowdin_cleanup_used = false`.
- No final shell-realized pair block, Hamiltonian data, IDA data, export, or artifact path is claimed.

Validation run:

- `julia --project=. test/nested/cartesian_pair_block_materialization_contract_runtests.jl` passed.
- `julia --project=. test/nested/cartesian_retained_unit_transform_contracts_runtests.jl` passed.
- `julia --project=. -e 'using GaussletBases; println("load ok")'` passed with `load ok`.
- `git diff --check` passed.

Deletion/shrinkage report:

- No old code was deleted in this pass.
- The active preflight path was simplified to avoid assuming materialized/indexable rich unit-pair storage; it now consumes the lightweight pair iterator directly. That bypasses the stale tuple/vector indexing assumption for the active path.
- No old PQS source-box helper became obsolete yet because this pass adds the first retained source-mode contraction but does not replace shell-realization or CCPM oracle paths.
- The new test coverage is genuinely new live-contract coverage for retained source-mode block contraction. It does not replace an older test.
- Remaining stale/duplicate surfaces to retire next:
  - old shell/support-row source-box reference helpers should remain oracle-only until retained source-mode pair-block contraction is wired into the next route layer;
  - source-space-only bridge/readiness wording can be narrowed once retained source-mode blocks become the active consumer for the next source-box-first PQS step.

-- repo-doer@macmini
