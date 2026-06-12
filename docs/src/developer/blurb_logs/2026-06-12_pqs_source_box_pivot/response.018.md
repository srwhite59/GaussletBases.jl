Audit complete. I did not change source code or tests.

Surfaces read:

- `docs/src/developer/pqs_source_box_operator_framework.md`
- `docs/src/developer/raw_product_source_retained_transform_policy.md`
- `docs/src/developer/projected_q_shell_policy.md`
- `src/CartesianContractedParentMetrics.jl`
- `src/cartesian_nested_faces.jl`
- `src/cartesian_pair_block_materialization/pqs_source_shell_bridge.jl`
- `src/cartesian_pair_block_materialization/pqs_source_final_readiness.jl`
- `src/cartesian_pair_block_materialization/one_body_placement_plan.jl`
- focused PQS integration test sections in `test/nested/pqs_projected_q_shell_local_layer_integration_runtests.jl` and `test/nested/bond_aligned_diatomic_high_order_recipe_opt_in_source_construction_integration_runtests.jl`

Target Card

1. Minimal final-basis object for a one-center PQS H1 probe

The minimal object should be a compact PQS shell-realized final-basis object, not another source-box operator wrapper. It needs to carry:

- raw source identity and source mode dimensions/order;
- boundary COMX-product selector from full source product modes to selected boundary modes;
- shell support indices/states;
- shell projection matrix from selected boundary source modes to shell rows;
- full-rank symmetric Lowdin cleanup matrix;
- final shell coefficient matrix on shell rows, `shell_projection * lowdin_cleanup`;
- realized overlap/isometry diagnostics against the shell support metric;
- final retained dimension and column range;
- nonclaim flags that it has not built IDA, Hamiltonian/export/artifacts, or driver data.

The existing old helper `_pqs_shell_realization_plan(...)` already has most of this shape as an oracle object. The production-facing object should be owned in the CPBM/PQS source-final boundary, consume CRPS/CPBM source facts, and avoid depending on CCPM as route authority.

2. Transform direction

The direction is:

```text
raw source product modes
  -> boundary COMX-product selection
  -> shell-row projection
  -> symmetric Lowdin cleanup
  -> final shell-realized PQS columns
```

With:

- `B`: selector from full source mode space to boundary source modes;
- `P`: shell projection matrix, shell support rows x boundary source modes;
- `L`: Lowdin cleanup, boundary source modes x final retained columns;
- `R = P * L`: final shell-row coefficient matrix.

Lowdin alone is not the raw-to-final transform. `L` only cleans the shell-projected boundary-mode metric. Applying `L' * O_boundary * L` to the existing full raw-box boundary operator is not a valid final shell-realized operator unless the operator has already been projected to the shell-row metric/operator boundary.

3. One-body transform formula

For a shell-realized final operator, the safe formula is:

```text
O_shell_boundary = P' * O_shell_support * P
O_final = L' * O_shell_boundary * L
```

Equivalently:

```text
O_final = R' * O_shell_support * R
```

where `O_shell_support` is the operator restricted to the owned shell support rows. The implementation should build this through factorized/source-owned operator pieces where possible; it should not fall back to generic support-local shell-row contraction as the PQS algorithm.

If the route wants to reuse CPBM source-box operator blocks directly, the missing piece is an exact compact source-space realization transform for the shell-zeroed final functions. The existing `_pqs_current_route_shell_realization_transform_fact(...)` correctly reports that this compact transform is not available and that source-box operator application is not ready.

4. Required overlap/isometry checks

Before an ordinary final H1 solve:

- `P' * S_shell * P` must be finite symmetric positive definite for full-rank cleanup;
- Lowdin cleanup must retain the expected rank, with no silent dropped directions for the first one-center probe;
- `R' * S_shell * R` must be identity within the same final-basis tolerance used elsewhere;
- final overlap symmetry and identity error must be reported;
- only after those checks can the H1 solve be ordinary Hermitian in the final basis.

No generalized final-basis solve should be introduced for this path.

5. Oracle/kernel references versus route authority

Usable as oracle or kernel reference:

- `_nested_projected_q_shell_boundary_comx_product_modes(...)`
- `_nested_projected_q_shell_descriptor_seed_coefficients(...)`
- `_nested_projected_q_shell_cleanup(...)`
- `_pqs_shell_realization_plan(...)`
- selected descriptor/metric prototype checks in `cartesian_nested_faces.jl`

Compatibility wrapper only:

- `_pqs_product_box_realization_plan(...)`, because it returns both raw and shell plans for old diagnostics.

Do not adopt as production route authority:

- `_pqs_current_route_safe_term_matrices(...)`; it explicitly uses support-local fallback/oracle paths for shell-realized PQS pairs.
- `_pqs_current_route_shell_realization_transform_fact(...)` as a final transform; it is useful negative metadata and currently reports `source_box_operator_application_ready = false`.
- old current-route safe-term/authority comparison tests as algorithm definitions.

6. Smallest next implementation pass

Add a narrow CPBM-owned PQS shell-realization final-basis object for a single raw source/PQS unit:

```text
pqs_source_shell_realization_final_basis(...)
```

or similarly named internal helper. It should consume:

- the raw source fact / source mode metadata already used by the pass-017 probe;
- the boundary retained rule;
- shell support indices/states;
- source-axis transform facts;
- the shell projection and Lowdin cleanup data.

The first implementation should only materialize the final basis and its overlap/isometry diagnostics. It should not yet transform `T`/`V` or solve H1 unless the exact shell-projected operator path is present in the same data.

If that object is available, the following pass can add final one-body matrix materialization:

```text
O_final = L' * (P' * O_shell_support * P) * L
```

for overlap, kinetic, and one uncharged electron-nuclear by-center term, preserving charge application at the H1 assembly stage.

7. Next blocker after this audit

The current blocker should be sharpened from:

```text
:missing_pqs_shell_realization_lowdin_final_basis_construction
```

to:

```text
:missing_pqs_shell_realization_final_basis_object
```

and then, after that object exists:

```text
:missing_pqs_shell_projected_one_body_operator_materialization
```

The second blocker is important because Lowdin cleanup alone must not be treated as a valid transform of the retained-source H1 block.

Deletion/shrinkage report:

- Nothing was deleted in this audit-only pass.
- No existing code surface became obsolete yet because there is still no production CPBM final-basis object for PQS shell realization.
- After the next implementation, the metadata-only CPBM bridge/readiness tests that only assert `:shell_realization_not_materialized` should shrink to one blocker smoke, and pass-015/pass-017 source-box H1 probes should be demoted to historical artifacts.
- Once a final PQS H1 path exists, old current-route safe-term tests that assert support-local fallback/oracle vocabulary should be reduced or quarantined. The likely first shrink targets are sections of `test/nested/pqs_projected_q_shell_local_layer_integration_runtests.jl` around current-route/component-route smoke and shell-realization metadata flags.
- The high-order source-construction integration test should keep compact checks for descriptor construction, boundary-mode selection, and Lowdin isometry, but it should not remain responsible for proving source-box operator readiness.

Recommended next blurb:

Implement the first CPBM-owned PQS shell-realization final-basis object for one source box. It should consume repo-owned source-axis transform facts, the PQS boundary retained rule, shell support rows, shell projection, and Lowdin cleanup; report final overlap/isometry diagnostics; and keep operator materialization blocked until the shell-projected one-body operator formula is implemented. Do not transform H1 with Lowdin alone, do not use current-route safe-term matrices as authority, and do not add IDA/RHF/global driver work.

Validation:

- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git diff --check`

-- repo-doer@macmini
