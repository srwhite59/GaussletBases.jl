# Algorithm Implementation Index

This is a navigation index for agents before Cartesian/PQS numerical work. It
does not create new algorithm authority. Use it to find existing implementations,
optimization lessons, and oracle/reference paths before writing new code.

## Term-First Coulomb Gaussian Contractions

Why to check:
Avoid rebuilding one dense 3D matrix per Coulomb Gaussian term. The repo already
records the faster shape: prepare one-dimensional factor data in term-first form
and keep the Coulomb expansion index as the short inner reduction.

Key docs:
- `docs/ordinary_cartesian_gg_term_first_optimization.md`
- `docs/ordinary_cartesian_gg_pairfactor_legacy_optimization.md`
- `docs/ordinary_cartesian_gg_scalar_legacy_optimization.md`
- `docs/ordinary_cartesian_qiu_white_ga_cross_legacy_optimization.md`

Source anchors:
- `src/ordinary_mapped_backends.jl`: `_mapped_coulomb_expanded_symmetric_matrix`
- `src/ordinary_cartesian_ida.jl`
- `src/ordinary_qw_raw_blocks.jl`
- `src/ordinary_qw_operator_assembly.jl`

Do-not-forget rule:
Use the Gaussian expansion term index as an inner reduction:

```julia
sum(c .* Fx[:, ix, jx] .* Fy[:, iy, jy] .* Fz[:, iz, jz])
```

The source usually implements this as an explicit loop to avoid broadcast
temporaries.

## Gaussian Factor And Pair-Factor Reuse

Why to check:
Avoid slow midpoint-grid or per-term reconstruction when mapped proxy/raw-layer
factor builders already exist.

Key docs:
- `docs/ordinary_cartesian_gg_pairfactor_legacy_optimization.md`
- `docs/ordinary_cartesian_gg_scalar_legacy_optimization.md`
- `docs/ordinary_cartesian_qiu_white_ga_cross_legacy_optimization.md`
- `docs/src/developer/pqs_manager_running_log.md`, Slice B Gaussian factor reuse
  entries

Source anchors:
- `src/ordinary_mapped_backends.jl`: `_mapped_legacy_proxy_layer`,
  `_mapped_legacy_proxy_scalar_data`, `mapped_ordinary_one_body_operators`
- `src/ordinary_coulomb.jl`: `coulomb_gaussian_expansion`,
  `gaussian_factor_matrices`
- `src/ordinary_cartesian_ida.jl`: `_pair_gaussian_factor_matrices`
- `src/ordinary_qw_raw_blocks.jl`: `_qwrg_cross_1d_blocks`

Do-not-forget rule:
Prefer analytic mapped proxy/raw-layer contraction and existing factor packets
over sampled midpoint-grid kernels or duplicated primitive formulas.

## CPB Local Blocks And Parent Axis Factors

Why to check:
Coordinate-product-box geometry, parent factor packets, and CPB local block
providers are intended to replace ad hoc route/report payload factor plumbing.

Key docs:
- `docs/src/developer/cartesian_coordinate_product_box_contract.md`
- `docs/src/developer/cartesian_parent_factors_and_cpb_kernels.md`

Source anchors:
- `src/CartesianParentAxisFactors.jl`: `parent_overlap_axis_factor_packet`
- `src/CartesianCPBBlockProviders.jl`: `cpb_interval_pair`,
  `cpb_overlap_axis_blocks`, `cpb_overlap_dense_block`,
  `cpb_electron_nuclear_by_center_local_block`

Do-not-forget rule:
CPB geometry stays pure. Parent factors and CPB providers own reusable local
operator blocks; route-global placement is a separate layer.

## PQS Terminal Basis Realization And Lowdin

Why to check:
Avoid confusing shell-row truncation, raw source operators, and final shell
realization. The terminal basis route uses previous-block projection and
shell-local Lowdin; no global Lowdin repair is allowed.

Key docs:
- `docs/src/developer/cartesian_hamiltonian_producer_design.md`
- `docs/src/developer/pqs_source_box_operator_framework.md`
- `docs/src/developer/pqs_near_term_final_basis_realization_plan.md`
- `JuliaStyle.md`, Lowdin guidance

Source anchors:
- `src/cartesian_final_basis_realization/pqs_terminal_basis_realization.jl`
- `src/cartesian_final_basis_realization/pqs_source_shell_final_basis.jl`
- `src/cartesian_final_basis_realization/pqs_complete_core_shell_final_basis.jl`
- `src/pqs_multilayer_complete_core_shell_h1.jl`

Do-not-forget rule:
Use `inv(sqrt(Symmetric(overlap)))` for symmetric Lowdin. Build from raw
source modes, then project and clean shell-locally. Do not substitute a global
cleanup.

## PQS Source-Box Retained Transforms

Why to check:
Boundary product-mode selection is the PQS retained rule. Shell realization is a
later step and should not be mixed into raw source-box retained-transform
authority.

Key docs:
- `docs/src/developer/raw_product_source_retained_transform_policy.md`
- `docs/src/developer/projected_q_shell_policy.md`
- `docs/src/developer/pqs_source_box_operator_framework.md`

Source anchors:
- `src/cartesian_route_core/retained_spaces.jl`
- `src/cartesian_contracted_parent_metrics/source_box_pair_shadow.jl`:
  `_pqs_raw_product_box_plan`, `_pqs_product_box_realization_plan`
- `src/cartesian_pair_block_materialization/pqs_source_safe_terms.jl`
- `src/cartesian_pair_block_materialization/pqs_source_one_body.jl`

Do-not-forget rule:
COMX boundary product modes define retained source modes. Support rows and
shell-local final coefficients are later realization details.

## One-Body Operators And Unit Nuclear Convention

Why to check:
Prevent sign, charge, and center-summing errors.

Key docs:
- `docs/src/developer/cartesian_hamiltonian_producer_design.md`
- `docs/src/developer/numerical_contracts.md`
- `docs/src/developer/pqs_near_term_final_basis_realization_plan.md`

Source anchors:
- `src/cartesian_final_basis_realization/pqs_terminal_one_body.jl`
- `src/ordinary_qw_operator_assembly.jl`: `assembled_one_body_hamiltonian`
- `src/pqs_multilayer_support_one_body.jl`:
  `pqs_multilayer_support_electron_nuclear_by_center_matrices`
- `src/cartesian_pair_block_materialization/one_body_global_electron_nuclear.jl`

Do-not-forget rule:
By-center nuclear matrices are uncharged unit attractions, `U_A = -1/r_A`.
Apply `Z_A` and sum centers only when forming the Hamiltonian.

## IDA And Density Pair-Factor Conventions

Why to check:
Avoid mixing raw/source weights, density-normalized pair factors, retained
weights, and final IDA weights.

Key docs:
- `docs/src/developer/numerical_contracts.md`
- `docs/src/developer/white_lindsey_low_order_density_density_builder_contract_2026-06-04.md`
- `docs/src/developer/raw_product_source_retained_transform_policy.md`
- `docs/src/developer/pqs_source_box_operator_framework.md`

Source anchors:
- `src/cartesian_nested_faces.jl`: `_nested_factorized_weight_aware_pair_terms`,
  `_nested_weight_aware_pair_terms`
- `src/ordinary_qw_raw_blocks.jl`: `_qwrg_fixed_block_interaction_matrix`
- `src/cartesian_contracted_parent_metrics/core.jl`:
  `_pqs_source_box_ida_factor_provenance`

Do-not-forget rule:
Carry raw pair numerators through projection/realization first. Normalize only
at the reviewed final retained/final-basis density boundary.

## Performance And Reuse Policy

Why to check:
Correct tiny tests can still hide unusable algorithms.

Key docs:
- `docs/src/developer/performance_review_contracts.md`
- `docs/src/developer/cartesian_route_retirement_ledger.md`
- `docs/code_bloat_and_wrong_contract_cleanup_note.md`
- `AGENTS.md`, performance review policy and anti-bloat gate

Do-not-forget rule:
Before nontrivial numerical coding, define scaling, identify reuse points,
measure representative time/allocation, and delete or quarantine stale paths
instead of preserving parallel implementations.

## High-Order Doside, Endcap, And Panel Donors

Why to check:
High-order notes contain useful donor algorithms and stress-test lessons, but
not all of them are current production authority.

Key docs:
- `docs/src/algorithms/cartesian_nested_endcap_panel_shared_shell.md`
- `docs/src/developer/high_order_endcap_panel_mainline_validation_2026-05-16.md`
- `docs/src/developer/high_order_doside_experimental_plan.md`
- `docs/src/developer/high_order_doside_physical_block_speed_plan.md`

Source anchors:
- `src/cartesian_nested_owned_units.jl`: `_nested_endcap_panel_owned_units`,
  `_nested_endcap_panel_shell_layer`
- `src/cartesian_high_order_doside_experimental.jl`
- `src/cartesian_high_order_doside_ida_experimental.jl`
- `src/cartesian_nested_faces.jl`: `_nested_doside_1d`

Do-not-forget rule:
Treat high-order doside/endcap/panel code as opt-in, experimental, or donor
material unless a current design explicitly promotes it.

## Migration And Oracle Paths

Why to check:
Old fixed-block, shell-row, route-shadow, and donor code can be valuable
oracles, but should not decide the current route structure.

Key docs:
- `docs/src/developer/cartesian_route_retirement_ledger.md`
- `docs/src/developer/cartesian/route_migration.md`
- `docs/src/developer/cartesian/feature_donor_inventory.md`
- `docs/src/developer/old_flat_cartesian_retirement_audit_2026-06-14.md`

Source anchors:
- `src/cartesian_nested_faces.jl`
- `src/ordinary_qw_raw_blocks.jl`
- `src/ordinary_qw_operator_assembly.jl`
- `src/cartesian_contracted_parent_metrics/`

Do-not-forget rule:
Old code may feed adapters or oracle comparisons. It should not silently become
new route authority.
