# Algorithm Implementation Index

This is a navigation index for agents before Cartesian/PQS numerical work. It
does not create new algorithm authority. Use it to find existing implementations,
optimization lessons, and oracle/reference paths before writing new code.

Status labels:
- **active reusable kernel**: current implementation code that may be called
  when its contract matches the new route.
- **active implementation surface**: current code named by approved design
  authority; use it through its owning boundary rather than treating it as new
  public API authority.
- **active donor pattern**: current code worth copying organizationally, but
  not necessarily callable from the new route.
- **consumer example only**: useful for seeing how a kernel is used; do not copy
  its orchestration contract.
- **oracle/reference only**: comparison or migration code, not production route
  authority.
- **retired/do not call**: historical or explicitly disallowed path.

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
- **active donor pattern**: `src/ordinary_mapped_backends.jl`,
  `_mapped_coulomb_expanded_symmetric_matrix`
- **consumer example only**: `src/ordinary_cartesian_ida.jl`,
  `_ordinary_cartesian_ida_from_gausslet_bundle`
- **retired/do not call**: `src/ordinary_mapped_backends.jl`,
  `_mapped_cartesian_one_body_matrix` term-by-term `kron` path
- **active donor pattern**: `src/ordinary_qw_raw_blocks.jl`
- **consumer example only**: `src/ordinary_qw_operator_assembly.jl`

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
- **active reusable kernel**: `src/ordinary_coulomb.jl`,
  `coulomb_gaussian_expansion`, `gaussian_factor_matrices`
- **active reusable local data**: PGDG intermediates such as
  `pgdg.gaussian_factor_terms` when already centered for the requested axis and
  center
- **active donor pattern**: `src/ordinary_mapped_backends.jl`:
  `_mapped_legacy_proxy_layer`,
  `_mapped_legacy_proxy_scalar_data`, `mapped_ordinary_one_body_operators`
- **active donor pattern**: `src/ordinary_cartesian_ida.jl`,
  `_pair_gaussian_factor_matrices`
- **active donor pattern**: `src/ordinary_qw_raw_blocks.jl`,
  `_qwrg_cross_1d_blocks`
- **retired/do not call**: `src/ordinary_qw_raw_blocks.jl`,
  `_qwrg_cross_1d_blocks_midpoint`

Do-not-forget rule:
Prefer analytic mapped proxy/raw-layer contraction and existing factor packets
over sampled midpoint-grid kernels or duplicated primitive formulas.

## Cartesian Gaussian Nuclear Raw Blocks

Why to check:
Residual Gaussian augmentation and Qiu-White route code both need exact
by-center Cartesian Gaussian nuclear `G-A` and `A-A` matrices. The neutral
owner prevents duplicate route-local nuclear loops while keeping raw analytic
formula ownership out of Residual Gaussian physics and Qiu-White route objects.

Key docs:
- `docs/src/developer/designs/cartesian_hamiltonian_producer/cartesian_gaussian_raw_blocks_nuclear.md`
- `docs/src/developer/designs/cartesian_hamiltonian_producer/registry.md`,
  `HP-CGRB-FILE-01`, `HP-CGRB-FN-01`, `HP-CGRB-FN-02`,
  `HP-CGAI-FN-01`, `HP-CGRB-WIRE-01`, `HP-CGRB-TEST-01`

Source anchors:
- **active implementation surface after extraction**:
  `src/cartesian_gaussian_raw_blocks/CartesianGaussianRawBlocks.jl`,
  `src/cartesian_gaussian_raw_blocks/nuclear_blocks.jl`
- **current extraction targets / parity consumers**:
  `src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl`,
  `src/ordinary_qw_raw_blocks.jl`,
  `src/ordinary_qw_operator_assembly.jl`

Do-not-forget rule:
The neutral raw-block owner returns uncharged by-center `U_A = -1/r_A`
parent-supplement and supplement-supplement matrices. It does not apply
physical nuclear charges, project through terminal blocks, transform into
Residual Gaussian bases, construct overlap/kinetic/moment blocks, create route
objects, or own artifacts/status/payloads.

## Cartesian Gaussian Non-Nuclear Raw Blocks

Why to check:
After nuclear axis-family reuse, Cr2 q4 exact-operator cost moved to
non-nuclear mixed/self Gaussian raw blocks and residual mixed-overlap setup.
The neutral owner now has a separate non-nuclear slice for exact overlap,
kinetic, coordinate-moment, and second-moment `G-A`/`A-A` matrices.

Key docs:
- `docs/src/developer/designs/cartesian_hamiltonian_producer/cartesian_gaussian_raw_blocks_non_nuclear.md`
- `docs/src/developer/designs/cartesian_hamiltonian_producer/registry.md`,
  `HP-CGRB-NN-FILE-01`, `HP-CGRB-NN-FN-01`,
  `HP-CGRB-NN-WIRE-01`, `HP-CGRB-NN-TEST-01`

Source anchors:
- **active implementation surface after extraction**:
  `src/cartesian_gaussian_raw_blocks/non_nuclear_blocks.jl`
- **current extraction targets / parity consumers**:
  `src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl`,
  `src/ordinary_qw_raw_blocks.jl`,
  `src/ordinary_qw_operator_assembly.jl`

Do-not-forget rule:
This owner is for non-nuclear `G-A`/`A-A` raw blocks only: overlap, kinetic,
`x`/`y`/`z`, and `x^2`/`y^2`/`z^2`. It does not optimize final-basis `G-G`
product matrices, project through terminal blocks, transform into Residual
Gaussian bases, change Qiu-White semantics, create route objects, or own
artifacts/status/payloads.

## CPB Local Blocks And Parent Axis Factors

Why to check:
Coordinate-product-box geometry, parent factor packets, and CPB local block
providers are intended to replace ad hoc route/report payload factor plumbing.

Key docs:
- `docs/src/developer/cartesian_coordinate_product_box_contract.md`
- `docs/src/developer/cartesian_parent_factors_and_cpb_kernels.md`

Source anchors:
- **active reusable kernel**: `src/CartesianParentAxisFactors.jl`,
  `parent_overlap_axis_factor_packet`
- **active reusable kernel**: `src/CartesianCPBBlockProviders.jl`,
  `cpb_interval_pair`,
  `cpb_overlap_axis_blocks`, `cpb_overlap_dense_block`,
  `cpb_electron_nuclear_by_center_local_block`

Do-not-forget rule:
CPB geometry stays pure. Parent factors and CPB providers own reusable local
operator blocks; route-global placement is a separate layer.

## PQS Terminal Basis Realization And Lowdin

Why to check:
Avoid confusing shell-row truncation, raw source operators, and final shell
realization. The terminal basis route uses full-box boundary source modes
restricted to shell-owned support followed by shell-local Lowdin; previous-block
projection and global Lowdin repair are forbidden. Cross-block overlap is
structurally zero for disjoint owned supports; cross-block kinetic, nuclear,
and IDA operator terms may still be nonzero and are assembled over terminal
block pairs.

Key docs:
- `docs/src/developer/designs/cartesian_hamiltonian_producer/`
- `docs/src/developer/pqs_source_box_operator_framework.md`
- `docs/src/developer/pqs_near_term_final_basis_realization_plan.md`
- `JuliaStyle.md`, Lowdin guidance

Source anchors:
- **active reusable kernel**:
  `src/cartesian_final_basis_realization/pqs_terminal_basis_realization.jl`
- **active donor pattern**:
  `src/cartesian_final_basis_realization/pqs_source_shell_final_basis.jl`
- **oracle/reference only**:
  `src/cartesian_final_basis_realization/pqs_complete_core_shell_final_basis.jl`
- **oracle/reference only**: `src/pqs_multilayer_complete_core_shell_h1.jl`

Do-not-forget rule:
Use `inv(sqrt(Symmetric(overlap)))` for symmetric Lowdin. Build from raw
full source-box boundary modes, restrict rows to
`support.support_indices` / `support.support_states`, and clean only within the
shell-owned support. Do not project against previous terminal blocks, grow
support onto previous regions, or substitute a global cleanup.

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
- **active reusable kernel**: `src/cartesian_route_core/retained_spaces.jl`
- **active donor pattern**:
  `src/cartesian_contracted_parent_metrics/source_box_pair_shadow.jl`:
  `_pqs_raw_product_box_plan`, `_pqs_product_box_realization_plan`
- **oracle/reference only**:
  `src/cartesian_pair_block_materialization/pqs_source_safe_terms.jl`
- **retired/do not call**:
  `src/cartesian_pair_block_materialization/pqs_source_one_body.jl`

Do-not-forget rule:
COMX boundary product modes define retained source modes. Support rows and
shell-local final coefficients are later realization details.

## One-Body Operators And Unit Nuclear Convention

Why to check:
Prevent sign, charge, and center-summing errors.

Key docs:
- `docs/src/developer/designs/cartesian_hamiltonian_producer/`
- `docs/src/developer/numerical_contracts.md`
- `docs/src/developer/pqs_near_term_final_basis_realization_plan.md`

Source anchors:
- **active implementation surface**:
  `src/cartesian_final_basis_realization/pqs_terminal_one_body.jl`,
  `assemble_terminal_product_operator!`
- **active implementation surface under `HP-R3UN-FN-01`**:
  `src/cartesian_final_basis_realization/pqs_terminal_one_body.jl`,
  `_accumulate_terminal_gaussian_sum!`, `_terminal_gaussian_sum_action`
- **consumer example only**: `src/ordinary_qw_operator_assembly.jl`,
  `assembled_one_body_hamiltonian`
- **oracle/reference only**: `src/pqs_multilayer_support_one_body.jl`:
  `pqs_multilayer_support_electron_nuclear_by_center_matrices`
- **retired/do not call**:
  deleted CPBM global electron-nuclear retained-matrix pilot

Do-not-forget rule:
By-center nuclear matrices are uncharged unit attractions, `U_A = -1/r_A`.
Apply `Z_A` and sum centers only when forming the Hamiltonian.
The `HP-R3UN-FN-01` optimization lane covers only terminal final-basis
unit-nuclear `U_GG` Gaussian-sum allocation reduction. It does not approve
neutral raw-block changes, terminal kinetic/moment `G-G` product changes,
route/setup cleanup, residual/MWG/IDA changes, artifacts, public API, or Cr2
workflow.

## R3 Terminal G-G Product Matrices

Why to check:
After neutral nuclear and non-nuclear `G-A`/`A-A` raw-block reuse, Cr2 q4
exact-operator allocation is dominated by terminal final-basis `G-G` product
matrices and unrelated route/stage setup. The approved G-G lane is narrow and
only covers kinetic, coordinate-moment, and second-moment product matrices used
by Residual Gaussian exact augmented operators.

Key docs:
- `docs/src/developer/designs/cartesian_hamiltonian_producer/r3_terminal_gg_product_matrices.md`
- `docs/src/developer/designs/cartesian_hamiltonian_producer/registry.md`,
  `HP-R3GG-FN-01`, `HP-R3GG-TEST-01`

Source anchors:
- **active implementation surface**:
  `src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl`,
  `_r3a_product_matrix`, `pqs_terminal_residual_gto_augmented_operators`
- **active reusable kernel**:
  `src/cartesian_final_basis_realization/pqs_terminal_one_body.jl`,
  `assemble_terminal_product_operator!`

Do-not-forget rule:
This lane is final-basis `G-G` only. Do not change Gaussian `G-A`/`A-A` raw
blocks, nuclear Gaussian sums, IDA/MWG interaction, residual selection or
transforms, route setup, public API, artifacts, status/report fields, or Cr2
workflow.

## IDA And Density Pair-Factor Conventions

Why to check:
Avoid mixing raw/source weights, density-normalized pair factors, retained
weights, and final IDA weights.

Key docs:
- `docs/src/developer/designs/cartesian_hamiltonian_producer/`
- `docs/src/developer/numerical_contracts.md`
- `docs/src/developer/white_lindsey_low_order_density_density_builder_contract_2026-06-04.md`
- `docs/src/developer/raw_product_source_retained_transform_policy.md`
- `docs/src/developer/pqs_source_box_operator_framework.md`

Source anchors:
- **active implementation surface**:
  `src/cartesian_final_basis_realization/pqs_terminal_ida.jl`,
  `assemble_terminal_ida_interaction!`
- **oracle/reference only**: `src/cartesian_nested_faces.jl`,
  `_nested_factorized_weight_aware_pair_terms`,
  `_nested_weight_aware_pair_terms`, `_nested_support_reference_pair_sum`
- **active donor pattern**: `src/ordinary_mapped_backends.jl`,
  `_mapped_coulomb_expanded_symmetric_matrix`
- **active donor pattern**: `src/ordinary_qw_raw_blocks.jl`,
  `_qwrg_fixed_block_interaction_matrix`
- **active donor pattern**: `src/cartesian_contracted_parent_metrics/core.jl`:
  `_pqs_source_box_ida_factor_provenance`
- **consumer example only**: `src/ordinary_cartesian_ida.jl`,
  `_ordinary_cartesian_ida_from_pair_factors`,
  `_ordinary_cartesian_ida_from_gausslet_bundle`
- **retired/do not call**: `src/ordinary_cartesian_ida.jl`,
  `_ordinary_cartesian_ida_from_pair_factors` term-by-term `kron` construction
  when a term-first mapped/terminal contraction is available

Do-not-forget rule:
Use PGDG `pair_factor_terms_raw` as the raw numerator source and keep the
Gaussian expansion index as the short inner reduction. Carry raw pair numerators
through projection/realization first. Normalize only at the reviewed final
retained/final-basis density boundary.

## Base Hamiltonian Construction And Materialization

Why to check:
Use the implemented Slice C2/D handoff without recreating wrapper payloads,
report mirrors, or artifact formats.

Key docs:
- `docs/src/developer/designs/cartesian_hamiltonian_producer/`
- `docs/src/developer/numerical_contracts.md`

Source anchors:
- **active reusable kernel**: `src/cartesian_ida_hamiltonian.jl`,
  `CartesianIDAHamiltonian`, `write_cartesian_ida_hamiltonian`
- **active implementation surface**:
  `src/pqs_source_box_low_order_materialization.jl`,
  `_pqs_source_box_route_driver_terminal_ida_hamiltonian`,
  `_pqs_source_box_route_driver_materialization`
- **active implementation surface**:
  `src/pqs_source_box_route_driver_helpers.jl`,
  `cartesian_materialization(report, terminal_basis_realization,
  materialization_inputs)`
- **approved driver workflow surface**: `bin/cartesian_ham_builder.jl`,
  `HP-DRV-FILE-01` / `HP-DRV-FN-01` / `HP-DRV-ATOM-FN-01`

Do-not-forget rule:
Requested base PQS materialization returns `CartesianIDAHamiltonian{Float64}`
directly; no-request materialization returns `nothing`. Do not add a materialized
wrapper, status mirror, report field, metadata-carried matrix, or new artifact
shape for the base handoff, except the approved `HP-R1-ART-01`
`producer_provenance/` keys stored in the final Hamiltonian file. The canonical
driver may write artifacts and print compact timings/summaries, but must call
approved producer surfaces rather than underscored route internals. That
provenance is for consumer tracking, not a staged algorithm input after initial
lattice/parent construction. The atom driver workflow is base-only through the
existing base facade; it does not approve supplemented atom Hamiltonians or
broader atom source support.

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
- **active donor pattern**: `src/cartesian_nested_owned_units.jl`,
  `_nested_endcap_panel_owned_units`,
  `_nested_endcap_panel_shell_layer`
- **oracle/reference only**: `src/cartesian_high_order_doside_experimental.jl`
- **oracle/reference only**:
  `src/cartesian_high_order_doside_ida_experimental.jl`
- **active donor pattern**: `src/cartesian_nested_faces.jl`, `_nested_doside_1d`

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
- **oracle/reference only**: `src/cartesian_nested_faces.jl`
- **oracle/reference only**: `src/ordinary_qw_raw_blocks.jl`
- **consumer example only**: `src/ordinary_qw_operator_assembly.jl`
- **oracle/reference only**: `src/cartesian_contracted_parent_metrics/`

Do-not-forget rule:
Old code may feed adapters or oracle comparisons. It should not silently become
new route authority.
