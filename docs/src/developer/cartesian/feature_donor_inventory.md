# Cartesian Feature Donor Migration

This note tracks old Cartesian feature lines that must be absorbed into the
thin Cartesian driver route before their old implementation lines can be
deleted.

## Current Policy

- The thin staged Cartesian driver is the target architecture.
- Old Cartesian lines are feature donors, not the target architecture.
- Do not delete a donor line until its feature is migrated, explicitly
  abandoned, or proven unused.
- Validation for this private route should use temporary driver ladders, not
  helper/schema `Test.jl` files.
- Donor migration should preserve numerical capability without resurrecting old
  report, readiness, probe, status, blocker, availability, or helper-payload
  schema machinery.

## Migration Table

| feature | donor files/functions | current surface | what thin route lacks | temporary ladder case | proposed thin-route destination | deletion condition | priority |
|---|---|---|---|---|---|---|---|
| Residual-GTO / MWG supplement materialization | `ordinary_qw_residuals.jl`; `ordinary_qw_raw_blocks.jl`; `ordinary_qw_operator_assembly.jl`; `cartesian_gto_probes.jl`; residual-GTO and mixed gausslet/GTO block helpers | mixed public/internal | Thin route can materialize a PQS Ham/Basis plus residual-GTO sidecar artifact for the H2 independent PQS route. It still lacks full mixed/GTO provider blocks, residual MWG supplemented values, and a full provider-block Hamiltonian. | `h2_pqs_q5_independent_source_box_r4_gto_materialized` in the `pqs_diatomic` line ladder; future cases should cover provider blocks and supplemented values. | CPB provider-block supplement stage plus compact supplement artifact contract. | Delete or shrink donor pieces after thin route builds mixed gausslet/GTO blocks, GTO/GTO blocks, residual representation, and supplemented values, or after the feature is intentionally abandoned. | P1 |
| Ham/JLD2 artifact contract and basis transfer/roundtrip | `cartesian_basis_representation.jl`; `cartesian_cross_overlap.jl`; `cartesian_representation_transfer.jl`; `cartesian_carried_spaces.jl`; current JLD2 writer patterns | public/mixed | Thin route now writes/reloads H2 PQS Ham/Basis plus residual-GTO sidecar artifacts. It still lacks a broader downstream consumer roundtrip and artifact contract for provider-block supplemented values. | `h2_pqs_q5_independent_source_box_r4_gto_materialized`; later add `he_wl_artifact_roundtrip` and downstream H2 PQS consumer roundtrip cases. | Compact driver artifact writer for basis, H1/H1-J/RHF/supplement/provenance plus transfer smoke. | Donor line can be deleted or retained narrowly once thin artifacts can be loaded and used by downstream consumers without old representation wrappers. | P1 |
| Hydrogenic-core / ESOI corrections | `ordinary_qw_corrections.jl`; `HydrogenicCoreProjectorCorrectionSpec`; `HydrogenicCoreBranchCorrectionSpec`; `HydrogenicCoreCorrectionSpec`; `include_esoi`; `apply_ordinary_cartesian_corrections`; `ordinary_cartesian_corrected_branch` | public | No thin-route correction stage for hydrogenic-core projector, branch-local correction, or ESOI-style correction after one-body assembly. | `h2_wl_branch_correction` or a smaller `he_wl_core_correction` case. | Thin route correction stage after H1 assembly and before optional branch/fragment materialization. | Delete donor implementation after thin route exposes equivalent correction semantics, or explicitly retire the public correction feature. | P2 |
| EGOI / density-density correction | `hamiltonian_corrections.jl`; `EGOIDensityDensityCorrectionResult`; `egoi_target_product_matrix`; `egoi_target_coulomb_matrix`; `egoi_density_density_correction`; `egoi_stationary_hamiltonian_correction`; `ordinary_cartesian_egoi_stationary_correction` | public | Thin H1-J/RHF path has no EGOI correction hook or driver-owned correction artifact. | `he_wl_egoi_correction` or `h2_wl_egoi_correction`. | Post-H1/H1-J correction stage in thin materialization, with compact correction summary and optional saved artifact. | Delete donor path only after thin route reproduces or intentionally drops the EGOI correction capability. | P2 |
| Branch / fragment Hamiltonian workflow | `ordinary_qw_operator_assembly.jl`; `ordinary_qw_corrections.jl`; branch nuclear-charge and fragment/counterpoise helpers | public/mixed | Thin route does not yet express branch nuclear-charge variants, fragment Hamiltonians, or counterpoise-style corrected branches. | `h2_wl_fragment_branch` with pure gausslet first; supplement branch later. | Route-level branch materialization mode using the same basis/operator artifact contract. | Delete donor path after branch/fragment workflow is represented by thin driver cases or abandoned. | P3 |
| High-order slab/endcap/panel geometry variants | `cartesian_high_order_doside_experimental.jl`; `cartesian_high_order_doside_ida_experimental.jl`; `cartesian_nested_owned_units.jl`; high-order pieces in `cartesian_nested_diatomic.jl`; `ordinary_qw_experimental_paths.jl` | internal/experimental | Thin route does not yet express high-order slab/endcap/panel geometry variants as route-owned shellification/lowering choices. | Future `h2_high_order_pure` or smaller high-order atomic/diatomic geometry ladder. | Route-owned geometry/shellification variant with typed support regions and retained-unit construction. | Delete donor paths once high-order concepts are ported, intentionally abandoned, or replaced by a clearer route-owned implementation. | P3 |
| Legacy nested fixed-source oracles | `cartesian_nested_faces.jl`; `cartesian_nested_atomic.jl`; `cartesian_nested_diatomic.jl`; `ordinary_qw_nested_frontends.jl`; `bond_aligned_diatomic_nested_fixed_source`; `one_center_atomic_full_parent_fixed_block` | exported/internal | Thin route has independent PQS and WL ladders, but some old fixed-source paths still serve as references, migration scaffolds, or public historical APIs. | Keep current WL/PQS line ladders; add a narrow oracle-comparison ladder only if a migration pass needs one. | Typed route-owned support/source construction, with old oracle use removed from active authority paths. | Delete or quarantine donor paths once no current migration or public workflow needs them. | P4 |
| QW carried-space / receipt wrappers | `cartesian_qw_operator_carried_spaces.jl`; `cartesian_qw_hybrid_representation.jl`; QW operator construction receipt helpers | mixed | Thin route uses compact staged objects and lacks the older carried-space receipt/consumer wrapper semantics. | Add only if a downstream consumer still requires a receipt-shaped object. | Compact driver artifact and consumer contract, not a large receipt wrapper. | Delete wrappers after consumers use thin artifacts directly, or mark a small retained adapter as intentional. | P4 |

## Suggested Priority Order

1. P1 residual-GTO / MWG supplement materialization: completed first vertical
   sidecar slice; remaining work is provider blocks, supplemented values, and
   a full provider-block Hamiltonian.
2. P1 Ham/JLD2 artifact contract and basis transfer/roundtrip: H2 PQS
   Ham/Basis plus residual-GTO sidecar artifacts write/reload; remaining work
   is artifact/consumer roundtrip and provider-block artifact coverage.
3. P2 hydrogenic-core / ESOI corrections.
4. P2 EGOI / density-density correction.
5. P3 branch / fragment Hamiltonian workflow.
6. P3 high-order slab/endcap/panel geometry variants.
7. P4 legacy nested fixed-source oracles, only if still needed as references.

## Guardrails

- Do not resurrect old status/readiness/probe/helper-schema tests.
- Do not add new helper-payload schema tests for donor migration.
- Each donor feature gets a tiny driver ladder case before or while migrating.
- A ladder case should check coarse semantic or numerical facts, not exact
  print lines or private payload field lists.
- Once a feature is migrated, delete the donor line or mark it intentionally
  retained with a clear reason.
- Do not treat WL/QW, residual-GTO, correction, branch, or high-order donor code
  as dead merely because the current thin route does not call it yet.
