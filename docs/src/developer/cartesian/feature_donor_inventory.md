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
- The base H2 PQS Hamiltonian producer now exists and materializes
  `CartesianIDAHamiltonian{Float64}` through the current driver harness. This
  inventory tracks donor features outside that recovered base path.
- Residual-GTO/MWG supplements, corrections, branch/fragment workflows, and
  high-order geometry are future roadmap lanes, not proof that the base PQS
  Hamiltonian endpoint is still blocked.

## Migration Table

| feature | donor files/functions | current surface | what thin route lacks | temporary ladder case | proposed thin-route destination | deletion condition | priority |
|---|---|---|---|---|---|---|---|
| Residual-GTO / MWG supplement materialization | `ordinary_qw_residuals.jl`; `ordinary_qw_raw_blocks.jl`; `ordinary_qw_operator_assembly.jl`; `cartesian_gto_probes.jl`; residual-GTO and mixed gausslet/GTO block helpers | mixed public/internal | The H2-specific residual-GTO materializer has been retired. The base PQS Hamiltonian endpoint now exists; what remains missing is generic final-basis supplement augmentation and consumer validation for mixed gausslet/GTO and MWG terms. | None yet. Future cases should cover a reviewed generic supplement endpoint and external consumer use. | Generic final-basis supplement augmentation paired with the public Cartesian IDA Hamiltonian boundary and any consumer-driven supplement/basis artifact contract. | Delete or shrink donor pieces after a public producer absorbs mixed gausslet/GTO blocks, GTO/GTO blocks, residual representation, MWG interaction kernels, and consumer-validated IDA output, or after the feature is intentionally abandoned. | P1 |
| Ham/JLD2 artifact contract and basis transfer/roundtrip | current owners: `cartesian_basis_representation.jl`, `cartesian_cross_overlap.jl`, `cartesian_representation_transfer.jl`, `CartesianParentGaussletBases.jl`, `CartesianContractedParents.jl`, and current artifact writers | public/mixed | Public Cartesian IDA Hamiltonian read/write and cross-overlap transfer exist. This still does not provide a general public basis/supplement artifact contract. | `tools/h2_pqs_base_hamiltonian_smoke.jl` covers the current base Hamiltonian artifact; future basis/provenance cases should wait for a downstream consumer need. | Extend current representation, overlap, transfer, parent, and artifact owners only when a real consumer requires it. | The orphaned `cartesian_carried_spaces.jl` adapter is not a donor; retire it under `HP-RETIRE-CARRIED-SPACE-*`. | P1 |
| Hydrogenic-core / ESOI corrections | `ordinary_qw_corrections.jl`; `HydrogenicCoreProjectorCorrectionSpec`; `HydrogenicCoreBranchCorrectionSpec`; `HydrogenicCoreCorrectionSpec`; `include_esoi`; `apply_ordinary_cartesian_corrections`; `ordinary_cartesian_corrected_branch` | public | No thin-route correction stage for hydrogenic-core projector, branch-local correction, or ESOI-style correction after one-body assembly. | `h2_wl_branch_correction` or a smaller `he_wl_core_correction` case. | Thin route correction stage after H1 assembly and before optional branch/fragment materialization. | Delete donor implementation after thin route exposes equivalent correction semantics, or explicitly retire the public correction feature. | P2 |
| EGOI / density-density correction | `hamiltonian_corrections.jl`; `EGOIDensityDensityCorrectionResult`; `egoi_target_product_matrix`; `egoi_target_coulomb_matrix`; `egoi_density_density_correction`; `egoi_stationary_hamiltonian_correction`; `ordinary_cartesian_egoi_stationary_correction` | public | Thin H1-J/RHF path has no EGOI correction hook or driver-owned correction artifact. | `he_wl_egoi_correction` or `h2_wl_egoi_correction`. | Post-H1/H1-J correction stage in thin materialization, with compact correction summary and optional saved artifact. | Delete donor path only after thin route reproduces or intentionally drops the EGOI correction capability. | P2 |
| Branch / fragment Hamiltonian workflow | `ordinary_qw_operator_assembly.jl`; `ordinary_qw_corrections.jl`; branch nuclear-charge and fragment/counterpoise helpers | public/mixed | Thin route does not yet express branch nuclear-charge variants, fragment Hamiltonians, or counterpoise-style corrected branches. | `h2_wl_fragment_branch` with pure gausslet first; supplement branch later. | Route-level branch materialization mode using the same basis/operator artifact contract. | Delete donor path after branch/fragment workflow is represented by thin driver cases or abandoned. | P3 |
| High-order slab/endcap/panel geometry variants | retiring `cartesian_high_order_doside_experimental.jl`, `cartesian_high_order_doside_ida_experimental.jl`, and `ordinary_qw_experimental_paths.jl`; preserved adjacent helpers in `cartesian_nested_owned_units.jl`, `cartesian_nested_faces.jl`, and `cartesian_nested_diatomic.jl` | retirement approved | The old experimental implementation did not establish a useful distorted-parent route. Endcap/panel concepts that survived are already owned by adjacent active code. | None. Future chain/high-order work requires a new design. | Preserved route-owned geometry helpers plus historical evidence. | Delete the approved cluster under `HP-RETIRE-QW-DONOR-*`; do not add adapters. | closed |
| Legacy nested fixed-source oracles | `cartesian_nested_faces.jl`; `cartesian_nested_atomic.jl`; `cartesian_nested_diatomic.jl`; `ordinary_qw_nested_frontends.jl`; `bond_aligned_diatomic_nested_fixed_source`; `one_center_atomic_full_parent_fixed_block` | exported/internal | Thin route has independent PQS and WL ladders, but some old fixed-source paths still serve as references, migration scaffolds, or public historical APIs. | Keep current WL/PQS line ladders; add a narrow oracle-comparison ladder only if a migration pass needs one. | Typed route-owned support/source construction, with old oracle use removed from active authority paths. | Delete or quarantine donor paths once no current migration or public workflow needs them. | P4 |
| QW carried-space / receipt wrappers | retired `cartesian_qw_operator_carried_spaces.jl`; preserved `cartesian_qw_hybrid_representation.jl`; orphaned `cartesian_carried_spaces.jl` approved for retirement | receipt cluster retired; final adapter pending | The receipt cluster is gone, and the remaining adapter has no committed caller. | None. | Current producer artifacts and preserved representation/overlap/transfer helpers. | Delete the adapter without stubs under `HP-RETIRE-CARRIED-SPACE-*`. | closed |

## Suggested Priority Order

1. P1 residual-GTO / MWG supplement materialization: the H2-specific
   residual-GTO materializer has been retired. Public `CartesianIDAHamiltonian`
   read/write remains live, and the H2 base PQS endpoint now materializes that
   artifact. Residual-GTO/MWG supplement work is future generic final-basis
   augmentation work.
2. P1 Ham/JLD2 artifact contract and basis transfer/roundtrip: the Ham side has
   a public Cartesian IDA writer/reader. The current base H2 PQS smoke covers
   the Hamiltonian artifact. Remaining work is downstream consumer coverage and
   any explicitly needed basis/provenance contract, not private H2 sidecar
   revival or private H1-J/RHF solver diagnostics.
3. P2 hydrogenic-core / ESOI corrections.
4. P2 EGOI / density-density correction.
5. P3 branch / fragment Hamiltonian workflow.
6. P4 legacy nested fixed-source oracles, only if still needed as references.

The old high-order/receipt cluster is no longer a migration priority. Its
retirement is approved by `HP-RETIRE-QW-DONOR-*`; surviving endcap/panel
helpers remain active adjacent code, and atomic chains remain a future target.

## Shared Kernel Status

- Weighted-Hadamard and 1D Gaussian axis-table kernels now live in neutral
  private helpers under `cartesian_gaussian_axis_integrals.jl`.
- Shared QW wrapper names remain only where surviving donor callers use them;
  the caller-free receipt-wrapper names are approved for retirement.
- Remaining donor dependencies include MWG/residual interaction kernels and
  larger 3D support/product contractions.

## Guardrails

- Do not resurrect old status/readiness/probe/helper-schema tests.
- Do not add new helper-payload schema tests for donor migration.
- Each donor feature gets a tiny driver ladder case before or while migrating.
- A ladder case should check coarse semantic or numerical facts, not exact
  print lines or private payload field lists.
- Once a feature is migrated, delete the donor line or mark it intentionally
  retained with a clear reason.
- Do not treat surviving WL/QW, residual-GTO, correction, branch, or high-order
  donor code as dead merely because the current thin route does not call it.
  The explicitly audited `HP-RETIRE-QW-DONOR-*` cluster is the bounded
  exception.
