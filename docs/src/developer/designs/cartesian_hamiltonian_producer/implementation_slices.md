# Historical Implementation Slices

Document role: compact historical implementation and decision index. This file
is not live status and does not independently authorize work. Current
permission and lifecycle live in `registry.md`; current numerical law lives in
the linked canonical subsystem documents. Git history preserves the former
2,050-line implementation ledger.

## Foundational Slice Disposition

| Slice | Final disposition | Canonical contract | Evidence |
| --- | --- | --- | --- |
| A - terminal basis realization | Implemented, then corrected to support-local shell realization with structurally disjoint owned rows | [Terminal basis and base assembly](terminal_basis_and_base_assembly.md) | `082c9cb8b`, `d2bf139c6`, `d6968d15b`; manager Passes 9-10 and later structural correction |
| B - final-basis one-body operators | Implemented blockwise product and term-first Gaussian-sum assembly | [Terminal basis and base assembly](terminal_basis_and_base_assembly.md) | `3d39a158c`, `ba322abfa`; manager Passes 13-24 |
| C1 - localized IDA matrix | Implemented blockwise final-weight-normalized IDA assembly | [Terminal basis and base assembly](terminal_basis_and_base_assembly.md) | `a33842fb8`; manager Pass 26 |
| C2 - `CartesianIDAHamiltonian` construction | Implemented through the existing Hamiltonian object and current staged base assembly | [Terminal basis and base assembly](terminal_basis_and_base_assembly.md) | manager Passes 27-28; current base/IDA tests |
| D - route-driver materialization handoff | Historical wrapper implemented, later retired and deleted | [Route-driver materialization retirement](route_driver_materialization_retirement.md) | `e2e164e9b`; completed `HP-RETIRE-DRV-MAT-*` |

## Durable Lessons

- Terminal blocks own disjoint support-local parent rows. Previous-block and
  recursive projection were rejected after they enlarged shell support and
  obscured ownership.
- Structural cross-block overlap is zero. Cross-block kinetic, nuclear, and
  IDA matrix elements are generally nonzero and remain assembled blockwise.
- Direct blocks stay implicit identity. Compact blocks carry only support-local
  coefficients and native final column ranges.
- Shell orthonormalization is local and uses the established symmetric Lowdin;
  a global parent/final overlap and global Lowdin are not producer paths.
- Slice B/C produce real final-basis matrices. They do not create payload,
  cache, report, or wrapper frameworks.
- The old Slice D `cartesian_materialization(...)` wrapper is not a compatibility
  surface and must not be restored.

## Later Implementation Index

The old ledger accumulated many independent subsystems after Slices A-D. Their
current contracts are listed here instead of repeated.

| Historical ledger topic | Current canonical owner |
| --- | --- |
| White-Lindsey terminal basis seam and compact retained correction | [White-Lindsey terminal basis realization](white_lindsey_terminal_basis_realization.md) |
| Residual Gaussian module migration, exact augmented operators, and MWG | [Residual Gaussian domain module](residual_gaussian_domain_module.md), [R3 compatibility history](r3_residual_gto_mwg_augmentation.md) |
| Numerical-complete residual supplement | [Numerical-complete residual basis](numerical_complete_residual_basis.md) |
| Compact Hamiltonian artifact manifest | [Cartesian Hamiltonian artifact manifest](cartesian_hamiltonian_artifact_manifest.md) |
| Cartesian Gaussian nuclear raw blocks | [Nuclear raw blocks](cartesian_gaussian_raw_blocks_nuclear.md) |
| Cartesian Gaussian overlap/kinetic/moment raw blocks | [Non-nuclear raw blocks](cartesian_gaussian_raw_blocks_non_nuclear.md) |
| R3 terminal `G-G` product-matrix optimization | [R3 terminal G-G product matrices](r3_terminal_gg_product_matrices.md) |
| R3 exact-operator allocation attribution | [R3 remaining exact-operator allocation audit](r3_remaining_exact_operator_allocation_audit.md) |
| R3 terminal unit-nuclear Gaussian-sum optimization | [R3 unit-nuclear U_GG Gaussian sum](r3_unit_nuclear_ugg_gaussian_sum.md) |
| Same-construction base K/U reuse and supplemented diatomics | [R3 homonuclear supplemented workflow](r3_homonuclear_diatomic_supplemented_workflow.md), [R3 usability workflow](r3_usability_supplemented_workflow.md) |
| Canonical driver and artifact-producing workflow | [Cartesian driver usability workflow](cartesian_driver_usability_workflow.md) |
| Driver atom input and hidden-`d` cleanup | [Cartesian driver atom workflow](cartesian_driver_atom_workflow.md) |
| R1 one-center base atoms | [R1 one-center base atoms](r1_one_center_base_atoms.md) |
| PQS complete-shell aspect source modes | [PQS aspect source modes](pqs_complete_shell_aspect_source_modes.md) |
| Atom/diatomic, PQS/WL, base/supplemented composition | [Nesting/supplement composition](nesting_supplement_composition_plan.md) |
| Public `ns` and direct-core odd-side parity | [Public ns direct-core parity](public_ns_core_side_parity.md) |
| Common shell geometry, thin slabs, and face products | [Common terminal shell decomposition](common_terminal_shell_decomposition.md) |
| Complete-core-shell RHF removal | [Complete-core-shell RHF retirement](complete_core_shell_rhf_retirement.md) |
| Route recipe family cleanup | [Nesting/supplement composition](nesting_supplement_composition_plan.md), registry `HP-ROUTE-RECIPE-*` |
| Route inventory cleanup | [Route inventory type-surface cleanup](route_inventory_type_surface_cleanup.md) |
| Raw product source-mode inventory cleanup | [Raw product source-mode cleanup](raw_product_source_mode_inventory_cleanup.md) |
| Contract-plan vector cleanup | [Contract-plan vector cleanup](contract_plan_vector_cleanup.md) |
| Route/stage type cleanup | [Route/stage type cleanup](route_stage_type_surface_cleanup.md) |
| Route/stage carrier cleanup | [Route/stage carrier cleanup](route_stage_carrier_cleanup.md) |
| Homonuclear z-axis supplemented workflow | [R3 homonuclear supplemented workflow](r3_homonuclear_diatomic_supplemented_workflow.md) |
| Protected basis, artifacts, EGOI, ladders, and additive references | [Protected basis](protected_localized_basis.md), [artifact](protected_localized_artifact.md), [EGOI](retained_gto_egoi.md), [ladder](protected_localized_ladder.md), [additive reference](protected_additive_reference_correction.md) |
| Atomic packets and screened Hartree | [Atomic packets](atomic_hf_reference_packets.md), [screened-Hartree formalism](screened_hartree_residual_density.md), [correction assembly](screened_hartree_correction_assembly.md) |

## Rejected And Retired Paths

### Recursive Terminal Projection

The early Slice A implementation projected each new shell against prior blocks
and allowed effective supports to grow. That approach is superseded. Current
terminal blocks remain on their authoritative owned rows; structural overlap
failures identify construction bugs rather than projection work.

### Global Terminal Repair

No global terminal overlap, global coefficient matrix, global Lowdin, or
cross-block overlap minimization is an accepted producer repair. Shell-local
identity and disjoint support are the relevant contracts.

### Slice D Wrappers

The route-report/materialization/save wrapper family and its dedicated tools
were deleted by `e2e164e9b`. Current production uses staged base construction,
direct `CartesianIDAHamiltonian` assembly, and separately owned artifact
writers. Historical wrapper names carry no authority.

### Direct Interaction Rotation

Later protected-basis measurements rejected treating a two-index IDA/MWG
interaction like a one-body matrix. `C' V C` is not an alternative interaction
construction. The current protected and numerical-complete contracts own their
respective interaction conventions.

## Evidence Lookup

Use the append-only `docs/src/developer/pqs_manager_running_log.md` for accepted
pass interpretation and `git log -- <path>` for source chronology. The
foundational implementation sequence is concentrated around manager Passes
8-33; later entries record structural correction, optimization, retirement,
and subsystem extraction. Do not copy those narratives back into active
contracts.
