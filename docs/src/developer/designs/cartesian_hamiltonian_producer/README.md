# Cartesian Hamiltonian Producer

This page is navigation and orientation for the Cartesian/PQS Hamiltonian
producer. It does not independently grant source authority: `current.md` owns
live status, `authority.toml` owns record-level permissions and surfaces,
`registry.md` is its generated human view, `invariants.md` owns
cross-subsystem guardrails, and linked subsystem pages own their numerical
contracts. On disagreement, follow the fail-closed rule in `invariants.md`.

## Orientation

The producer constructs orthonormal Cartesian working bases and their
one-body/IDA Hamiltonians through one staged pipeline:

1. normalize explicit atom or bond-aligned diatomic inputs;
2. build the mapped parent lattice and route-independent terminal supports;
3. realize those supports as direct identity blocks or PQS/White-Lindsey
   contractions;
4. assemble exact one-body operators and the final-basis IDA interaction;
5. optionally add Gaussian supplements through the Residual Gaussian boundary;
6. optionally form protected members, in-memory reference corrections, or
   transfer sidecars without changing the underlying interaction convention.

The stable conceptual boundaries are:

- [terminal basis and base assembly](terminal_basis_and_base_assembly.md) for
  support ownership, realization, exact operators, and base Hamiltonians;
- [composition](nesting_supplement_composition_plan.md) and
  [Residual Gaussians](residual_gaussian_domain_module.md) for optional
  supplements and MWG augmentation;
- [protected-localized bases](protected_localized_basis.md),
  [reference corrections](screened_hartree_correction_assembly.md), and
  [representation transfer](external_gto_orbital_import.md) for opt-in internal
  consumers;
- [artifact provenance](cartesian_hamiltonian_artifact_manifest.md) and the
  [canonical driver](cartesian_driver_usability_workflow.md) for persisted
  producer output.

Read [current status](current.md) for implemented facilities, active lanes,
physics targets, and blockers. Read [invariants](invariants.md) for
cross-subsystem guardrails. For work on a specific ID, use its generated
[registry entry](registry.md) and linked canonical contract.

## Documentation Map

Normal startup is the compact set in `AGENTS.md`. The list below is a
task-specific contract index; read only the documents relevant to the assigned
ID or subsystem.

- [current.md](current.md)
- [registry.md](registry.md)
- [invariants.md](invariants.md)
- [Terminal basis and base assembly](terminal_basis_and_base_assembly.md) for
  the implemented foundational terminal objects, support-local realization,
  blockwise one-body/IDA matrices, and Hamiltonian construction boundary
- [Historical implementation slices](implementation_slices.md) for the compact
  non-authoritative Slice A/B/C/D disposition and later subsystem index
- [Residual Gaussian domain module](residual_gaussian_domain_module.md) for
  current RG algorithm authority
- [Residual Gaussian orthogonality robustness](residual_gaussian_orthogonality_robustness.md)
  for the current merge/identity checks, production cutoff precedence, and
  superseded tolerance history
- [Default-off direct-G residual injection](residual_gaussian_injection_hybrid.md)
  for the implemented preservation-only compatibility path and its completed
  measurement history
- [Protected-localized basis convention](protected_localized_basis.md)
  for compact-first protected-original replacement, exact localized one-body
  operators, inherited-site `Vee`, and the rejection of direct `C' V C`
- [Retained-GTO local-product EGOI](retained_gto_egoi.md) for the completed
  measurement and approved-pending retained-original `s1+s2`, local-product,
  `M2` correction helper
- [Protected-localized ladder bundles](protected_localized_ladder.md) for the
  implemented same-parent member, cross-overlap, native restart, manifest,
  and target-evaluation contract
- [Occupied-first injection geometry](occupied_first_injection.md)
  for mandatory occupied protection, physical capture validation, optional
  supplement selection, and the current unwired-consumer boundary
- [Protected additive atomic reference correction](protected_additive_reference_correction.md)
  for the internal occupied-union, additive packet `P0/J0/E0`, and native
  protected-localized screened-Hartree correction lane
- [Numerical-complete residual Gaussian basis](numerical_complete_residual_basis.md)
  for the implemented internal path that preserves `G`, appends the full
  numerical supplement complement, and validates packet capture afterward
- [Reference Hartree numerics](reference_hartree_numerics.md) for implemented
  neutral exact `GG/GA/AA` kernels and protected fixed/localized transforms
- [Rho0 and reference-density correction history](rho0_reference_density_matrix.md)
  for superseded fixed-`P0` experiments and the deferred XPAIR question
- [Terminal shellification due diligence](terminal_shellification_due_diligence.md)
  for the derived basis review report and advisory warning contract
- [PQS complete-shell aspect source modes](pqs_complete_shell_aspect_source_modes.md)
  for the separate `(q,q,L)` source-policy lane
- [Semantic per-shell PQS source-q overrides](pqs_semantic_shell_q_overrides.md)
  for the approved-pending private owner-balanced refinement diagnostic
- [Cartesian Gaussian raw blocks - nuclear slice](cartesian_gaussian_raw_blocks_nuclear.md)
  for the neutral uncharged nuclear raw-block owner
- [Cartesian Gaussian raw blocks - non-nuclear slice](cartesian_gaussian_raw_blocks_non_nuclear.md)
  for the neutral overlap/kinetic/moment raw-block owner
- [R3 terminal G-G product matrices](r3_terminal_gg_product_matrices.md)
  for implemented exact terminal kinetic and first/second moment products
- [R3 remaining exact-operator allocation audit](r3_remaining_exact_operator_allocation_audit.md)
  for the measurement-only decision after terminal `G-G` workspace reuse
- [R3 unit-nuclear U_GG Gaussian sum](r3_unit_nuclear_ugg_gaussian_sum.md)
  for implemented exact uncharged by-center terminal `U_GG` assembly
- [R3 same-construction base reuse](r3_same_construction_base_reuse.md)
  for trusted base kinetic/unit-nuclear reuse, validation, fallbacks, and
  canonical-driver call-site behavior
- [Screened Hartree residual-density formalism](screened_hartree_residual_density.md)
  for the durable physics contract: keep `Vnuc_G` Galerkin and use IDA/MWG
  only for residual density fluctuations `q - q0`
- [Atomic HF reference packets](atomic_hf_reference_packets.md)
  for the canonical determinant, density-fit, potential-fit, identity, and
  validation contract of reusable one-center packet references
- [Screened Hartree correction assembly](screened_hartree_correction_assembly.md)
  for the implemented internal `Delta_J0`/`C` correction API built from
  represented references and same-basis `V_IDA`
- [External GTO orbital import](external_gto_orbital_import.md)
  for the representation-transfer facility that imports explicit external AO
  orbitals into an orthonormal final basis by `C_F = <F|G> C_G`, including the
  implemented protected-member composition and standalone native `S_LG` sidecar
- [PQS/WL mapping `s_factor`](pqs_mapping_s_factor.md)
  for the expert mapping-strength scalar that preserves default behavior while
  allowing CR2-style scans
- [Producer-wide Coulomb accuracy](coulomb_accuracy_policy.md)
  for the expert compact/standard/high presets and one-expansion
  construction/provenance
  contract
- [Cartesian driver usability workflow](cartesian_driver_usability_workflow.md)
  for the compact artifact-producing canonical driver lane
- [R1 one-center base atoms](r1_one_center_base_atoms.md)
  for explicit origin-centered all-electron one-center base atom inputs
- [Cartesian driver atom workflow](cartesian_driver_atom_workflow.md)
  for explicit origin-centered base atom driver inputs
- [R3 homonuclear z-axis diatomic supplemented workflow](r3_homonuclear_diatomic_supplemented_workflow.md)
  for the explicit homonuclear diatomic molecule-scope relaxation
- [White-Lindsey terminal basis realization](white_lindsey_terminal_basis_realization.md)
  for the narrow terminal-basis seam needed by `nesting = :wl`
- [Nesting/supplement composition](nesting_supplement_composition_plan.md)
  for the implemented 2 x 2 x 2 composition matrix and dependency order
- [Public ns direct-core side parity](public_ns_core_side_parity.md)
  for deriving direct nucleus-centered core side from public `ns` rather than
  route-local `q`
- [Common terminal shell decomposition](common_terminal_shell_decomposition.md)
  for keeping direct core and shell-owned support construction shared before
  PQS/WL retained-realization geometry diverges
- [Mapped-COMX source span](mapped_comx_source_span.md)
  for the mainline protected-`P2` plus mapped Chebyshev source-span option at
  the existing doside/COMX seam and the high-order consumer validation
  contract, plus the terminal-basis wiring that consumes carried materialized
  axis facts and the compact `source_span` driver selector
- [Cartesian Hamiltonian artifact manifest](cartesian_hamiltonian_artifact_manifest.md)
  for compact JLD2 sidecar groups describing matrix-order basis rows and public
  recipe provenance, plus the narrow construction-native source-mode
  provenance seam
- [Route/stage metadata contract](route_stage_metadata_contract.md)
  for current vector-backed inventory, plan, deterministic-order, and compact
  stage-carrier semantics. The five former cleanup pages are historical
  pointers to this contract.
- [Completed complete-core-shell RHF retirement record](complete_core_shell_rhf_retirement.md)
  records deletion of the stale RHF payload stack
- [Completed route-driver materialization retirement record](route_driver_materialization_retirement.md)
  records deletion of the old route-driver materialization/report/save wrapper
  workflow and stale tool/test pressure
- [Algorithm implementation index](../../algorithm_implementation_index.md)

## Contracts, Amendments, And Completed Records

The task-specific map above is the contract index. The schema-v3
[machine authority](authority.toml) owns each ID's lifecycle, grant, exact
surfaces, dependencies, and canonical references. The generated
[registry](registry.md) and marked `AGENTS.md` whitelist are checked views;
linked pages own behavior. Surrounding policy and contracts may restrict
machine authority but cannot broaden it. Do not copy their specifications into
this navigation page.

The [atomic cutover record](authority_atomic_cutover_plan.md) preserves the
reviewed transition digest, activation boundary, validation, and whole-commit
rollback rule. Earlier shadow/candidate/transition files were removed rather
than retained as a second authority system.

Completed implementation plans remain as short historical pointers only where
their names are still useful for old links. Completed retirement records remain
available for restoration guardrails, not as active implementation authority.

Candidate amendments:

- Translated atoms, Cr2-specific workflow, public supplemented
  workflow/export, basis/supplement-realism beyond explicit supplied
  labels/files, and broad driver diagnostics remain candidate-only until
  separately approved. Supplemented atoms and supplemented WL z-axis diatomics
  are already implemented internal composition cells.

Historical material:

- [history/cartesian_hamiltonian_producer_design_2026-06_full.md](history/cartesian_hamiltonian_producer_design_2026-06_full.md)
  preserves the full June 2026 design document.
- [reviews/README.md](reviews/README.md) preserves design review rounds and reconciliations.

Historical files are useful for context and audit trails, but they are not
normal startup reading and should not override the compact current authority.
Recursive/previous-block projection review, spike, and performance material is
stale: terminal supports are owned local rows, cross-block overlap is
structurally zero, and post-`d2bf139c` Be2 measurements supersede earlier
recursive-projection Be2 measurements.

Strategic context:

- [Cartesian long-range roadmap](../../roadmaps/cartesian_long_range_roadmap.md)
  sequences future public-producer, unification, supplement, high-order, and
  Cr2 validation work. It is not implementation authority.
