# Cartesian Hamiltonian Producer

This page is navigation and orientation for the Cartesian/PQS Hamiltonian
producer. It does not independently grant source authority: `current.md` owns
live status, `authority.toml` owns record-level permissions and surfaces,
`registry.md` is its generated human view, `invariants.md` owns
cross-subsystem guardrails, and linked subsystem pages own their numerical
contracts. On disagreement, follow the fail-closed rule in `invariants.md`.

## Orientation

Status: the foundational
[terminal-basis and base-assembly contract](terminal_basis_and_base_assembly.md)
is implemented for the internal base path. It owns support-local terminal
realization, blockwise exact one-body assembly, localized IDA, and direct
`CartesianIDAHamiltonian` construction. The historical Slice D route-driver
wrapper is retired. The exported R1 public base
producer is implemented for the atom/diatomic scope recorded in
`r1_public_base_producer.md`, with a narrow explicit origin-centered
one-center all-electron atom relaxation recorded in
`r1_one_center_base_atoms.md`. Residual Gaussian basis,
exact augmented operators, and residual MWG/IDA interaction now belong to the
internal `CartesianResidualGaussians` module. The implemented R3 usability
lane provides a non-exported supported facade for H2 and internal/performance-supported
Be2 artifacts. The neutral Cartesian Gaussian raw-block owner is implemented
for the [uncharged by-center nuclear slice](cartesian_gaussian_raw_blocks_nuclear.md)
and the exact [non-nuclear overlap/kinetic/moment slice](cartesian_gaussian_raw_blocks_non_nuclear.md).
The terminal `G-G` product and by-center unit-nuclear `U_GG` optimization
layers are implemented, including trusted same-construction base-block reuse
in supported supplemented callers.
The canonical Cartesian driver is implemented so the standard
driver can directly produce base and supported supplemented Hamiltonian
artifacts from visible public `system`, `basis`, and optional `supplement`
contracts. `HP-DRV-INV-FN-01` implements a bounded terminal-region inventory
summary
in the canonical driver output so users can see region kind, support rows,
final columns, compression, and identity-vs-compact realization without
turning the driver into a route diagnostic dump. `HP-DRV-SHELLDD-FN-01`
implements
a standard terminal due-diligence report so consumers can review derived
system/geometry facts, parent axes and 1D centers, gausslet/IDA weight stats,
dimension accounting, shell-by-shell physical aspect ratios, source-mode
shapes, retained counts, and warning flags before interpreting energies or
residual/injection behavior. `HP-PQS-ASPECTSHELL-FN-01` implements the
separate
construction policy that selects explicit aspect-aware `(q,q,L)` source
modes for PQS shared complete shells.
`HP-PQS-SHELLQ-OVERRIDE-*` separately approves a pending private diagnostic
that raises transverse source `q` for selected semantic atom-local or shared
complete shells while preserving route `q`, parent geometry, support, and
ordinary defaults. It is limited to ordinary-span numerical-complete
additive-reference composition and adds no mapped-COMX, public, or artifact
control.
The implemented [protected-localized basis convention](protected_localized_basis.md)
owns compact-main replacement, localized `L`, exact `H1_L`, and inherited-site
`Vee_L`. The implemented
[protected-localized artifact contract](protected_localized_artifact.md) owns
its opt-in persistence, native-order locality metadata, compatibility, and
readback failure behavior.
The retained-GTO EGOI audit is complete;
[`HP-RG-PROTECT-EGOI-*`](retained_gto_egoi.md) records the approved but
unimplemented retained-original-GTO local-product helper and validation
contract. It is not an artifact or solver workflow.
The same-parent ladder audit is complete. The implemented
[protected-localized ladder bundle](protected_localized_ladder.md) owns the
opt-in directory manifest, protected member references, exact final-basis
cross-overlap sidecars, optional native-order restarts, and bounded summaries.
`HP-RG-PROTECT-ADDREF-*` implements a separate internal in-memory lane that
makes the full union of placed atomic packet occupied spaces mandatory in the
protected basis, preserves the original per-packet blocks for additive `P0`,
and assembles native-order screened-Hartree `Delta_J0/C` without changing the
protected artifact or inherited-site `Vee_L` convention.
`HP-RG-NUMCOMP-*` separately implements an internal opt-in lane that preserves
`G`, retains the owner-local numerical supplement
complement at `eta_num = 1e-10`, validates packet occupied capture after
construction, and assembles the same in-memory additive correction in native
unlocalized `[G,R_num]` order. It does not change the production `1e-6`
residual cutoff.
The supported supplemented
workflow now accepts explicit homonuclear
z-axis diatomics without element-specific branches. The implemented compact
artifact manifest adds matrix-order labels and recipe provenance without
changing matrix keys or reader behavior. Its construction-native source-mode
seam remains partial: terminal shell/mode and retained boundary-seed facts are
written, while direct/residual relations and ray/radial labels remain
unavailable. Broad public
API/export, translated atoms, Cr2-specific production workflow, ECP, broad
EGOI workflow, RHF, and solver work remain deferred.
The bounded 2 x 2 x 2 matrix over geometry, nesting, and supplement state is
implemented for origin-centered atoms and homonuclear z-axis diatomics.
[Nesting/supplement composition](nesting_supplement_composition_plan.md) owns
the shared producer path, public `ns` and derived route-local `q`, the WL
diatomic `ns < 4` guard, and retained-support saturation. One-center
`radius` is physical box-extent authority, while
[public-`ns` direct-core parity](public_ns_core_side_parity.md) limits
oddization to direct nucleus-centered core blocks.
`HP-PQS-MAP-SFACTOR-*` implements a narrow expert mapping-strength scalar:
`s_factor`, default `1.0`, with one-center
`effective_s = s_factor * sqrt(Z * core_spacing)` and explicit provenance.
This is the only approved public mapping-strength knob and does not revive
public `d` or `parent_mapping_d`.
`HP-PQS-COULOMB-ACCURACY-*` implements the compact45 and high135 producer
tiers. The fixed analytic K60 `:standard` tier and canonical-driver exposure
are approved but not yet implemented; the default remains `:compact` and the
driver currently does not expose the policy. Once implemented, `:standard`
is the recommended accuracy/cost opt-in while `:high` remains the reference
tier. One
resolved `CoulombGaussianExpansion` must govern parent/PGDG, base IDA,
residual-GTO, and MWG construction, with one Hamiltonian-wide provenance
summary. Its narrow stability amendment replaces cancellation-prone analytic
Gaussian determinant/weighted-variance forms; it does not add a scaled PGDG
carrier or terminal-contraction redesign.
The WL z-axis diatomic compact retained-basis correction is implemented under
`HP-WLDIAT-COMPACT-FN-01`; boundary strata are compact products rather than
full-support identity rows.
The follow-up WL boundary-stratum parity cleanup is implemented under
`HP-WLDIAT-PARITY-FN-01`: direct nucleus-centered core blocks keep odd-side
centering, but boundary shells retain the requested shell count such as
`ns = 4 -> 4^3 - 2^3 = 56`.
Common terminal shell decomposition is implemented as a shared route-family-free
first step under `HP-COMP-SHELLGEOM-FN-01`: PQS and White-Lindsey differ only
after
common shell/core support regions exist.
For z-axis diatomics, `HP-COMP-SHELLGEOM-DIAT-FN-01` requires PQS and
White-Lindsey to call the common shellifier with the same first-step arguments;
central-gap/contact and shared-shell ownership are common shell geometry.
`HP-COMP-THINSLAB-FN-01` implements the matching lowering guardrail for
midpoint slabs
and outer-mismatch slab stacks: PQS and White-Lindsey use the same compact
slab lowering, with unit-slice scale `ns x ns x 1`, and neither may retain
those slabs as full identity rows.
`HP-COMP-THINSLAB-META-FN-01` implements the matching metadata/scaffold
inventory
cleanup so terminal-shellification summaries no longer claim midpoint,
outer-mismatch, or angular z-extension slabs are direct identity sectors.
`HP-COMP-FACEPROD-FN-01` implements the neutral terminal face-product helper seam
needed by that lowering: WL facets and thin slabs reuse the same
module-internal face-like product block builder instead of putting shared
coefficient assembly in a PQS- or WL-owned terminal file.
`HP-COMP-ANGBOX-FN-01` implements the shellification side of that correction:
z-axis diatomic shared-shell growth may emit planned angular z-extension slab
stacks when the ordinary index-layer shell body underreaches the physical
outer-nucleus angular target. Lowering those slabs remains governed by
`HP-COMP-THINSLAB-FN-01`.
Mapped-COMX source spans, terminal consumption of carried axis facts, and the
compact `source_span = :ordinary | :mapped_comx` driver choice are
implemented under the six `HP-MCOMX` source IDs. Ordinary remains default;
mapped-COMX is
PQS-only and high-order remains a consumer/benchmark lane.
The old complete-core-shell RHF payload stack was removed by `28e9b2c84` under
`HP-RETIRE-CCS-RHF-*`. The route-driver materialization/report/save wrappers
were removed by `e2e164e9b`, and the dangling ladder runners by `77fa2700b`.
Those IDs are historical deletion records; the canonical producer path is the
staged human-facing driver plus `CartesianIDAHamiltonian` artifacts.
The implemented [reference Hartree numerics](reference_hartree_numerics.md)
remain neutral exact `GG/GA/AA` and protected-transform infrastructure. Old
row-gauge, fixed-`P0`, FAPP, anchor, and corrected-Hamiltonian experiments are
compressed in
[rho0 and reference-density correction history](rho0_reference_density_matrix.md).
Live screened-Hartree and additive-reference policy belongs to their separate
canonical pages; deferred `HP-RHO0-XPAIR-AUDIT-01` is not a current source
lane or blocker.

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
