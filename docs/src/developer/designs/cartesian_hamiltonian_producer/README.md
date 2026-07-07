# Cartesian Hamiltonian Producer

This directory is the compact current authority for the Cartesian/PQS
Hamiltonian producer.

Status: Slice A, Slice B, Slice C1, Slice C2, and Slice D base handoff are
implemented for the internal base PQS Hamiltonian path. R1 public base
producer implementation is approved for the H/H2 scope recorded in
`r1_public_base_producer.md`, with a narrow explicit origin-centered
one-center all-electron atom relaxation recorded in
`r1_one_center_base_atoms.md`. Residual Gaussian basis,
exact augmented operators, and residual MWG/IDA interaction now belong to the
internal `CartesianResidualGaussians` module. The R3 usability lane approves
only a non-exported supported facade for H2 and internal/performance-supported
Be2 artifacts. A neutral Cartesian Gaussian raw-block owner is approved for the
uncharged by-center nuclear slice and the narrow non-nuclear
overlap/kinetic/moment slice. A narrow R3/RG terminal `G-G` product-matrix
optimization lane is approved separately. A narrow unit-nuclear `U_GG`
Gaussian-sum allocation lane is approved after the remaining-allocation audit.
The canonical Cartesian driver usability lane is approved so the standard
driver can directly produce base and supported supplemented Hamiltonian
artifacts from visible public `system`, `basis`, and optional `supplement`
contracts. `HP-DRV-INV-*` approves a bounded terminal-region inventory summary
in the canonical driver output so users can see region kind, support rows,
final columns, compression, and identity-vs-compact realization without
turning the driver into a route diagnostic dump. `HP-DRV-SHELLDD-*` approves
a standard terminal due-diligence report so consumers can review derived
system/geometry facts, parent axes and 1D centers, gausslet/IDA weight stats,
dimension accounting, shell-by-shell physical aspect ratios, source-mode
shapes, retained counts, and warning flags before interpreting energies or
residual/injection behavior. `HP-PQS-ASPECTSHELL-*` separately approves the
future source-policy lane that may change PQS complete shells from cubic
`(q,q,q)` source modes to explicit aspect-aware `(q,q,L)` source modes.
`HP-RG-PROTECT-ART-*` approves a narrow opt-in protected-localized injection
Hamiltonian artifact variant so the accepted `L`, `H1_L`, inherited-site
`Vee_L` convention can be read back without in-memory reconstruction.
`HP-RG-PROTECT-ARTLOC-*` approves native-order row-locality metadata for that
artifact: diagonal position-expectation centers, deterministic z-order
permutations, and per-row sector labels without changing matrix order.
`HP-RG-PROTECT-EGOI-*` approves the first retained-original-GTO local-product
EGOI helper as an internal in-memory correction primitive, not an artifact or
solver workflow. `HP-RG-PROTECT-LADDER-XFER-AUDIT-01` approves only a
same-parent protected-localized ladder transfer measurement audit, using
final-basis cross overlaps to move occupied orbitals between comparable
Cr2 ladder bases. `HP-RG-PROTECT-LADDER-BUNDLE-*` approves an opt-in
directory-bundle facility for this recurring ladder workflow: protected
Hamiltonian member artifacts, exact cross-overlap sidecars, optional restart
sidecars, manifest/provenance, and bounded summaries.
The supported supplemented
workflow now accepts explicit homonuclear
z-axis diatomics without element-specific branches. A compact artifact manifest
lane is approved for JLD2 sidecar groups that make written Hamiltonians
self-describing without changing matrix keys or reader behavior. Broad public
API/export, translated atoms, Cr2-specific production workflow, ECP, broad
EGOI workflow, RHF, and solver work remain deferred.
The documented 2 x 2 x 2 matrix over `geometry`, `nesting`, and supplement
state now has approved implementation authority under the current
origin-centered atom and homonuclear z-axis diatomic constraints.
The approved composition IDs are `HP-COMP-WLDIAT-*`, `HP-COMP-BASEDIAT-*`,
`HP-COMP-SUPPWL-*`, and `HP-COMP-SUPPATOM-*`.
The one-center atom parent-sizing correction is approved under
`HP-COMP-ATOMBOX-*`: atom `radius`/driver `padding` is physical box extent
authority. Public size-parameter normalization is approved under
`HP-COMP-NS-*`: user-facing input should name `ns` as the requested
cube/source/nesting size, while route-local `q` is derived from `ns` and
`nesting`. WL z-axis diatomic `ns` cleanup is approved under
`HP-COMP-WLNS-*`: `nesting = :wl`, `Natom = 2`, and normalized `ns < 4`
should reject early, while working WL diatomic `ns` ranges may saturate the
final retained support when the physical parent extent dominates.
`HP-PQS-MAP-SFACTOR-*` approves a narrow expert mapping-strength scalar:
`s_factor`, default `1.0`, with one-center
`effective_s = s_factor * sqrt(Z * core_spacing)` and explicit provenance.
This is the only approved public mapping-strength knob and does not revive
public `d` or `parent_mapping_d`.
The WL z-axis diatomic compact retained-basis correction is approved under
`HP-WLDIAT-COMPACT-*`: the current mechanical boundary-stratum identity path is
not the intended compact WL retained basis and must not be used as the final
PQS/WL comparison story.
The follow-up WL boundary-stratum parity cleanup is approved under
`HP-WLDIAT-PARITY-*`: direct nucleus-centered core blocks keep odd-side
centering, but boundary shells retain the requested shell count such as
`ns = 4 -> 4^3 - 2^3 = 56`.
Public-`ns` direct core side parity is approved under `HP-COMP-NSCORE-*`:
direct nucleus-centered core blocks use the oddized public `ns` side, while
boundary retained sizes keep route-local construction.
Common terminal shell decomposition is approved as a shared route-family-free
first step under `HP-COMP-SHELLGEOM-*`: PQS and White-Lindsey differ only after
common shell/core support regions exist.
For z-axis diatomics, `HP-COMP-SHELLGEOM-DIAT-*` requires PQS and
White-Lindsey to call the common shellifier with the same first-step arguments;
central-gap/contact and shared-shell ownership are common shell geometry.
`HP-COMP-THINSLAB-*` adds the matching lowering guardrail for midpoint slabs
and outer-mismatch slab stacks: PQS and White-Lindsey use the same compact
slab lowering, with unit-slice scale `ns x ns x 1`, and neither may retain
those slabs as full identity rows.
`HP-COMP-THINSLAB-META-*` approves the matching metadata/scaffold inventory
cleanup so terminal-shellification summaries no longer claim midpoint,
outer-mismatch, or angular z-extension slabs are direct identity sectors.
`HP-COMP-FACEPROD-*` approves the neutral terminal face-product helper seam
needed by that lowering: WL facets and future thin slabs should reuse the same
module-internal face-like product block builder instead of putting shared
coefficient assembly in a PQS- or WL-owned terminal file.
`HP-COMP-ANGBOX-*` approves the shellification side of that correction:
z-axis diatomic shared-shell growth may emit planned angular z-extension slab
stacks when the ordinary index-layer shell body underreaches the physical
outer-nucleus angular target. Lowering those slabs remains governed by
`HP-COMP-THINSLAB-*`.
Mapped-COMX source spans are approved as a narrow option at the existing
doside/COMX source-span seam under `HP-MCOMX-*`; high-order remains a
consumer/benchmark lane, not the owner of the installed facility.
`HP-MCOMX-TERM-*` approves only terminal-basis wiring so carried materialized
mapped-COMX axis facts can define PQS shell seed coefficients. `HP-MCOMX-DRV-*`
approves only a compact `source_span` driver/facade construction choice for
`:ordinary` versus `:mapped_comx`.
The old complete-core-shell RHF payload stack is approved for retirement under
`HP-RETIRE-CCS-RHF-*`; the current CR2-facing producer path is the canonical
driver plus `CartesianIDAHamiltonian` artifacts, not the stale RHF payload
workflow.
The old route-driver materialization/report/save wrapper workflow is approved
for retirement/quarantine under `HP-RETIRE-DRV-MAT-*`; the canonical producer
path is the staged human-facing driver plus `CartesianIDAHamiltonian` artifacts.
The row-gauge `rho0/Galerkin` measurement lane remains historical evidence
only. The current screened-correction target is the reference-density-matrix
path in `rho0_reference_density_matrix.md`: fixed `P0`, exact Hartree side,
direct-only approximate Hartree anchor `Delta_J0`/`C0_J`, and later
small-system corrected behavior. The older full-interaction anchor is
superseded for Hartree-correction interpretation. H/Be/Be2 evidence now makes
the remaining blocker exchange/direct pairing, not unrepresented `P0`;
`HP-RHO0-XPAIR-AUDIT-01` is the measurement-only next lane. Artifact, solver,
public workflow, Cr2 production, and exact exchange correction remain
deferred.

Agents should read first:

- [current.md](current.md)
- [registry.md](registry.md)
- [invariants.md](invariants.md)
- [implementation_slices.md](implementation_slices.md)
- [Residual Gaussian domain module](residual_gaussian_domain_module.md) for
  current RG algorithm authority
- [Residual Gaussian orthogonality robustness](residual_gaussian_orthogonality_robustness.md)
  for the narrow final residual identity-check robustness, tolerance, and
  cutoff-policy lanes
- [Residual Gaussian injection hybrid memo](residual_gaussian_injection_hybrid.md)
  for the optional near-gausslet injection proposal, protected-original
  compact-main design, and the narrow staged geometry source prototype
- [Reference-density-matrix IDA correction](rho0_reference_density_matrix.md)
  for the successor to row-gauge rho0/Galerkin probes: fixed `P0`, exact and
  approximate reference Fock/energy matching, and the measurement-only
  `HP-RHO0-REFDENS-AUDIT-01` lane
- [Terminal shellification due diligence](terminal_shellification_due_diligence.md)
  for the derived basis review report and advisory warning contract
- [PQS complete-shell aspect source modes](pqs_complete_shell_aspect_source_modes.md)
  for the separate `(q,q,L)` source-policy lane
- [Cartesian Gaussian raw blocks - nuclear slice](cartesian_gaussian_raw_blocks_nuclear.md)
  for the neutral uncharged nuclear raw-block owner
- [Cartesian Gaussian raw blocks - non-nuclear slice](cartesian_gaussian_raw_blocks_non_nuclear.md)
  for the neutral overlap/kinetic/moment raw-block owner
- [R3 terminal G-G product matrices](r3_terminal_gg_product_matrices.md)
  for the narrow final-basis `G-G` product-matrix optimization lane
- [R3 remaining exact-operator allocation audit](r3_remaining_exact_operator_allocation_audit.md)
  for the measurement-only decision after terminal `G-G` workspace reuse
- [R3 unit-nuclear U_GG Gaussian sum](r3_unit_nuclear_ugg_gaussian_sum.md)
  for the narrow terminal final-basis unit-nuclear `U_GG` optimization lane
- [PQS/WL mapping `s_factor`](pqs_mapping_s_factor.md)
  for the expert mapping-strength scalar that preserves default behavior while
  allowing CR2-style scans
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
- [Nesting/supplement composition plan](nesting_supplement_composition_plan.md)
  for the target 2 x 2 x 2 composition matrix and dependency order
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
- [Route inventory type-surface cleanup](route_inventory_type_surface_cleanup.md)
  for the first retained-unit route inventory cleanup lane
- [Raw product source-mode inventory cleanup](raw_product_source_mode_inventory_cleanup.md)
  for the vector-backed raw product source-mode inventory cleanup lane
- [Contract-plan vector cleanup](contract_plan_vector_cleanup.md)
  for the terminal-lowering and retained-unit transform contract-plan vector
  cleanup lane
- [Route/stage type-surface cleanup](route_stage_type_surface_cleanup.md)
  for the Be2 q5 compile-attributed route/stage compatibility-inventory
  cleanup lane
- [Route/stage carrier cleanup](route_stage_carrier_cleanup.md)
  for the post-cleanup route/stage carrier and plan-signature cleanup lane
- [Complete-core-shell RHF retirement](complete_core_shell_rhf_retirement.md)
  for the narrow deletion lane for the stale RHF payload stack
- [Route-driver materialization workflow retirement](route_driver_materialization_retirement.md)
  for retiring the old route-driver materialization/report/save wrapper
  workflow and stale tool/test pressure
- [Algorithm implementation index](../../algorithm_implementation_index.md)

Approved amendments:

- [R1 public base producer](r1_public_base_producer.md) defines the approved
  minimal public base Hamiltonian producer surface for first origin-centered H
  and z-axis H2 implementation.
- [R1 one-center base atoms](r1_one_center_base_atoms.md)
  relaxes the one-center base producer from H-only to explicit origin-centered
  all-electron atoms, while requiring atoms and diatomics to share the same
  producer workflow after geometry/shellification normalization.
- [R3 residual-GTO/MWG augmentation](r3_residual_gto_mwg_augmentation.md)
  records implementation history and compact artifact provenance for the first
  H2 endpoint. Current residual Gaussian algorithm authority is the domain
  module page below.
- [R3 usability supplemented workflow](r3_usability_supplemented_workflow.md)
  approves only a non-exported supported facade that constructs H2 and
  internal/performance-supported Be2 supplemented artifacts from system, base
  basis, and supplement specs.
- [Residual Gaussian domain module](residual_gaussian_domain_module.md)
  is the canonical current RG algorithm contract for residual basis selection,
  exact augmented operators, matched-width Gaussian descriptors, and residual
  IDA interaction blocks.
- [Residual Gaussian orthogonality robustness](residual_gaussian_orthogonality_robustness.md)
  approves only a narrow robust final residual identity validation pass for
  small floating-point overshoots with healthy owner-local and merge metrics,
  plus the default final identity tolerance and residual occupation cutoff
  policy.
- [Cartesian Gaussian raw blocks - nuclear slice](cartesian_gaussian_raw_blocks_nuclear.md)
  approves only the neutral owner for exact uncharged by-center Gaussian
  nuclear `G-A`/`A-A` raw blocks shared by Residual Gaussian and Qiu-White
  consumers.
- [Cartesian Gaussian raw blocks - non-nuclear slice](cartesian_gaussian_raw_blocks_non_nuclear.md)
  approves only the neutral owner for exact non-nuclear Gaussian
  overlap/kinetic/coordinate-moment/second-moment `G-A`/`A-A` raw blocks.
- [R3 terminal G-G product matrices](r3_terminal_gg_product_matrices.md)
  approves only the terminal final-basis `G-G` product matrices used by
  residual-Gaussian exact augmented operators.
- [R3 remaining exact-operator allocation audit](r3_remaining_exact_operator_allocation_audit.md)
  approves only a measurement-only audit of the post-`954c86cd` Cr2 exact
  augmented-operator allocation. It does not approve a unit-nuclear, route
  setup, raw-block setup, or other source lane.
- [R3 unit-nuclear U_GG Gaussian sum](r3_unit_nuclear_ugg_gaussian_sum.md)
  approves only terminal final-basis unit-nuclear `U_GG` Gaussian-sum
  allocation reduction under `CartesianFinalBasisRealization`.
- [Cartesian driver usability workflow](cartesian_driver_usability_workflow.md)
  approves only the compact canonical `bin/cartesian_ham_builder.jl` workflow:
  visible defaults, optional trusted input file, command-line overrides,
  visible public contract construction, the public `nesting = :pqs` / `:wl`
  construction-family input, compact run-level hooks, visible physics-level
  construction stages, coarse timing/summary, and artifact production through
  approved producer surfaces. The staged producer surface remains
  driver-facing through `src/cartesian_base_hamiltonian.jl`; lower operator
  owners may be factored only to expose coarse product/moment, unit-nuclear,
  and electron-electron stage timings. This is not route-stage/report
  authority.
- [White-Lindsey terminal basis realization](white_lindsey_terminal_basis_realization.md)
  approves only the narrow seam that lets the existing
  `:white_lindsey_low_order` route produce the same
  `CartesianTerminalBasisRealization` consumed by the staged Hamiltonian path,
  and records the separate WL diatomic compact retained-basis correction.
- [Public ns direct-core side parity](public_ns_core_side_parity.md)
  approves only deriving the direct nucleus-centered core side from public
  `ns`, with oddization limited to direct core identity blocks and no change to
  boundary retained-count construction.
- [Common terminal shell decomposition](common_terminal_shell_decomposition.md)
  approves only route-family-free common shell/core region decomposition and a
  narrow audit/cleanup surface before PQS/WL retained-realization geometry.
- [Mapped-COMX source span](mapped_comx_source_span.md)
  approves only the mainline doside source-span option, compact provenance,
  narrow wiring, terminal-basis consumption of materialized axis facts, compact
  driver/facade selection, and bounded real-atom validation gates. It does not
  import high-order scaffolding, add a CRPS numerical builder, change source
  defaults, or change artifacts.
- [Complete-core-shell RHF retirement](complete_core_shell_rhf_retirement.md)
  approves only deleting the stale `pqs_multilayer_complete_core_shell_rhf.jl`
  stack and root include, with no replacements, adapters, new status/payload
  objects, or workflow changes.
- [Route-driver materialization workflow retirement](route_driver_materialization_retirement.md)
  approves only retiring the old `cartesian_materialization`,
  `cartesian_print_summary`, `cartesian_print_details`, and `cartesian_save`
  wrapper workflow plus stale tool/test/docs pressure. It does not change the
  canonical driver, staged producer functions, artifacts, or numerical kernels.
- [Cartesian driver atom workflow](cartesian_driver_atom_workflow.md)
  approves only explicit origin-centered one-center base atom driver inputs
  through the existing base facade. Current validation remains the
  origin-centered H endpoint for that base-only driver lane. Supplemented atoms
  are governed separately by `HP-COMP-SUPPATOM-*`.
- [R3 homonuclear z-axis diatomic supplemented workflow](r3_homonuclear_diatomic_supplemented_workflow.md)
  approves explicit homonuclear two-center z-axis diatomic supplemented inputs
  through the existing facade and canonical driver, with no element-specific
  defaults or Cr2-specific branch.
- [Nesting/supplement composition plan](nesting_supplement_composition_plan.md)
  records the target three-choice contract:
  `geometry = atom | z-axis diatomic`, `nesting = :pqs | :wl`, and
  `supplement = off | on`. The WL z-axis diatomic base cell is approved under
  `HP-COMP-WLDIAT-*`, and the base homonuclear z-axis diatomic validation
  relaxation is approved under `HP-COMP-BASEDIAT-*`. Supplemented WL z-axis
  diatomics are approved under `HP-COMP-SUPPWL-*`; supplemented one-center
  atoms are approved under `HP-COMP-SUPPATOM-*`. WL diatomic `ns` early
  rejection and retained-support saturation wording are approved under
  `HP-COMP-WLNS-*`; the compact WL diatomic retained-basis correction is
  approved under `HP-WLDIAT-COMPACT-*`.
- [Cartesian Hamiltonian artifact manifest](cartesian_hamiltonian_artifact_manifest.md)
  approves only compact JLD2 sidecar groups for matrix-order basis labels and
  uniform recipe provenance, and separately approves a compact source-mode
  provenance seam for optional source-shell/source-mode manifest groups. It
  does not change the Hamiltonian object, matrix keys, public inputs, or
  `read_cartesian_ida_hamiltonian`.
- [Route inventory type-surface cleanup](route_inventory_type_surface_cleanup.md)
  approves only replacing runtime-keyed retained-unit and pair-family
  `NamedTuple` route inventories in `src/pqs_source_box_route_driver_helpers.jl`
  with vector/table or stable-dictionary storage. Raw product source-mode tuple
  inventories were deferred from that lane.
- [Raw product source-mode inventory cleanup](raw_product_source_mode_inventory_cleanup.md)
  approves replacing `RawProductBoxPlan` source-mode and source-mode-column
  tuple inventories with vector-backed storage while preserving deterministic
  source-mode ordering, retained-rule behavior, and manifest source-mode
  provenance.
- [Contract-plan vector cleanup](contract_plan_vector_cleanup.md)
  approves replacing terminal-lowering and retained-unit transform plan-level
  contract tuple inventories with vector-backed storage while preserving
  accessors, iteration order, summaries, and behavior. Per-contract
  `source_cpbs` tuples remain out of scope.
- [Route/stage type-surface cleanup](route_stage_type_surface_cleanup.md)
  approves a narrow cleanup of Be2 q5 compile-attributed route/stage
  compatibility inventories in `src/pqs_source_box_route_driver_helpers.jl`
  and `src/cartesian_terminal_shellification_geometry.jl`. It does not approve
  driver changes, artifact changes, numerical changes, route diagnostics, or
  precompile/sysimage work.
- [Route/stage carrier cleanup](route_stage_carrier_cleanup.md)
  approves a narrow follow-up cleanup of stage-carried shellification,
  route-skeleton, support-plan, retained-rule-plan, and terminal-plan carriers
  in the route helper and complete-core-shell path. It preserves driver,
  artifact, route semantic, and numerical behavior.

Candidate amendments:

- Translated atoms, Cr2-specific workflow, public supplemented
  workflow/export, basis/supplement-realism beyond explicit supplied
  labels/files, and broad driver diagnostics remain candidate-only until
  separately approved. Supplemented atoms and supplemented WL z-axis diatomics
  are approved under `HP-COMP-SUPPATOM-*` and `HP-COMP-SUPPWL-*`.

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
