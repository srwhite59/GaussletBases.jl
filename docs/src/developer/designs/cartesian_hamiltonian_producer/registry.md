# Registry

Source work requires an approved or implemented source-bearing entry, an
explicit source-maintenance permission or equivalent legacy grant, and the
active `AGENTS.md` whitelist. A completed validation contract remains active
only for the exact test/validation maintenance its entry explicitly grants;
completed evidence without that permission authorizes no work. Measurement-
only entries do not authorize production source edits. Candidate or rejected
entries do not authorize implementation.

This registry is the authority lookup, not the algorithm manual. An entry's
explicit permission and lifecycle govern when present; surrounding section
headings are navigation only. During this transition, a legacy entry without
an explicit `Status:` remains source-bearing only when it is also on the
active `AGENTS.md` whitelist and its prose grants that exact surface.
Permission and lifecycle are distinct: an entry may be
measurement-only, design-only, validation-only, or source-bearing, and may be
approved, implemented, superseded, retired, suspended, or rejected. Closed,
superseded, rejected, and completed-retirement entries remain here as records
without remaining on the active source whitelist in `AGENTS.md`; a historical
implementation may remain preservation-only only when its entry says so and it
is still whitelisted. Numerical formulas, behavioral invariants, and rationale
belong in the linked canonical subsystem document. Every ID has its own
heading plus explicit lifecycle, permission, and canonical or historical
ownership. The generated
[registry/whitelist shadow](registry_whitelist_shadow.toml) preserves these
records and raw `AGENTS.md` whitelist membership for exact parity checks. It is
explicitly non-authoritative and authorization-incomplete; this Markdown
registry and `AGENTS.md` remain the authority until a separately reviewed
atomic cutover. Regenerate the shadow only with
`julia --project=docs docs/check_cartesian_authority_shadow.jl --write`.

## Approved And Implemented

### HP-OBJ-01 — `CartesianTerminalBasisBlock`

Lifecycle: implemented. Permission: source maintenance.

Owner/canonical: `CartesianFinalBasisRealization`;
[terminal basis and base assembly](terminal_basis_and_base_assembly.md).

Source: `src/cartesian_final_basis_realization/pqs_terminal_basis_realization.jl`.

Validation: `test/driver_public/cartesian_base_hamiltonian_runtests.jl` and
downstream `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl`.

Dependencies: authoritative terminal support and retained records.

Permission: preserve the exact support-local block object and direct-identity
versus compact-coefficient semantics.

Non-goals: global coefficients/overlap, route metadata, reports, artifacts, or
support growth.

### HP-OBJ-02 — `CartesianTerminalBasisRealization`

Lifecycle: implemented. Permission: source maintenance.

Owner/canonical: `CartesianFinalBasisRealization`;
[terminal basis and base assembly](terminal_basis_and_base_assembly.md).

Source: `src/cartesian_final_basis_realization/pqs_terminal_basis_realization.jl`.

Validation: `test/driver_public/cartesian_base_hamiltonian_runtests.jl` and
terminal inventory/due-diligence consumers.

Dependencies: `HP-OBJ-01` blocks and disjoint terminal support ownership.

Permission: preserve the exact realization fields, native block order, final
dimension, and structural-overlap interpretation. `max_cross_overlap` remains
legacy object shape, not a physical repair signal.

Non-goals: global basis matrices, route-stage state, algorithm inputs, or new
provenance fields outside separately approved artifact authority.

### HP-FILE-01 — terminal realization file

Lifecycle: implemented. Permission: source maintenance.

Owner/canonical: `CartesianFinalBasisRealization`;
[terminal basis and base assembly](terminal_basis_and_base_assembly.md).

Source: `src/cartesian_final_basis_realization/pqs_terminal_basis_realization.jl`.

Validation: current base H/H2 and downstream augmented H2 gates.

Dependencies: `HP-OBJ-01`, `HP-OBJ-02`, and typed terminal records.

Permission: maintain the implemented terminal object and PQS realization owner.

Non-goals: new modules, public exports, shell geometry, artifacts, or workflow.

### HP-FN-00 — block-local terminal shell realization

Lifecycle: implemented. Permission: source maintenance.

Owner/canonical: `CartesianFinalBasisRealization`;
[terminal basis and base assembly](terminal_basis_and_base_assembly.md).

Source: `src/cartesian_final_basis_realization/pqs_terminal_basis_realization.jl`.

Validation: current base H/H2 gates; structural correction commits
`d2bf139c6` and `d6968d15b`.

Dependencies: retained source-mode contracts and authoritative shell support.

Permission: maintain support-local shell seed, Gram/Lowdin realization, and
deterministic sign canonicalization.

Non-goals: recursive projection, support growth, global Lowdin, or retained
policy changes.

### HP-FN-01 — terminal basis realizer

Lifecycle: implemented. Permission: source maintenance.

Owner/canonical: `CartesianFinalBasisRealization`;
[terminal basis and base assembly](terminal_basis_and_base_assembly.md).

Source: `src/cartesian_final_basis_realization/pqs_terminal_basis_realization.jl`.

Validation: `test/driver_public/cartesian_base_hamiltonian_runtests.jl` and
downstream augmented H2 validation.

Dependencies: `HP-FN-00`, terminal support/retained records, transform
contracts, and mapped axis bundles.

Permission: maintain the live `pqs_terminal_basis_realization(...)` signature
and direct/shell/compact-slab dispatch recorded in the canonical contract.

Non-goals: removed `cross_atol` plumbing, global overlap repair, route
reconstruction, or public API.

### HP-FN-02 — structural terminal support checks

Lifecycle: implemented. Permission: source maintenance.

Owner/canonical: `CartesianFinalBasisRealization`;
[terminal basis and base assembly](terminal_basis_and_base_assembly.md).

Source: `src/cartesian_final_basis_realization/pqs_terminal_basis_realization.jl`.

Validation: structural support checks in the live realizer and base/terminal
consumer gates; implementation commit `d6968d15b`.

Dependencies: orthonormal parent rows and authoritative disjoint supports.

Permission: maintain exact support equality, duplicate-row rejection,
pairwise disjointness, and shell-local identity validation.

Non-goals: computed cross-block residuals, previous-block projection, or
removal of the legacy object field without separate authority.

### HP-WIRE-01 — terminal-basis stage integration

Lifecycle: implemented. Permission: source maintenance.

Owner/canonical: `cartesian_transforms`;
[terminal basis and base assembly](terminal_basis_and_base_assembly.md).

Source: `src/pqs_source_box_route_driver_helpers.jl`.

Validation: `test/driver_public/cartesian_base_hamiltonian_runtests.jl`.

Dependencies: typed terminal plans and `HP-FN-01`.

Permission: pass current PQS support, retained, transform, and bundle objects
directly into the terminal realizer.

Non-goals: report reconstruction, new stage fields, WL policy, artifacts, or
public workflow changes.

## Implemented White-Lindsey Terminal Basis

Canonical contract:
[White-Lindsey terminal basis realization](white_lindsey_terminal_basis_realization.md).

### HP-WLTERM-FILE-01 — WL terminal realization file

Lifecycle: implemented. Permission: source maintenance.

Canonical contract:
[White-Lindsey terminal basis realization](white_lindsey_terminal_basis_realization.md).

Owner/source:
`src/cartesian_final_basis_realization/white_lindsey_terminal_basis_realization.jl`
and its include in `CartesianFinalBasisRealization.jl`.

Permission: maintain the existing WL-specific sibling that returns the shared
`CartesianTerminalBasisRealization`. No new module, basis object, route result,
artifact, report, or export is authorized.

### HP-WLTERM-FN-01 — WL low-order terminal basis realization

Lifecycle: implemented. Permission: source maintenance.

Owner/canonical: White-Lindsey terminal basis realization;
[White-Lindsey terminal basis realization](white_lindsey_terminal_basis_realization.md).

Source:

- `src/cartesian_final_basis_realization/pqs_terminal_basis_realization.jl`;
- `src/cartesian_final_basis_realization/white_lindsey_terminal_basis_realization.jl`.

Dependency: the neutral face-product coefficient helper remains separately
owned by `HP-COMP-FACEPROD-*`; it is not an additional source surface under
this ID.

Permission: realize direct identity blocks and compact WL facet/edge/corner,
boundary-stratum, and thin-slab products on authoritative owned supports while
preserving retained/transform order and block-local identity checks.

Non-goals: old WL H1/H1+J materialization, route or shell policy, artifacts,
public API, raw blocks, RG/MWG/IDA, solvers, or Cr2 workflow.

### HP-WLTERM-WIRE-01 — WL route helper terminal-basis wiring

Lifecycle: implemented. Permission: source maintenance.

Canonical contract:
[White-Lindsey terminal basis realization](white_lindsey_terminal_basis_realization.md).

Owner/source: `src/pqs_source_box_route_driver_helpers.jl`.

Permission: pass native `:white_lindsey_low_order` support, retained, and
transform records into the WL realizer without changing PQS behavior or
restoring old WL materialization.

### HP-WLTERM-TEST-01 — WL terminal-basis validation

Lifecycle: completed validation contract. Permission: validation maintenance.

Owner/canonical: White-Lindsey terminal basis realization;
[White-Lindsey terminal basis realization](white_lindsey_terminal_basis_realization.md).

Evidence: accepted bounded PQS parity, WL atom/H2 artifact/readback, and
supplemented endpoint smokes recorded in the manager log. No dedicated
committed WL fixture is owned by this ID.

### HP-COMP-WLDIAT-FN-01 — WL z-axis diatomic base terminal records

Lifecycle: implemented. Permission: source maintenance.

Owner/canonical: White-Lindsey terminal basis realization;
[White-Lindsey terminal basis realization](white_lindsey_terminal_basis_realization.md).

Source:

- `src/pqs_source_box_diatomic_complete_core_shell.jl`;
- `src/cartesian_terminal_shellification_geometry.jl`;
- `src/cartesian_terminal_lowering/selection.jl`;
- `src/cartesian_terminal_lowering/region_contracts.jl`;
- `src/pqs_source_box_route_driver_helpers.jl`;
- `src/cartesian_final_basis_realization/pqs_terminal_basis_realization.jl`;
- `src/cartesian_final_basis_realization/white_lindsey_terminal_basis_realization.jl`;
- `src/cartesian_final_basis_realization/CartesianFinalBasisRealization.jl`;
- `src/cartesian_base_hamiltonian.jl`.

Permission: maintain native WL z-axis diatomic terminal records and the shared
base Hamiltonian path, including truthful route provenance
`:z_axis_diatomic_wl_base`.

Non-goals: driver special cases, parallel Hamiltonian construction, old WL
materialization, artifacts/readers, RG/MWG/IDA policy, solver/ECP, or Cr2.

### HP-COMP-WLDIAT-TEST-01 — WL z-axis diatomic base validation

Lifecycle: completed validation contract. Permission: validation maintenance.

Owner/canonical: White-Lindsey terminal basis realization;
[White-Lindsey terminal basis realization](white_lindsey_terminal_basis_realization.md).

Evidence: bounded PQS/WL atom and H2 artifact/readback, route-provenance, and
downstream supplemented smokes recorded in the manager log.

### HP-WLDIAT-COMPACT-FN-01 — WL diatomic compact retained basis

Lifecycle: implemented. Permission: source maintenance.

Canonical contract:
[White-Lindsey terminal basis realization](white_lindsey_terminal_basis_realization.md).

Owner/source:

- `src/cartesian_shellification/terminal_geometry.jl`;
- `src/cartesian_terminal_lowering/region_contracts.jl`;
- `src/cartesian_retained_units/lower_contract_units.jl`;
- `src/cartesian_retained_unit_transform_contracts/unit_contracts.jl`;
- `src/cartesian_final_basis_realization/white_lindsey_terminal_basis_realization.jl`;
- `src/pqs_source_box_route_driver_helpers.jl`, for narrow route wiring only.

Dependency: the neutral face-product helper remains separately governed by
`HP-COMP-FACEPROD-*`; it is not an additional source surface under this ID.

Permission: maintain compact products of one-dimensional contractions for WL
boundary units. Identity realization remains valid only for true direct/core
units; support rows are not themselves retained functions.

### HP-WLDIAT-COMPACT-TEST-01 — WL compact-basis validation

Lifecycle: completed validation contract. Permission: validation maintenance.

Owner/canonical: White-Lindsey terminal basis realization;
[White-Lindsey terminal basis realization](white_lindsey_terminal_basis_realization.md).

Evidence: bounded H2/Be2 WL base and supplemented artifact/readback,
finite/symmetric operators, retained-count, and PQS parity smokes.

### HP-WLDIAT-PARITY-FN-01 — WL boundary retained-count parity

Lifecycle: implemented. Permission: source maintenance.

Canonical contract:
[White-Lindsey terminal basis realization](white_lindsey_terminal_basis_realization.md).

Owner/source:
`src/cartesian_final_basis_realization/white_lindsey_terminal_basis_realization.jl`.

Permission: keep odd-side enforcement only for direct nucleus-centered cores;
WL boundary products retain the requested count. Canonical examples remain
`ns=4 -> 56` and `ns=5 -> 98` boundary columns.

### HP-WLDIAT-PARITY-TEST-01 — WL parity validation

Lifecycle: completed validation contract. Permission: validation maintenance.

Owner/canonical: White-Lindsey terminal basis realization;
[White-Lindsey terminal basis realization](white_lindsey_terminal_basis_realization.md).

Evidence: bounded `ns=4/5` retained-count, artifact/readback,
finite/symmetric-operator, supplemented, and PQS endpoint smokes.

Family-wide non-goals: public input changes, shell ownership changes, direct
core changes, artifact/schema changes, old WL materialization, route reports,
raw-block or RG/MWG/IDA changes, solvers/ECP, and Cr2 production claims.

## Composition Input: Public ns Direct-Core Side

Canonical contract:
[Public ns direct-core side parity](public_ns_core_side_parity.md).

### HP-COMP-NSCORE-FN-01 — public ns direct-core side

Lifecycle: implemented. Permission: source maintenance.

Owner/canonical: public `ns` direct-core side parity;
[Public ns direct-core side parity](public_ns_core_side_parity.md).

Source:

- `src/cartesian_base_hamiltonian.jl`;
- `src/pqs_source_box_route_driver_helpers.jl`.

Dependencies: `HP-COMP-NS-FN-01` public-size normalization and the existing
direct-core route construction.

Permission: preserve
`direct_core_side = isodd(ns) ? ns : ns + 1` for direct nucleus-centered
identity blocks only. Boundary retained construction remains route-local.

Evidence: implementation commit `e41aba6eb`; manager running-log Pass 160.

Non-goals: driver/public-input changes, shellification or terminal-lowering
redesign, boundary oddization, artifact/manifest changes, old WL
materialization, fixtures, or Cr2 workflow.

### HP-COMP-NSCORE-TEST-01 — public ns direct-core validation

Lifecycle: completed validation contract. Permission: validation maintenance.

Canonical contract:
[Public ns direct-core side parity](public_ns_core_side_parity.md).

Evidence and test surfaces:

- `docs/src/developer/pqs_manager_running_log.md`, Pass 160;
- `test/driver_public/cartesian_base_hamiltonian_runtests.jl`, current
  downstream public-base regression.

The accepted gate covered bounded same-`ns` PQS/WL atom dimensions,
even-`ns` parity, provenance, readback, and a diatomic smoke. No dedicated
ID-specific committed test file or Cr2 fixture was added.

## Common Shellification And Compact Thin-Slab Family

Canonical contract:
[Common terminal shell decomposition](common_terminal_shell_decomposition.md).

This family owns route-neutral first-step shell geometry, native angular-z slab
classification, shared face-product realization, and compact thin-slab
lowering. It does not own PQS aspect-aware complete-shell source dimensions,
which have a separate canonical contract and IDs.

### HP-COMP-SHELLGEOM-FN-01 — common shell decomposition

Status: implemented.

Canonical contract:
[common terminal shell decomposition](common_terminal_shell_decomposition.md).

Owner/source:
`src/cartesian_shellification/terminal_geometry.jl`, with narrow caller
plumbing in `src/pqs_source_box_route_driver_helpers.jl`.

Permission: maintain route-family-free direct-core/shell regions, ordering,
coverage, and owned support before PQS/WL retained construction diverges.

### HP-COMP-SHELLGEOM-TEST-01 — common shell validation

Status: completed validation evidence.

Permission: none. The completed evidence does not authorize test maintenance
or a committed fixture.

Owner/canonical: common shell decomposition;
[common terminal shell decomposition](common_terminal_shell_decomposition.md).

Evidence: bounded same-input PQS/WL shell/support comparisons and downstream
artifact/endpoint smokes recorded in the manager log. No committed fixture is
owned by this ID.

### HP-COMP-SHELLGEOM-DIAT-FN-01 — diatomic common shellifier entry

Status: implemented.

Canonical contract:
[common terminal shell decomposition](common_terminal_shell_decomposition.md).

Owner/source:
`src/cartesian_shellification/terminal_geometry.jl` and narrow
`src/pqs_source_box_route_driver_helpers.jl` caller plumbing.

Permission: feed the same public `ns`, direct-core side, centers, bond axis,
and parent facts into common z-axis diatomic shellification before family
lowering. Central-gap/contact redesign is not approved.

### HP-COMP-SHELLGEOM-DIAT-TEST-01 — diatomic shellifier-entry validation

Status: completed validation evidence.

Permission: none. The completed evidence does not authorize test maintenance
or a committed Cr2 gate.

Owner/canonical: common diatomic shell geometry;
[common terminal shell decomposition](common_terminal_shell_decomposition.md).

Evidence: bounded PQS/WL same-function/same-argument and shell-region parity;
no Cr2 committed gate.

### HP-COMP-OUTERMM-FN-01 — outer-mismatch-only correction

Status: superseded; no permission.

Permission: none.

Owner/history:
[common terminal shell decomposition](common_terminal_shell_decomposition.md).

Superseded by `HP-COMP-THINSLAB-FN-01`. Do not restore a separate
outer-mismatch path.

### HP-COMP-OUTERMM-TEST-01 — outer-mismatch-only validation

Status: superseded; no permission.

Permission: none.

Owner/history:
[common terminal shell decomposition](common_terminal_shell_decomposition.md).

Historical evidence is subsumed by common thin-slab validation.

### HP-COMP-ANGBOX-FN-01 — angular-balanced diatomic shellification

Status: implemented.

Canonical contract:
[common terminal shell decomposition](common_terminal_shell_decomposition.md).

Owner/source:
`src/cartesian_shellification/terminal_geometry.jl`.

Permission: emit native ordered `:angular_z_extension_slab` stacks so the
ordinary shell body plus planned axial extensions realizes the physical
outer-nucleus angular target. It does not change real-shell retained policy or
central-gap/contact ownership.

### HP-COMP-ANGBOX-TEST-01 — angular shellification validation

Status: completed validation evidence.

Permission: none. The completed evidence does not authorize test maintenance
or a committed Cr2 fixture.

Owner/canonical: angular extension shell geometry;
[common terminal shell decomposition](common_terminal_shell_decomposition.md).

Evidence: ignored H2/Be2/Cr2-style geometry inventories and accepted
shellification Passes 179/186; no production Cr2 claim.

### HP-COMP-FACEPROD-FN-01 — neutral terminal face-product helper

Status: implemented.

Canonical contract:
[common terminal shell decomposition](common_terminal_shell_decomposition.md).

Owner/source:

- `src/cartesian_final_basis_realization/terminal_face_product_blocks.jl`;
- module include in
  `src/cartesian_final_basis_realization/CartesianFinalBasisRealization.jl`;
- consumers in PQS and White-Lindsey terminal realization.

Permission: one route-neutral face/face-stack coefficient assembly over fixed
normal-axis indices. It is not a new terminal-basis policy.

### HP-COMP-FACEPROD-TEST-01 — face-product validation

Status: completed validation evidence.

Permission: none. The completed evidence does not authorize test maintenance
or a dedicated committed fixture.

Owner/canonical: neutral terminal face products;
[common terminal shell decomposition](common_terminal_shell_decomposition.md).

Evidence: exact White-Lindsey facet coefficient parity and compact slab reuse;
no committed fixture.

### HP-COMP-THINSLAB-FN-01 — common compact thin-slab lowering

Status: implemented.

Canonical contract:
[common terminal shell decomposition](common_terminal_shell_decomposition.md).

Owner/source:

- `src/cartesian_terminal_lowering/selection.jl`;
- `src/cartesian_terminal_lowering/region_contracts.jl`;
- `src/cartesian_retained_units/lower_contract_units.jl`;
- `src/cartesian_retained_unit_transform_contracts/unit_contracts.jl`;
- PQS/WL terminal realization through the neutral face-product owner;
- native shellification/caller support in the established terminal owners.

Permission: midpoint, outer-mismatch, and angular-z-extension slabs lower as
compact face stacks for both PQS and WL, never as full identity CPBs. Real
shells remain family-specific.

### HP-COMP-THINSLAB-TEST-01 — compact thin-slab validation

Status: completed validation evidence.

Permission: none. The completed evidence does not authorize test maintenance
or a committed Cr2 gate.

Owner/canonical: common compact thin-slab lowering;
[common terminal shell decomposition](common_terminal_shell_decomposition.md).

Evidence: PQS/WL H2/Be2 artifact/readback, compact retained counts, exact WL
facet parity, and H2 residual-GTO endpoint replay; no committed Cr2 gate.

### HP-COMP-THINSLAB-META-FN-01 — thin-slab inventory metadata

Status: implemented.

Canonical contract:
[common terminal shell decomposition](common_terminal_shell_decomposition.md).

Owner/source:
`src/cartesian_terminal_shellification_geometry.jl`.

Permission: describe native compact slab kinds consistently in internal
inventory/scaffold summaries. It does not materialize coefficients or create
artifact/report payloads.

### HP-COMP-THINSLAB-META-TEST-01 — thin-slab inventory validation

Status: completed validation evidence.

Permission: none. The completed evidence does not authorize test maintenance,
coefficient work, or report/artifact payloads.

Owner/canonical: terminal shellification inventory metadata;
[common terminal shell decomposition](common_terminal_shell_decomposition.md).

Evidence: focused inventory checks confirm midpoint, outer-mismatch, and
angular-z-extension slabs no longer advertise direct identity lowering.

Family-wide non-goals: public inputs, artifact schema, route-report
frameworks, RG/MWG/IDA, Hamiltonian semantics, solver/ECP, real-shell
PQS/WL convergence, or Cr2 production claims.
## Mapped-COMX Source-Span Facility

Canonical contract:
[Mapped-COMX source span](mapped_comx_source_span.md).

Ordinary remains the default. The installed mapped option is PQS-only and
remains opt-in after the bounded He `ns=5` physics limitation.

### HP-MCOMX-FILE-01 — mapped-COMX source ownership

Status: implemented.

Permission: source maintenance in the existing owners only; no new production
file or second COMX owner.

Owner/canonical: mapped-COMX source span;
[mapped-COMX source span](mapped_comx_source_span.md).

Source owners:
`src/cartesian_nested_faces.jl` and narrow existing PQS source-axis plumbing.
No new production file was created under `CartesianRawProductSources`.

### HP-MCOMX-OBJ-01 — mapped source specification

Status: implemented internal construction specification.

Owner/canonical: mapped-COMX source span;
[Mapped-COMX source span](mapped_comx_source_span.md).

Permission: fixed protected-P2, mapped-Chebyshev, lambda/no-sqrt-J, and
physical-localization facts. No public export or general tuning object.

### HP-MCOMX-FN-01 — mapped source-span construction

Status: implemented.

Canonical contract:
[mapped-COMX source span](mapped_comx_source_span.md).

Owner/source:
`src/cartesian_nested_faces.jl`.

Permission: construct normalized-local-coordinate mapped enrichment before the
existing physical-coordinate COMX cleanup. Ordinary behavior remains default.

### HP-MCOMX-WIRE-01 — PQS axis-transform wiring

Status: implemented.

Canonical contract:
[mapped-COMX source span](mapped_comx_source_span.md).

Owner/source:
`src/cartesian_pair_block_materialization/pqs_source_axis_transforms.jl`.

Permission: pass the internal source-span choice into the existing doside seam
and return ordinary carried `AxisSourceTransformFact` objects.

### HP-MCOMX-TEST-01 — mapped source validation

Status: completed validation evidence.

Permission: none. The completed evidence does not authorize test maintenance,
a committed fixture, or default-promotion work.

Owner/canonical: mapped-COMX source span;
[mapped-COMX source span](mapped_comx_source_span.md).

Evidence: `ns=5/6/7` rank/protected-span/coordinate/orthogonality checks and
bounded H/He/H2 comparisons. He `ns=5` blocks default promotion.

### HP-MCOMX-TERM-FN-01 — terminal shell-seed consumption

Status: implemented.

Canonical contract:
[mapped-COMX source span](mapped_comx_source_span.md).

Owner/source:

- `src/cartesian_final_basis_realization/pqs_terminal_basis_realization.jl`;
- narrow module import/include support in
  `src/cartesian_final_basis_realization/CartesianFinalBasisRealization.jl`.

Permission: validate and consume materialized carried axis facts as the PQS
shell seed while preserving ordinary fallback, boundary selection, support
restriction, Lowdin, and canonicalization.

### HP-MCOMX-TERM-TEST-01 — terminal seam validation

Status: completed validation evidence.

Permission: none. The completed evidence does not authorize test maintenance
or a new committed test surface.

Owner/canonical: mapped-COMX terminal consumption;
[mapped-COMX source span](mapped_comx_source_span.md).

Evidence: carried-fact coefficient parity/difference checks plus ordinary and
supplemented bounded endpoints; no committed Cr2 gate.

### HP-MCOMX-DRV-FN-01 — canonical source-span selector

Status: implemented.

Canonical contract:
[mapped-COMX source span](mapped_comx_source_span.md).

Owner/source:

- `bin/cartesian_ham_builder.jl`;
- `src/cartesian_base_hamiltonian.jl`;
- narrow `src/pqs_source_box_route_driver_helpers.jl` propagation.

Permission: expose `source_span = :ordinary | :mapped_comx`, default ordinary,
and reject mapped-COMX with White-Lindsey. This is not a diagnostic route
switch.

### HP-MCOMX-DRV-TEST-01 — driver source-span validation

Status: completed validation evidence.

Permission: none. The completed evidence does not authorize test maintenance
or a new committed test surface.

Owner/canonical: mapped-COMX driver selection;
[mapped-COMX source span](mapped_comx_source_span.md).

Evidence: ordinary default parity, mapped PQS artifact/readback/provenance, and
bounded endpoint smokes.

Facility-wide non-goals: default promotion, automatic tuning, another COMX
path, White-Lindsey support, route/report payloads, artifact schema,
Hamiltonian/raw-block/RG/MWG/IDA/EGOI/screening/solver changes, or Cr2
production claims.

## Implemented Composition And Input Family

Canonical contracts:

- [Nesting/supplement composition](nesting_supplement_composition_plan.md);
- [R1 one-center base atoms](r1_one_center_base_atoms.md);
- [Public ns direct-core side parity](public_ns_core_side_parity.md).

### HP-COMP-BASEDIAT-FN-01 — base homonuclear z-axis diatomics

Lifecycle: implemented. Permission: source maintenance.

Owner/canonical: nesting/supplement composition;
[Nesting/supplement composition](nesting_supplement_composition_plan.md).

Source: `src/cartesian_base_hamiltonian.jl`.

Dependencies: the shared base producer and the implemented nesting-specific
terminal-basis paths.

Permission: validate explicit equal-symbol/equal-charge neutral all-electron
diatomics at two finite distinct z-axis centers and send them through the
existing PQS/WL base path.

Evidence: implementation commit `095a89d41`; manager running-log Pass 139.

Non-goals: heteronuclear/general geometry, inferred element data, driver,
supplement, route/shellification, artifact, solver/ECP, or Cr2-specific work.

### HP-COMP-BASEDIAT-TEST-01 — base diatomic validation

Lifecycle: completed validation contract. Permission: validation maintenance.

Owner/canonical: nesting/supplement composition;
[Nesting/supplement composition](nesting_supplement_composition_plan.md).

Evidence and test surfaces:

- `docs/src/developer/pqs_manager_running_log.md`, Pass 139;
- `test/driver_public/cartesian_base_hamiltonian_runtests.jl`, current
  downstream base-facade regression.

The accepted gate covered H2 and bounded Be2 PQS/WL artifact/readback plus
invalid geometry, charge, and electron-count rejection. No dedicated
ID-specific committed test file was added.

### HP-COMP-SUPPWL-FN-01 — supplemented WL z-axis diatomics

Lifecycle: implemented. Permission: source maintenance.

Owner/canonical: nesting/supplement composition;
[Nesting/supplement composition](nesting_supplement_composition_plan.md).

Source owners:

- `src/cartesian_base_hamiltonian.jl`;
- `src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl`, shared
  residual-GTO/MWG consumer.

Dependencies: implemented WL terminal realization and the common supplemented
Hamiltonian path.

Permission: a WL terminal basis enters the same residual-GTO, exact augmented
one-body, residual MWG/IDA, assembly, and artifact path as PQS.

Evidence: implementation commit `4cfb47ace`; manager running-log Pass 141.

Non-goals: a WL-specific supplement algorithm, old WL H1/H1+J revival, driver,
route/lowering, RG/MWG convention, artifact, solver/ECP, or Cr2 changes.

### HP-COMP-SUPPWL-TEST-01 — supplemented WL validation

Lifecycle: completed validation contract. Permission: validation maintenance.

Owner/canonical: nesting/supplement composition;
[Nesting/supplement composition](nesting_supplement_composition_plan.md).

Evidence and test surfaces:

- `docs/src/developer/pqs_manager_running_log.md`, Pass 141;
- `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl`, current
  downstream supplemented-path regression.

The accepted gate covered H2 PQS/WL and bounded Be2 WL artifact/readback,
finite/symmetric operators, and invalid-system rejection. No dedicated
ID-specific committed test file was added.

### HP-COMP-SUPPATOM-FN-01 — supplemented one-center atoms

Lifecycle: implemented. Permission: source maintenance.

Owner/canonical: nesting/supplement composition;
[Nesting/supplement composition](nesting_supplement_composition_plan.md).

Source owners:

- `src/cartesian_base_hamiltonian.jl`;
- `bin/cartesian_ham_builder.jl`;
- `src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl`, shared
  residual-GTO/MWG consumer.

Dependencies: the one-center base contract and common supplemented producer.

Permission: origin-centered all-electron atoms use the existing atomic
supplement loader and the same PQS/WL residual-GTO/MWG path as supported
diatomics.

Evidence: implementation commit `2e9818c90`; manager running-log Pass 142.

Non-goals: atom-only Hamiltonian/materialization paths, translated atoms,
element defaults, new driver controls, route/lowering, RG/MWG convention,
artifact, solver/ECP, or Cr2 changes.

### HP-COMP-SUPPATOM-TEST-01 — supplemented atom validation

Lifecycle: completed validation contract. Permission: validation maintenance.

Owner/canonical: nesting/supplement composition;
[Nesting/supplement composition](nesting_supplement_composition_plan.md).

Evidence and test surfaces:

- `docs/src/developer/pqs_manager_running_log.md`, Pass 142;
- `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl`, current
  downstream supplemented-path regression.

The accepted gate covered base H, supplemented H under PQS/WL, bounded Be,
H2 regression smokes, translated-atom rejection, and supplement-count
rejection. No dedicated ID-specific committed test file was added.

### HP-COMP-ATOMBOX-FN-01 — one-center atom physical extent

Lifecycle: implemented. Permission: source maintenance.

Source: `src/cartesian_base_hamiltonian.jl`.

Canonical owner:
[R1 one-center base atoms](r1_one_center_base_atoms.md).

Dependencies: White-Lindsey atomic mapping, public `radius`, spacing policy,
and public-`ns` direct-core minimum.

Permission: derive atom parent counts from physical `basis.radius` and the
existing mapping/spacing policy; `ns` remains resolution/nesting input.

Evidence: implementation commit `18d683575`; manager running-log Pass 146.

Non-goals: driver or diatomic sizing changes, translated atoms, broad parent
redesign, route/RG/MWG/IDA, artifact, solver/ECP, or Cr2 changes.

### HP-COMP-ATOMBOX-TEST-01 — atom physical-extent validation

Lifecycle: completed validation contract. Permission: validation maintenance.

Owner/canonical: one-center physical parent extent;
[R1 one-center base atoms](r1_one_center_base_atoms.md).

Evidence and test surfaces:

- `docs/src/developer/pqs_manager_running_log.md`, Pass 146;
- `test/driver_public/cartesian_base_hamiltonian_runtests.jl`, current
  downstream one-center base regression.

The accepted gate covered atom base/supplemented PQS/WL construction, bounded
Be, radius sensitivity, and H2/Be2 regression smokes. No dedicated
ID-specific committed test file was added.

### HP-COMP-NS-FN-01 — public ns and derived q

Lifecycle: implemented. Permission: source maintenance.

Owner/canonical: nesting/supplement composition;
[Nesting/supplement composition](nesting_supplement_composition_plan.md).

Source:

- `src/cartesian_base_hamiltonian.jl`;
- `bin/cartesian_ham_builder.jl`.

Dependencies: the supported PQS/WL composition matrix.

Permission: normalize public `ns`, derive `q = ns` for PQS and
`q = ns - 2` for WL, validate any legacy `q` compatibility, and record the
existing compact provenance.

Evidence: implementation commit `5f0185a16`; manager running-log Pass 145.

Non-goals: changing physical extent, route construction, shellification,
terminal lowering, artifact format/readback, driver stages, solver/ECP, or Cr2
workflow.

### HP-COMP-NS-TEST-01 — public ns validation

Lifecycle: completed validation contract. Permission: validation maintenance.

Owner/canonical: nesting/supplement composition;
[Nesting/supplement composition](nesting_supplement_composition_plan.md).

Evidence and test surfaces:

- `docs/src/developer/pqs_manager_running_log.md`, Pass 145;
- `test/driver_public/cartesian_base_hamiltonian_runtests.jl`, current
  downstream public-input/provenance regression.

The accepted gate covered atom/diatomic PQS/WL base construction, bounded
supplemented smokes, legacy-`q` compatibility, inconsistent input rejection,
and provenance. No dedicated ID-specific committed test file was added.

### HP-COMP-WLNS-FN-01 — WL diatomic ns guard

Lifecycle: implemented. Permission: source maintenance.

Owner/canonical: nesting/supplement composition;
[Nesting/supplement composition](nesting_supplement_composition_plan.md).

Source: `src/cartesian_base_hamiltonian.jl`.

Dependencies: `HP-COMP-NS-FN-01` normalization and WL complete-shell
diatomic construction.

Permission: reject normalized WL z-axis diatomic `ns < 4` before route
construction and preserve retained-support saturation as valid behavior.

Evidence: implementation commit `50327dc1e`; manager running-log Pass 148.

Non-goals: changing working WL retained support, driver compatibility,
shellification/lowering, artifacts, RG/MWG/IDA, solver/ECP, or Cr2 workflow.

### HP-COMP-WLNS-TEST-01 — WL diatomic ns validation

Lifecycle: completed validation contract. Permission: validation maintenance.

Owner/canonical: nesting/supplement composition;
[Nesting/supplement composition](nesting_supplement_composition_plan.md).

Evidence and test surfaces:

- `docs/src/developer/pqs_manager_running_log.md`, Pass 148;
- `test/driver_public/cartesian_base_hamiltonian_runtests.jl`, current
  downstream WL base regression.

The accepted gate covered early `ns = 3` rejection, working `ns = 4` base
and supplemented smokes, and PQS regressions. No dedicated ID-specific
committed test file was added.

## Approved Base Assembly Completion

### HP-FN-03 — blockwise one-body assembly

Lifecycle: implemented. Permission: source maintenance.

Owner/canonical: `CartesianFinalBasisRealization`;
[terminal basis and base assembly](terminal_basis_and_base_assembly.md).

Source: `src/cartesian_final_basis_realization/pqs_terminal_one_body.jl`.

Validation: `test/driver_public/cartesian_base_hamiltonian_runtests.jl` and
`test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl`.

Dependencies: `HP-OBJ-01`, `HP-OBJ-02`, reusable one-dimensional product and
Gaussian factor data.

Permission: maintain block-pair product assembly and the file-local term-first
Gaussian-sum accumulator within the current workspace bound.

Non-goals: K/U payloads, persistent caches, reports, stage fields, or
one-body orchestration APIs.

### HP-FN-04 — localized IDA assembly

Lifecycle: implemented. Permission: source maintenance.

Owner/canonical: `CartesianFinalBasisRealization`;
[terminal basis and base assembly](terminal_basis_and_base_assembly.md).

Source: `src/cartesian_final_basis_realization/pqs_terminal_ida.jl`.

Validation: `test/driver_public/cartesian_base_hamiltonian_runtests.jl` and
downstream augmented H2 validation; implementation commit `a33842fb8`.

Dependencies: terminal blocks, one-dimensional raw pair terms, positive final
IDA weights, and the producer-wide Coulomb expansion.

Permission: maintain blockwise final-weight-normalized localized IDA assembly
and symmetry validation.

Non-goals: Hamiltonian construction, interaction rotation, four-index tensors,
artifacts, route wiring, or pair payload/cache objects.

## Reference-Density Measurements And Internal Facilities

### HP-PQS-SCREEN-HARTREE-AUDIT-01 - protected-GTO screened Hartree residual-density audit

Status: completed historical H/Be/Be2 measurement; not active source
authority.

Permission: none.

Owner: Cartesian Hamiltonian producer screened-Hartree formalism.

Canonical contract:

- [Screened Hartree residual-density formalism](screened_hartree_residual_density.md)

Committed source/test surfaces: none. The completed work used ignored
measurement probes; evidence is retained in manager Passes 317-319.

Dependencies: exact/Galerkin nuclear attraction, represented occupied
determinants, Galerkin reference Hartree fields, and same-basis IDA/MWG.

Exclusions: source or workflow authority, artifacts, solver integration,
exchange, EGOI, row-gauge substitutions, and Cr2 production claims.

### HP-PQS-SCREEN-HARTREE-NE-AUDIT-01 - Ne screened Hartree endpoint measurement

Status: historical and operationally superseded by the fitted-cloud,
packet-driven practical path. The direct occupied-pair construction remains an
oracle, not a normal workflow.

Permission: none.

Owner: Cartesian Hamiltonian producer screened-Hartree formalism.

Canonical contract:

- [Screened Hartree residual-density formalism](screened_hartree_residual_density.md)

Committed source/test surfaces: none. Evidence and the superseding endpoint
record are retained in manager Passes 320 and 326.

Dependencies: the completed small-system audit and the Ne closed-shell
pure-GTO reference determinant.

Exclusions: active measurement authority, source/workflow changes, broad
first-row claims, exchange, EGOI, and Cr/Cr2 production claims.

### HP-PQS-SCREEN-HARTREE-NE-FITCLOUD-AUDIT-01 - Ne fitted-cloud measurement

Status: completed historical measurement; not active source authority.

Permission: none.

Owner: Cartesian Hamiltonian producer screened-Hartree formalism.

Canonical contract:

- [Screened Hartree residual-density formalism](screened_hartree_residual_density.md)

Committed source/test surfaces: none. Evidence is retained in manager
Passes 321-326.

Dependencies: the represented Ne occupied determinant and the
determinant-density oracle.

Exclusions: a tunable model cloud, fit Gaussians as protected orbitals, public
workflow, artifacts, solver integration, exchange, EGOI, and Cr/Cr2 claims.

### HP-PQS-SCREEN-HARTREE-POTFIT-AUDIT-01 - fitted-potential measurement

Status: completed historical measurement; its durable packet semantics are now
implemented under the atomic reference packet contract.

Permission: none.

Owner: historical screened-Hartree measurement; durable fit semantics are
owned by the `CartesianReferenceDensity` atomic packet subsystem.

Canonical contracts:

- [Screened Hartree residual-density formalism](screened_hartree_residual_density.md)
- [Atomic HF reference packets](atomic_hf_reference_packets.md)

Committed source/test surfaces: none for this historical audit. Evidence is
retained in manager Passes 327-329; packet source ownership begins with
Passes 330-331.

Dependencies: a determinant-defined reference, its validated density fit, and
the density-fit Galerkin `J0_G` oracle.

Exclusions: fitted potential terms as a density or protected orbitals, source
authority from this audit ID, public workflow, corrected artifacts, solver
integration, exchange, EGOI, and Cr/Cr2 production claims.

### HP-PQS-ATOMREF-PACKET-FN-01 — reusable atomic HF reference packets

Lifecycle: implemented. Permission: source maintenance.

Owner/canonical: `CartesianReferenceDensity`;
[atomic HF reference packets](atomic_hf_reference_packets.md).

Source:

- `src/cartesian_reference_density/CartesianReferenceDensity.jl`;
- `src/cartesian_reference_density/atomic_hf_reference_packets.jl`;
- include/qualified wiring in `src/GaussletBases.jl`;
- narrow packet-consumer rejection in
  `src/cartesian_reference_density/screened_hartree_correction.jl`;
- existing mixed-Hartree helpers only as canonical packet evaluation needs.

Permission: maintain converged one-center determinant packets, exact packet
self-integrity, exact owner/order/placement mapping, numerical owner-local
overlap equivalence at `1e-10`, ordinary density and radial-potential fits,
read/write validation, and explicit fit/provenance diagnostics. Density fits
own `E0`; potential fits approximate `J0`. Polished legacy packets reject.

Non-goals: corrected artifacts, public defaults, solvers, exchange, EGOI,
row-gauge rho0/P0, Cr2 claims, inferred occupancy, or fitted terms as orbitals.

### HP-PQS-ATOMREF-PACKET-TEST-01 — atomic packet validation

Lifecycle: completed validation contract. Permission: validation maintenance.

Owner/canonical: `CartesianReferenceDensity`;
[atomic HF reference packets](atomic_hf_reference_packets.md).

Tests:

- `test/nested/cartesian_atomic_hf_reference_packet_runtests.jl`;
- narrow consumer/embedding checks in
  `test/nested/cartesian_screened_hartree_correction_runtests.jl`;
- vendored-basis identity/parser checks in `test/misc/runtests.jl`.

Permission: maintain packet roundtrip, convergence, determinant/fit,
fingerprint, embedding-equivalence, compact-Coulomb-role, and malformed-input
coverage. The scientific body of `data/legacy/BasisSets` is not test authority
to rewrite that data.

### HP-PQS-ATOMREF-POTMOM-FN-01 - retired determinant-moment polish

Lifecycle: retired. Permission: none.

Owner/history:
[atomic HF reference packets](atomic_hf_reference_packets.md).

This same-day false start was approved in commit `9739c22a6` and implemented
as the moment-polish portion of manager Pass 353. It adjusted fitted-potential
coefficients against determinant moments on a fixed separation grid to force
one padded Be2 energy-consistency value below `1e-8 Ha`. That is not a generic
atomic or molecular fitting principle and must not be retained through an
adapter.

Durable packet, screened-Hartree, and additive-reference behavior is owned by
their active IDs; the ordinary radial fit supersedes this source behavior. No
adapter, molecule-trained fit, or polished-packet consumption is authorized.

### HP-PQS-ATOMREF-POTMOM-TEST-01 - retired polish validation

Lifecycle: retired historical evidence. Permission: none.

Owner/history:
[atomic HF reference packets](atomic_hf_reference_packets.md).

Evidence: manager Passes 351 and 353. The moment and padded-Be2 checks explain
the false start; they are not current packet or endpoint validation.

### HP-REP-XGTO-IMPORT-FN-01 — external GTO orbital import

Lifecycle: implemented. Permission: source maintenance.

Owner/canonical: `src/cartesian_external_gto_import.jl`;
[external GTO orbital import](external_gto_orbital_import.md).

Source: `src/cartesian_external_gto_import.jl`; existing include/export wiring
in `src/GaussletBases.jl`.

Dependencies: exact overlap/handoff kernels in `src/cartesian_gto_probes.jl`;
the cross-overlap-only convention in
`src/cartesian_representation_transfer.jl`; existing protected member data in
`src/cartesian_protected_ladder_bundle.jl` for the internal composition.

Validation: `test/nested/cartesian_external_gto_import_runtests.jl`.

Permission: maintain explicit packet validation, `S_FG*C_G` import,
spin-resolved capture, and direct unorthonormalized protected `S_LG*C_G`
composition.

Non-goals: packet-file/PySCF readers, solver orthonormalization, generalized
final metrics, Hamiltonian/interaction transforms, physics-policy changes, or
Cr2 endpoint claims.

### HP-REP-XGTO-IMPORT-TEST-01 — external GTO orbital import validation

Lifecycle: implemented validation. Permission: test maintenance.

Owner/canonical: external representation transfer;
[external GTO orbital import](external_gto_orbital_import.md).

Test: `test/nested/cartesian_external_gto_import_runtests.jl`.

Dependencies: `HP-REP-XGTO-IMPORT-FN-01` and repo-owned synthetic packets.

Permission: maintain packet identity/order/`S_GG`, restricted/spin-resolved
import, capture, rotation-invariance, and malformed-input checks.

Non-goals: PySCF, molecular endpoint, solver, energy, or Cr2 fixtures.

### HP-REP-XGTO-PROTECT-SIDECAR-FN-01 — protected external-GTO representation sidecar

Lifecycle: implemented internal facility. Permission: source maintenance.

Owner/canonical: `src/cartesian_external_gto_import.jl`;
[external GTO orbital import](external_gto_orbital_import.md), protected
representation sidecar section.

Source: `src/cartesian_external_gto_import.jl`.

Dependencies: exact final/GTO handoff in `src/cartesian_gto_probes.jl`, the
existing protected member contract in `src/cartesian_protected_ladder_bundle.jl`,
and existing protected-artifact/Coulomb serializers.

Validation: `test/nested/cartesian_external_gto_import_runtests.jl`.

Permission: maintain the standalone native-order, final-by-external v1 sidecar,
exact `S_LG`, direct spin imports, packet/member fingerprints, and metric-aware
capture diagnostics.

Non-goals: protected/ladder schema changes, raw `G_L/A_L`, H1/Vee transforms,
lossy overlap storage, public workflow, solver state, PySCF, or Cr2 claims.

### HP-REP-XGTO-PROTECT-SIDECAR-TEST-01 — protected external-GTO sidecar validation

Lifecycle: implemented validation. Permission: test maintenance.

Owner/canonical: protected external representation persistence;
[external GTO orbital import](external_gto_orbital_import.md).

Test: `test/nested/cartesian_external_gto_import_runtests.jl`.

Dependencies: `HP-REP-XGTO-PROTECT-SIDECAR-FN-01`, the protected artifact
reader, and centralized Coulomb-summary validation.

Permission: maintain exact key/identity, roundtrip, saved-overlap reimport,
rectangular capture, packet/member/artifact, and tamper-rejection checks.

Non-goals: Cr2-sized fixtures, solver/HF behavior, protected/ladder schema
tests, or external dependencies.

### HP-PQS-SCREEN-HARTREE-CORR-FN-01 - internal screened-Hartree correction

Lifecycle: implemented. Permission: source maintenance.

Owner/canonical: `CartesianReferenceDensity`;
[screened Hartree correction assembly](screened_hartree_correction_assembly.md).

Source: `src/cartesian_reference_density/screened_hartree_correction.jl`, its
module wiring, and narrow packet validation/field evaluation.

Permission: consume represented converged references and same-basis `V_IDA`,
`J0_G`, and `E0_G`; return in-memory
`Delta_J0 = J0_G - Diagonal(V_IDA*q0)` and
`C = 0.5*q0'V_IDA*q0 - 0.5*E0_G`; preserve strict representation,
finiteness, symmetry, convergence, and derivative/algebra failures while
reporting ordinary fitted-potential energy inconsistency.

Non-goals: physical kinetic/nuclear H1 changes, public workflow, corrected
artifacts, solvers, exchange, EGOI, row-gauge substitutes, source/interaction
transforms, `C' V C`, or Cr2 claims.

### HP-PQS-SCREEN-HARTREE-CORR-TEST-01 - correction validation

Lifecycle: completed validation contract. Permission: validation maintenance.

Owner/canonical: `CartesianReferenceDensity`;
[screened Hartree correction assembly](screened_hartree_correction_assembly.md).

Test: `test/nested/cartesian_screened_hartree_correction_runtests.jl`.

Permission: maintain packet/reference consistency, same-basis, anchor,
derivative, symmetry/finiteness, fitted-potential reporting, and malformed
input coverage without adding physics endpoint assertions.

## Approved Final Base Construction And Historical Handoff

### HP-FN-05 — final Hamiltonian construction

Lifecycle: implemented. Permission: source maintenance.

Owner/canonical: `CartesianIDAHamiltonian` and the staged base producer;
[terminal basis and base assembly](terminal_basis_and_base_assembly.md).

Sources:

- `src/cartesian_ida_hamiltonian.jl`;
- `src/cartesian_base_hamiltonian.jl`.

Validation:

- `test/ida/cartesian_ida_hamiltonian_runtests.jl`;
- `test/driver_public/cartesian_base_hamiltonian_runtests.jl`.

Dependencies: `HP-FN-03` one-body matrices, `HP-FN-04` IDA matrix, validated
charges/positions/electron counts, and the existing Hamiltonian object.

Permission: construct the existing `CartesianIDAHamiltonian` directly and
assemble `H1 = K + sum_A Z_A U_A` through its current accounting helpers.

Non-goals: wrapper result objects, route reports, new artifact schemas, solver
workflow, or restoration of retired Slice D materialization.

### HP-WIRE-02 — historical direct materialization Hamiltonian handoff

Status: retired by `e2e164e9b`. This entry is historical and no longer
authorizes source work.

Permission: none.

Owner/history:
[route-driver materialization retirement](route_driver_materialization_retirement.md).

Historically approved and implemented Slice D wrapper boundary:

```julia
cartesian_materialization(
    report,
    terminal_basis_realization,
    materialization_inputs,
)::Union{Nothing,CartesianIDAHamiltonian{Float64}}
```

This old route-driver wrapper workflow was retired under
`HP-RETIRE-DRV-MAT-*`. Current canonical producer work uses the staged
driver-facing producer functions and `CartesianIDAHamiltonian` artifact path;
do not restore `cartesian_materialization` or add compatibility callers.

The call site passes `transforms.terminal_basis_realization` directly. The
terminal basis must not be embedded in `cartesian_report`, reconstructed from
summaries, or recovered by passing the full `transforms` stage.

Return contract:

- no request returns `nothing`;
- requested base PQS materialization returns `CartesianIDAHamiltonian{Float64}`;
- `save_ham_artifact = true` writes with `write_cartesian_ida_hamiltonian` and
  still returns the same Hamiltonian;
- no materialization wrapper, `result_kind`, `materialized`, status mirror, or
  `ida_hamiltonian` field is approved.

## Implemented R1 Public Base Producer

### HP-R1-FILE-01 — public base producer source file

Lifecycle: implemented. Permission: source maintenance.

Owner/canonical: top-level `GaussletBases` public API;
[R1 public base producer](r1_public_base_producer.md).

Source: `src/cartesian_base_hamiltonian.jl`; export/include wiring in
`src/GaussletBases.jl`.

Validation: `test/driver_public/cartesian_base_hamiltonian_runtests.jl`.

Dependencies: implemented terminal/base assembly and existing Cartesian IDA
Hamiltonian object/writer.

Permission: maintain the public base producer in its current owner file and
the existing `cartesian_base_hamiltonian` export.

Non-goals: new files, exports, result types, driver behavior, supplements,
corrections, solver work, or artifact expansion.

### HP-R1-FN-01 — public base Hamiltonian producer facade

Lifecycle: implemented exported facade. Permission: source maintenance.

Owner/canonical: top-level `GaussletBases` public API;
[R1 public base producer](r1_public_base_producer.md).

Source: `src/cartesian_base_hamiltonian.jl` and export wiring in
`src/GaussletBases.jl`.

Validation: `test/driver_public/cartesian_base_hamiltonian_runtests.jl`.

Dependencies: `HP-R1-CORE-FN-01`, public `ns` normalization, atom/diatomic
scope contracts, mapping/Coulomb policies, and terminal/base assembly.

Permission: maintain
`cartesian_base_hamiltonian(system; basis, hamfile=nothing)` with plain
`NamedTuple` inputs and direct `CartesianIDAHamiltonian{Float64}` return.

Non-goals: public route/method selectors, wrapper payloads, general geometry,
heteronuclear/ECP inputs, supplements, corrections, solver work, or driver
input ownership.

### HP-R1-CORE-FN-01 — unified core-spacing producer contract

Lifecycle: implemented. Permission: source maintenance.

Owner/canonical: base input normalization;
[R1 public base producer](r1_public_base_producer.md).

Source: `src/cartesian_base_hamiltonian.jl`.

Validation: core-spacing and deprecated-`d` checks in
`test/driver_public/cartesian_base_hamiltonian_runtests.jl`.

Dependencies: one-center mapping and explicit atom/diatomic basis inputs.

Permission: maintain `core_spacing` as the single public near-nucleus physical
scale and atom-only compatibility `d == core_spacing`.

Non-goals: public `parent_mapping_d`, hidden element defaults, automatic
tuning, mapping-strength policy, ECP/solver behavior, or artifact changes.

### HP-PQS-MAP-SFACTOR-FN-01 — expert mapping `s_factor` keyword

Lifecycle: implemented. Permission: source maintenance.

Owner/canonical: parent mapping construction and provenance;
[PQS/WL mapping `s_factor`](pqs_mapping_s_factor.md).

Source:

- `src/mappings.jl`;
- `src/pqs_source_box_route_driver_helpers.jl`;
- `src/cartesian_base_hamiltonian.jl`;
- `bin/cartesian_ham_builder.jl` for the implemented expert input;
- `src/cartesian_protected_ladder_bundle.jl` for recipe provenance only.

Permission: maintain finite positive `s_factor`, default `1.0`, one-center
`effective_s = s_factor*sqrt(Z*core_spacing)`, the analogous per-center
multicenter combined-inverse-sqrt input, and explicit standard/effective
provenance.

Non-goals: element defaults, automatic tuning, revived public mapping internals,
solver/EGOI/rho0, protected-interaction, residual/injection, or Cr2-specific
policy.

### HP-PQS-MAP-SFACTOR-TEST-01 — mapping `s_factor` validation

Lifecycle: completed validation contract. Permission: validation maintenance.

Owner/canonical: mapping-factor validation;
[PQS/WL mapping `s_factor`](pqs_mapping_s_factor.md).

Evidence: default parity, nonunit one-center mapping/provenance, multicenter
fit/provenance, invalid-input, artifact/readback, and protected-recipe smokes.

### HP-PQS-COULOMB-ACCURACY-FN-01 - producer-wide Coulomb accuracy policy

Status: compact/high producer policy implemented; fixed standard tier and
narrow canonical-driver exposure approved for implementation.

Permission: source maintenance, including only completion of the already
approved Standard60 fingerprint/provenance and canonical-driver exposure.

Owner/canonical: Cartesian producer-wide Coulomb policy;
[Coulomb accuracy policy](coulomb_accuracy_policy.md).

Purpose: let expert consumers select a fixed compact, standard, or high
Coulomb Gaussian expansion while preserving one internally consistent
approximation from parent/PGDG construction through base IDA and
residual-GTO/MWG augmentation.

Approved producer input:

```text
coulomb_accuracy = :compact | :standard | :high
default = :compact
```

Exact presets:

```text
:compact -> doacc=false, 45 terms, del=0.6, s=0.5, c=0.03, maxu=27.0
:standard -> doacc=false, 60 terms, del=1.0,
             s=0.34257593251905827, c=0.042605721927199074, maxu=60.0
:high    -> doacc=true, 135 terms, del=1.0, s=0.16, c=0.01, maxu=135.0
```

The standard preset is the fixed analytic K60 midpoint quadrature. Its
canonical little-endian Float64 coefficient-then-exponent SHA-256 fingerprint is
`2de3ec44fc3d6b11ea26b7551e6b5ddef8bb2de1898fe0702d65f91cbf6c0f3a`.
The canonical contract owns its exact operation order and evidence. `doacc`
is legacy compatibility metadata, not preset identity.

Only those names are user-facing. This ID does not approve user inputs for
`doacc`, `del`, `s`, `c`, `maxu`, coefficients, exponents, or custom
expansion objects.

Canonical-driver amendment: `bin/cartesian_ham_builder.jl` may expose the same
`coulomb_accuracy` symbol with visible default `:compact`, accept it through
the existing trusted input-file/`key=value` allowlist, validate only
`:compact | :standard | :high`, add it to `common_basis`, and print it in the
existing basis contract summary. The driver passes the symbol to the producer
and must not resolve or inspect the expansion itself.

The option is route-family-neutral wherever current PQS and White-Lindsey
construction paths share parent/base/supplemented machinery. Neither route may
re-resolve a different expansion.

Approved construction behavior:

- resolve one existing `CoulombGaussianExpansion` before
  `cartesian_parent(...)` and parent-axis PGDG factor construction;
- carry that object in the current base working-basis construction instead of
  copying its scalar summary into route/stage records;
- build every parent-axis PGDG packet from its exponents;
- use the same expansion for base unit-nuclear attraction, base IDA
  electron-electron assembly, residual-GTO mixed/self Coulomb blocks,
  augmented unit-nuclear blocks, and residual MWG interaction;
- validate exact term-count and exponent-order parity at PGDG/base/augmentation
  boundaries and fail on mismatch;
- replace the MWG blanket rejection of explicit expansions with a parity check
  against the parent PGDG expansion;
- preserve deletion of the caller-free `_cartesian_base_ida_hamiltonian(...)`
  helper; no private replacement may independently select compact accuracy.

Stable-formula amendment:

- `GaussianAnalyticIntegrals.gaussian_factor` may replace its
  cancellation-prone weighted variance with the algebraically equivalent
  pairwise weighted squared-distance identity;
- `GaussianAnalyticIntegrals.gaussian_pair_factor` must form
  `D = alpha_a*alpha_b + 2g*(alpha_a + alpha_b)` directly and use the
  corresponding center-difference exponent rather than subtracting two
  `O(g^2)` terms and two nearly equal exponent terms;
- `CartesianGaussianRawBlocks._factor_axis_integral` may use the same
  pairwise weighted-distance identity for its three Gaussian factors;
- moderate-exponent numerical meaning, polynomial moments, normalization, and
  prefactors remain unchanged;
- do not clamp, take absolute values, truncate tight terms, or add a scaled/log
  carrier to hide cancellation.

Existing ordinary/Qiu-White callers may inherit the stable shared-kernel
evaluation, but this does not approve route-specific rewiring, cleanup,
defaults, or new Qiu-White validation infrastructure.

The bounded audit found H/H2 high-expansion factor scales below `2.6`,
finite/symmetric base and nuclear matrices, and finite supplemented
residual-GTO/MWG assembly in about `36.1 s`. This evidence removes the
reported need for a new PGDG carrier or terminal-contraction design.

The expansion must not become a field of `CartesianIDAHamiltonian`. That
object remains the finished numerical Hamiltonian. One compact expansion object
may cross construction stages; one compact summary may cross artifact
boundaries.

Approved artifact behavior:

- new base and supplemented artifacts write one Hamiltonian-wide
  `coulomb_expansion/` summary with `policy`, `doacc`, `term_count`,
  `del`, `s`, `c`, `maxu`, and the coefficient/exponent `fingerprint`;
- supplemented artifacts do not write separate base and augmentation
  expansion policies;
- ordinary matrix-only Cartesian readback may remain unchanged;
- protected-localized artifacts and protected ladder members/manifests preserve
  and expose the summary on readback;
- legacy artifacts without the summary remain readable where already readable;
  compact/high summaries without fingerprints remain legacy-readable with the
  fingerprint unavailable, but no fingerprint or missing policy is inferred
  as `:standard` or `:high`;
- a summary claiming `:standard` without the exact K60 fingerprint fails;
- new summaries must validate against the named deterministic preset.

Atomic HF reference packets are a deliberate exception to the one-Hamiltonian
summary because they record separately evaluated reference objects. Packet
writer/readback must record packet pure-GTO RHF as `:high`,
density-fit/self-energy evaluation as `:compact`, and the
fitted-potential broad-tail scaffold as `:compact`. Those packet-local
compact evaluations remain allowed under the measured Cr screened
scalar-constant error of about `0.0402 mHa`; they do not authorize mixed
compact/high construction inside one produced Hamiltonian.

Approved source surface:

```text
src/cartesian_base_hamiltonian.jl
src/pqs_source_box_route_driver_helpers.jl
src/pqs_source_box_low_order_materialization.jl
src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl
src/cartesian_residual_gaussians/mwg_interaction.jl
src/cartesian_ida_hamiltonian.jl
src/cartesian_protected_ladder_bundle.jl
src/cartesian_reference_density/atomic_hf_reference_packets.jl
src/ordinary_coulomb.jl
src/GaussianAnalyticIntegrals.jl
src/cartesian_gaussian_raw_blocks/nuclear_blocks.jl
bin/cartesian_ham_builder.jl
```

`src/cartesian_ida_hamiltonian.jl` is approved only for compact expansion
summary serialization/readback shared by current artifact owners. No new source
file, struct, public export, canonical driver input other than the exact policy
symbol above, general report object, or route-stage field cloud is approved.

Forbidden:

- changing the producer default to `:high`;
- custom expansion parameters or coefficient/exponent input;
- other canonical driver inputs or CLI changes;
- ordinary Qiu-White, legacy, or experimental path cleanup;
- scaled/log PGDG carriers, new stage objects, or terminal contraction changes;
- shellification, terminal realization, retained selection, mapping,
  residual-selection, injection, EGOI, or screened-Hartree formula changes;
- solver/HF/MP2-NO workflow;
- protected-localized interaction or ladder-transfer semantic changes;
- Cr/Cr2-specific producer branches or committed endpoint claims.

Failure rule: stop without a source commit if a single carried expansion cannot
serve parent/PGDG, base, residual-GTO, and MWG construction without a new broad
carrier, if MWG cannot prove PGDG exponent parity, if legacy provenance would
need to be guessed, or if files outside the approved surface are required.

Target source growth is at most about 300 added lines, offset where practical
by deleting independent compact selectors and the uncalled private base helper.
The stable-formula amendment itself should remain below about 60 added source
lines and must not introduce a new carrier, cache, status object, or module.

### HP-PQS-COULOMB-ACCURACY-TEST-01 - Coulomb accuracy validation

Status: approved validation authority.

Permission: validation maintenance for implemented compact/high behavior and
the already-approved Standard60/driver completion only.

Owner/canonical: Coulomb policy validation;
[Coulomb accuracy policy](coulomb_accuracy_policy.md).

Approved committed test surfaces:

```text
test/driver_public/cartesian_base_hamiltonian_runtests.jl
test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl
test/nested/cartesian_atomic_hf_reference_packet_runtests.jl
test/core/runtests.jl
test/docs/cartesian_ham_builder_policy_runtests.jl
```

Approved validation:

- `git diff --check`;
- package load;
- omitted policy and explicit `:compact` produce equal matrices in a bounded
  base case;
- bounded `:standard` construction records the exact K60 parameters, 60-term
  count, fingerprint, and parent/PGDG exponent parity;
- bounded `:high` base construction records the exact 135-term summary and
  parent/PGDG exponent parity;
- bounded White-Lindsey base construction confirms the policy reaches the
  shared parent/base path;
- bounded supplemented construction proves one expansion reaches base,
  residual-GTO exact blocks, augmented unit-nuclear work, and MWG;
- protected-localized member and ladder manifest write/readback preserve the
  summary;
- atomic packet roundtrip preserves separate RHF, density/self-energy, and
  potential-tail expansion provenance;
- at least one bounded compact/standard/high comparison reports term-dependent
  timing and allocation;
- a small core analytic-kernel test compares compact/standard/high exponent
  cases to a BigFloat oracle, preserves moderate-exponent parity, reproduces
  the exact standard fingerprint, and covers finite nonnegative s-type factors
  plus translated-center cancellation cases;
- every endpoint-style probe used for interpretation includes terminal
  due-diligence review.

The existing docs policy test may check the driver input/default/validation/
forwarding/printing contract without executing a Hamiltonian. Bounded
canonical-driver validation must separately show omitted-versus-explicit
compact matrix/artifact parity, one accepted standard request with exact
K60 provenance, and one accepted high request with finite/symmetric matrices
and `:high`/`135`-term artifact provenance. Use temporary or ignored outputs;
do not add a committed driver fixture or endpoint value.

High supplemented and protected-ladder checks may remain ignored bounded probes
if committed execution would materially increase routine test time. No new
committed test file, solver run, Cr/Cr2 endpoint assertion, or production
energy claim is approved.

### HP-R1-WIRE-01 — report-free base producer wiring

Lifecycle: implemented. Permission: source maintenance.

Owner/canonical: staged base composition in
`src/cartesian_base_hamiltonian.jl`;
[R1 public base producer](r1_public_base_producer.md).

Source: `src/cartesian_base_hamiltonian.jl`.

Validation: `test/driver_public/cartesian_base_hamiltonian_runtests.jl` and
downstream terminal/base endpoint gates.

Dependencies: `cartesian_base_working_basis`, staged product/unit-nuclear/Vee
construction, `cartesian_base_hamiltonian_assembly`, and one carried Coulomb
expansion.

Permission: maintain the report-free staged path used by the public facade and
human-facing driver.

Non-goals: retired `_cartesian_base_ida_hamiltonian`, pair/assembly reports,
duplicate Hamiltonian construction, exposed stages, status payloads, or new
artifact behavior.

### HP-ROUTE-RECIPE-FN-01 — family-selective route recipes

Lifecycle: implemented by `c6307a16d`. Permission: source maintenance.

Owner/canonical: route recipe construction;
[nesting/supplement composition](nesting_supplement_composition_plan.md).

Source: `src/pqs_source_box_route_driver_helpers.jl` and narrow base route
normalization in `src/cartesian_base_hamiltonian.jl`.

Permission: build only the selected `:pqs_source_box` or
`:white_lindsey_low_order` subrecipe and leave inactive family vocabulary
absent or `nothing` without merging the algorithms.

Non-goals: driver inputs, numerical/shell policies, artifacts, reports,
materialization revival, or route-stage redesign.

### HP-ROUTE-RECIPE-TEST-01 — route recipe validation

Lifecycle: completed validation contract. Permission: validation maintenance.

Owner/canonical: route recipe construction;
[nesting/supplement composition](nesting_supplement_composition_plan.md).

Evidence: bounded family-selective recipe, atom/H2 base, and supplemented
artifact/readback smokes accepted with implementation `c6307a16d`.

### HP-ROUTE-INV-FN-01 — retained-unit route inventory

Lifecycle: implemented by `c985723c7`. Permission: source maintenance.

Owner/canonical: `src/pqs_source_box_route_driver_helpers.jl`;
[route/stage metadata](route_stage_metadata_contract.md).

Permission: maintain vector-backed ordered retained-unit and pair-family rows
with label lookup; labels remain data rather than concrete type parameters.

Non-goals: recipe or shell policy, numerical behavior, public/driver inputs,
artifacts, reports, compatibility shapes, or Cr2 workflow.

### HP-ROUTE-INV-TEST-01 — route inventory validation

Lifecycle: completed validation contract. Permission: validation maintenance.

Owner/canonical: route inventory;
[route/stage metadata](route_stage_metadata_contract.md).

Evidence: focused type-shape/order checks and bounded base/supplemented
artifact gates accepted with `c985723c7`.

### HP-RAW-SRCMODE-FN-01 — raw product source-mode inventory

Lifecycle: implemented by `34cf8f106`. Permission: source maintenance.

Owner/canonical: raw-product source planning;
[route/stage metadata](route_stage_metadata_contract.md).

Source:

- `src/cartesian_raw_product_sources/records.jl`;
- `src/cartesian_raw_product_sources/source_mode_indices.jl`;
- `src/cartesian_raw_product_sources/summaries.jl`;
- narrow storage-change consumers only in
  `src/cartesian_retained_unit_transform_contracts/unit_contracts.jl`,
  `src/cartesian_final_basis_realization/pqs_terminal_basis_realization.jl`,
  and `src/cartesian_base_hamiltonian.jl`.

Permission: maintain vector-backed mode/column inventories while preserving
fixed `NTuple{3,Int}` coordinates, deterministic mode order, retained-rule
association, and manifest source provenance.

Non-goals: pair/source-box redesign, numerical behavior, public/driver inputs,
artifact schema, compatibility tuple shapes, or Cr2 workflow.

### HP-RAW-SRCMODE-TEST-01 — raw source-mode validation

Lifecycle: completed validation contract. Permission: validation maintenance.

Owner/canonical: raw product source modes;
[route/stage metadata](route_stage_metadata_contract.md).

Evidence: focused mode/order/retained-rule checks plus bounded base,
supplemented, R3, and manifest validation accepted with `34cf8f106`.

### HP-CONTRACT-VEC-FN-01 — vector-backed contract plans

Lifecycle: implemented by `5938ddbc5`. Permission: source maintenance.

Owner/canonical: terminal lowering and retained-unit transform contracts;
[route/stage metadata](route_stage_metadata_contract.md).

Source:

- `src/cartesian_terminal_lowering/contracts.jl`;
- `src/cartesian_terminal_lowering/selection.jl`;
- `src/cartesian_terminal_lowering/summaries.jl`;
- `src/cartesian_retained_unit_transform_contracts/records.jl`;
- `src/cartesian_retained_unit_transform_contracts/unit_contracts.jl`;
- `src/cartesian_retained_unit_transform_contracts/summaries.jl`;
- narrow storage-change consumers only in
  `src/pqs_source_box_route_driver_helpers.jl`,
  `src/cartesian_base_hamiltonian.jl`, and
  `src/cartesian_final_basis_realization/pqs_terminal_basis_realization.jl`.

Permission: maintain vector-backed available/selected lowering contracts and
retained-unit transform contracts with unchanged accessor and order semantics.
Per-contract `source_cpbs` and fixed mathematical tuples remain outside this
cleanup.

Non-goals: route or shell policy, numerical behavior, public/driver inputs,
artifact schema, compatibility tuple shapes, or Cr2 workflow.

### HP-CONTRACT-VEC-TEST-01 — contract-plan validation

Lifecycle: completed validation contract. Permission: validation maintenance.

Owner/canonical: terminal lowering and transform contracts;
[route/stage metadata](route_stage_metadata_contract.md).

Evidence: lowering/transform order checks and bounded base, supplemented, and
R3 gates accepted with `5938ddbc5`.

### HP-ROUTE-STAGE-TYPE-FN-01 — route/stage type surface

Lifecycle: implemented by `118a639bf`. Permission: source maintenance.

Owner/canonical: route helpers and terminal shellification inventory;
[route/stage metadata](route_stage_metadata_contract.md).

Source: `src/pqs_source_box_route_driver_helpers.jl` and
`src/cartesian_terminal_shellification_geometry.jl`.

Permission: preserve compact vector-backed route/shellification summaries and
narrow stage returns without duplicate lowering-plan ownership.

Non-goals: route or shell policy, numerical behavior, public/driver inputs,
artifacts, reports, precompile machinery, or Cr2 workflow.

### HP-ROUTE-STAGE-TYPE-TEST-01 — type-surface validation

Lifecycle: completed validation contract. Permission: validation maintenance.

Owner/canonical: route/stage type surfaces;
[route/stage metadata](route_stage_metadata_contract.md).

Evidence: focused shape/order scans and bounded base/supplemented gates
accepted with `118a639bf`.

### HP-ROUTE-STAGE-CARRIER-FN-01 — route/stage carriers

Lifecycle: implemented by `8c3df2ad9`. Permission: source maintenance.

Owner/canonical: route helpers, active complete-core-shell support planning,
and terminal realization;
[route/stage metadata](route_stage_metadata_contract.md).

Source: `src/pqs_source_box_route_driver_helpers.jl`,
`src/pqs_source_box_diatomic_complete_core_shell.jl`, and, only when directly
required to slim the approved carrier path,
`src/cartesian_final_basis_realization/pqs_terminal_basis_realization.jl`.

Permission: keep only live compact plans/realizations/summaries across stage
boundaries. This does not retire route skeletons or pair/assembly/report
stages.

Non-goals: route/stage or tool retirement, numerical or shell policy,
public/driver changes, artifacts, reports, precompile machinery, or Cr2.

### HP-ROUTE-STAGE-CARRIER-TEST-01 — carrier validation

Lifecycle: completed validation contract. Permission: validation maintenance.

Owner/canonical: route/stage carriers;
[route/stage metadata](route_stage_metadata_contract.md).

Evidence: terminal support/shellification/lowering order and bounded
base/supplemented/R3 validation accepted with `8c3df2ad9`.

### HP-R1-ART-01 — public base producer artifact provenance

Lifecycle: implemented. Permission: source maintenance.

Owner/canonical: base artifact composition;
[R1 public base producer](r1_public_base_producer.md).

Source: `src/cartesian_base_hamiltonian.jl`; minimal writer/readback in
`src/cartesian_ida_hamiltonian.jl`.

Validation: artifact/provenance section of
`test/driver_public/cartesian_base_hamiltonian_runtests.jl`.

Dependencies: existing version-1 Cartesian IDA artifact and separately owned
manifest/Coulomb summaries.

Permission: maintain the fixed `producer_provenance/` keys and truthful
route/size/mapping/system values listed in the canonical R1 contract.

Non-goals: new artifact kinds or matrix keys, separate provenance files,
public provenance readers, wrappers, algorithm consumption of sidecars, or
manifest schema expansion.

### HP-R1-TEST-01 — public base producer endpoint test/example

Lifecycle: implemented validation contract. Permission: validation
maintenance.

Owner/canonical: R1 public base producer;
[R1 public base producer](r1_public_base_producer.md).

Source: no production source permission.

Validation: `test/driver_public/cartesian_base_hamiltonian_runtests.jl`.

Dependencies: exported facade, existing Hamiltonian writer/readback,
`producer_provenance/`, and current compact/high Coulomb source behavior.

Permission: maintain the standalone atom/H2 endpoint, malformed-input,
deprecated-`d`, geometry, matrix, Coulomb, artifact, and provenance checks.

Non-goals: default test-suite wiring, private stage/report assertions, driver
input tests, supplemented/solver/Cr2 gates, or new source behavior.

## Implemented R1 One-Center Base Atoms

Canonical contract:
[R1 one-center base atoms](r1_one_center_base_atoms.md).

### HP-R1-ATOM-FN-01 — explicit one-center all-electron base atom facade

Lifecycle: implemented. Permission: source maintenance.

Canonical contract:
[R1 one-center base atoms](r1_one_center_base_atoms.md).

Owner/source: `src/cartesian_base_hamiltonian.jl`.

Permission: maintain explicit origin-centered, neutral, all-electron atom
validation in the existing `cartesian_base_hamiltonian(system; basis,
hamfile)` facade. Charge, electron counts, spin sectors, basis, and ECP behavior
must never be inferred from the atom label.

Non-goals: translated atoms, element defaults, supplements/corrections, ECP,
solver workflow, API redesign, artifact expansion, or atom-only construction.

### HP-R1-ATOM-WIRE-01 — one-center atom shared workflow wiring

Lifecycle: implemented. Permission: source maintenance.

Canonical contract:
[R1 one-center base atoms](r1_one_center_base_atoms.md).

Owner/source: `src/cartesian_base_hamiltonian.jl`.

Permission: map explicit charge and `core_spacing` into the private atomic
mapping, derive physical parent extent from `radius`, and use the same terminal,
one-body, IDA, Hamiltonian, writer, and provenance machinery as the supported
base producer. Atom routes remain `:one_center_pqs_base` or
`:one_center_wl_base` according to nesting.

### HP-R1-ATOM-TEST-01 — one-center base atom validation

Lifecycle: completed validation contract. Permission: validation maintenance.

Owner/canonical: one-center base atom validation;
[R1 one-center base atoms](r1_one_center_base_atoms.md).

Validation: `test/driver_public/cartesian_base_hamiltonian_runtests.jl` plus
accepted bounded non-H atom smokes recorded in the manager log.

Permission: maintain H regression, malformed atom/basis input rejection,
finite/symmetric operator, mapping/provenance, and artifact/readback checks.

Non-goals: committed non-H reference energies, translated/ECP/supplemented
fixtures, driver tests, solvers, or Cr2 gates.

## Implemented R3 Compatibility And Endpoint History

### HP-R3-OBJ-01 — residual-GTO augmentation object

Lifecycle: implemented compatibility alias. Permission: source maintenance.

Owner/canonical: `CartesianResidualGaussians`, with terminal compatibility;
[R3 compatibility history](r3_residual_gto_mwg_augmentation.md) and
[Residual Gaussian domain module](residual_gaussian_domain_module.md).

Source: `src/cartesian_residual_gaussians/residual_basis.jl` and
`src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl`.

Validation: `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl`.

Dependencies: `HP-RG-OBJ-01` and the owner-local residual basis contract.

Permission: maintain `CartesianTerminalResidualGTOAugmentation` only as the
live alias of `CartesianResidualGaussianBasis` required by compatibility
callers.

Non-goals: a second object schema, R3-named numerical fields, payload/status
objects, artifacts, or public API.

### HP-R3-FN-01 — residual-basis construction

Lifecycle: implemented historical/compatibility surface. Permission: source
maintenance.

Owner/canonical: `CartesianResidualGaussians`;
[Residual Gaussian domain module](residual_gaussian_domain_module.md).

Source: `src/cartesian_residual_gaussians/residual_basis.jl` and the composition
seam in `src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl`.

Validation: `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl`.

Dependencies: `HP-RG-FN-01`, exact mixed overlap, and owner identity.

Permission: maintain delegation to owner-local residual selection and one final
inter-owner merge.

Non-goals: global candidate selection, numerical-complete reinterpretation,
new cutoffs, or duplicate R3 basis logic.

### HP-R3-FN-02 — exact augmented one-body and moment assembly

Lifecycle: implemented historical/compatibility surface. Permission: source
maintenance.

Owner/canonical: `CartesianResidualGaussians`;
[Residual Gaussian domain module](residual_gaussian_domain_module.md).

Source: `src/cartesian_residual_gaussians/augmented_operators.jl` and raw-block
composition in
`src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl`.

Validation: `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl`.

Dependencies: `HP-RG-FN-02`, exact raw `[G,A]` blocks, and the residual
transform.

Permission: maintain exact augmented kinetic, by-center unit-nuclear,
coordinate, and second-moment transformations.

Non-goals: MWG approximation, new raw kernels, artifact fields, or interaction
rotation.

### HP-R3-FN-03 — residual MWG/IDA and in-memory Hamiltonian

Lifecycle: implemented compatibility entry point. Permission: source
maintenance.

Owner/canonical: terminal composition over `CartesianResidualGaussians`;
[R3 compatibility history](r3_residual_gto_mwg_augmentation.md).

Source: `src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl`
and `src/cartesian_residual_gaussians/mwg_interaction.jl`.

Validation: `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl`.

Dependencies: `HP-RG-FN-03`, `HP-RG-FN-04`, same-construction base inputs,
and the producer-owned Coulomb expansion.

Permission: maintain
`pqs_terminal_residual_gto_augmented_hamiltonian(...)` and direct return of the
existing `CartesianIDAHamiltonian{Float64}`.

Non-goals: duplicate MWG math, wrapper results, post-hoc opaque-Hamiltonian
augmentation, solver work, or public API.

### HP-R3-ART-01 — compact supplemented artifact provenance

Lifecycle: implemented internal artifact compatibility surface. Permission:
source maintenance.

Owner/canonical: `CartesianFinalBasisRealization` workflow;
[Cartesian Hamiltonian artifact manifest](cartesian_hamiltonian_artifact_manifest.md).

Source:
`src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl`.

Validation: the artifact/readback section of
`test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl`.

Dependencies: existing Cartesian IDA artifact, `HP-HAM-MANIFEST-*`, and a
validated residual construction.

Permission: maintain the compact `supplement_provenance/` group and existing
writer composition recorded in the artifact contract.

Non-goals: residual transforms, dense moments, matched-Gaussian arrays, full
spectra, new artifact kinds, public readers, or RG ownership of persistence.

### HP-R3-TEST-01 — residual-GTO/MWG compatibility endpoint

Lifecycle: implemented validation contract. Permission: validation
maintenance.

Owner/canonical: R3 compatibility family;
[R3 compatibility history](r3_residual_gto_mwg_augmentation.md).

Source: no production source permission.

Validation: `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl`.

Dependencies: the current RG domain functions, terminal compatibility seams,
and existing Hamiltonian/artifact readers.

Permission: maintain the standalone H2 residual geometry, exact-operator,
independent interaction, in-memory Hamiltonian, and artifact compatibility
gate.

Non-goals: normal `Pkg.test` wiring, Be2/Cr2 committed fixtures, private
status vocabulary, or new source behavior.

## Cartesian Hamiltonian Artifact Manifest

### HP-HAM-MANIFEST-FN-01 — compact Hamiltonian artifact manifest

Lifecycle: implemented. Permission: source maintenance.

Owner/canonical: ordinary Cartesian Hamiltonian artifact composition;
[artifact manifest](cartesian_hamiltonian_artifact_manifest.md).

Source:

```text
src/cartesian_base_hamiltonian.jl
src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl
src/cartesian_ida_hamiltonian.jl
```

Validation:

```text
test/driver_public/cartesian_base_hamiltonian_runtests.jl
test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl
docs/src/developer/pqs_manager_running_log.md (Pass 113)
```

Dependencies: minimal `CartesianIDAHamiltonian` writer/reader, R1 producer,
and supplemented residual-GTO/MWG writer.

Permission: maintain matrix-order `final_basis_labels/` and
`recipe_provenance/` on existing facade artifacts.

Non-goals: matrix/reader changes, public manifest API, inferred labels, dense
transforms, solver/consumer fields, protected artifacts, or new composition.

### HP-HAM-MANIFEST-TEST-01 — artifact manifest validation

Lifecycle: validation completed; tracked regression coverage is partial.
Permission: validation maintenance.

Owner/canonical: ordinary artifact validation;
[artifact manifest](cartesian_hamiltonian_artifact_manifest.md).

Test/evidence:

```text
test/driver_public/cartesian_base_hamiltonian_runtests.jl
test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl
docs/src/developer/pqs_manager_running_log.md (Pass 113)
```

Dependencies: `HP-HAM-MANIFEST-FN-01` and the unchanged ordinary reader.

Permission: maintain base/supplemented readback and selected provenance-root
checks. The accepted direct-JLD2 evidence owns detailed field/status coverage.

Non-goals: a new committed test file, public reader, schema dump, Cr2 fixture,
or solver/consumer workflow.

### HP-HAM-MANIFEST-SRC-FN-01 — source-mode provenance seam

Lifecycle: partially implemented. Permission: source maintenance and the
remaining construction-native subset only.

Owner/canonical: terminal source provenance at the ordinary artifact boundary;
[artifact manifest](cartesian_hamiltonian_artifact_manifest.md).

Implemented source:

```text
src/cartesian_base_hamiltonian.jl
```

Remaining approved source surfaces:

```text
src/cartesian_terminal_lowering/contracts.jl
src/cartesian_terminal_lowering/region_contracts.jl
src/cartesian_raw_product_sources/CartesianRawProductSources.jl
src/cartesian_raw_product_sources/records.jl
src/cartesian_raw_product_sources/source_mode_indices.jl
src/cartesian_retained_units/CartesianRetainedUnits.jl
src/cartesian_retained_units/records.jl
src/cartesian_retained_units/lower_contract_units.jl
src/cartesian_retained_unit_transform_contracts/CartesianRetainedUnitTransformContracts.jl
src/cartesian_retained_unit_transform_contracts/records.jl
src/cartesian_retained_unit_transform_contracts/unit_contracts.jl
src/cartesian_final_basis_realization/CartesianFinalBasisRealization.jl
src/cartesian_final_basis_realization/pqs_terminal_basis_realization.jl
```

Validation/evidence:

```text
test/driver_public/cartesian_base_hamiltonian_runtests.jl
test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl
docs/src/developer/pqs_manager_running_log.md (Passes 115-116)
```

Dependencies: terminal retained-rule/raw-product facts, terminal basis matrix
order, and `HP-HAM-MANIFEST-FN-01`.

Permission: maintain the existing `source_mode_provenance` carrier and add
only missing construction-native relations or labels on the listed surfaces.

Non-goals: inferred ray/radial/locality policy, coefficients, weights, spans,
dense transforms, route reports, new stage objects, reader/driver changes, or
consumer-specific schema.

### HP-HAM-MANIFEST-SRC-TEST-01 — source-mode provenance seam validation

Lifecycle: validation completed for the implemented subset; tracked regression
coverage is partial. Permission: validation maintenance.

Owner/canonical: terminal source-provenance validation;
[artifact manifest](cartesian_hamiltonian_artifact_manifest.md).

Test/evidence:

```text
test/driver_public/cartesian_base_hamiltonian_runtests.jl
test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl
docs/src/developer/pqs_manager_running_log.md (Passes 115-116)
```

Dependencies: `HP-HAM-MANIFEST-SRC-FN-01` and the existing reader.

Permission: maintain native source-group presence, status, no-inference, and
readback checks for the implemented subset.

Non-goals: a new committed test file, public reader, broad route/report test,
Cr2 fixture, or solver/consumer workflow.

### HP-NEST-ART-FN-01 — nesting artifact-truth cleanup

Lifecycle: implemented. Permission: source maintenance.

Owner/canonical: ordinary artifact nesting and route provenance;
[artifact manifest](cartesian_hamiltonian_artifact_manifest.md).

Source:

```text
src/cartesian_base_hamiltonian.jl
src/cartesian_final_basis_realization/CartesianFinalBasisRealization.jl
```

Validation/evidence:

```text
test/driver_public/cartesian_base_hamiltonian_runtests.jl
test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl
docs/src/developer/pqs_manager_running_log.md (Pass 134)
```

Dependencies: public nesting normalization and the separately authorized
PQS/WL composition cells.

Permission: maintain truthful `nesting` and route labels in producer/recipe
provenance and the nesting-neutral final-basis module description.

Non-goals: new composition support, route/shell changes, matrix/reader changes,
driver/API changes, diagnostics, or Cr2 workflow.

### HP-NEST-ART-TEST-01 — nesting artifact-truth validation

Lifecycle: validation completed; tracked regression coverage is partial.
Permission: validation maintenance.

Owner/canonical: artifact nesting/route validation;
[artifact manifest](cartesian_hamiltonian_artifact_manifest.md).

Test/evidence:

```text
test/driver_public/cartesian_base_hamiltonian_runtests.jl
test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl
docs/src/developer/pqs_manager_running_log.md (Pass 134)
```

Dependencies: `HP-NEST-ART-FN-01`; supplemented WL validation remains owned by
`HP-COMP-SUPPWL-TEST-01`.

Permission: maintain small artifact/readback provenance checks for truthful
nesting and route labels.

Non-goals: new fixtures, public reader/schema dump, composition validation,
solver workflow, or Cr2 tests.

## Implemented R3 Usability Supplemented Workflow

### HP-R3U-FILE-01 — supplemented workflow source and validation files

Lifecycle: implemented. Permission: source and validation maintenance.

Owner/canonical: base-Hamiltonian composition;
[R3 usability supplemented workflow](r3_usability_supplemented_workflow.md).

Source: `src/cartesian_base_hamiltonian.jl` and the compatibility/artifact seam
in `src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl`.

Validation: `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl`.

Dependencies: implemented base producer, RG domain module, and existing
Hamiltonian/artifact owners.

Permission: maintain the non-exported facade in the existing source and test
files.

Non-goals: new files, public export, driver policy, result payloads, artifact
schema expansion, or solver work.

### HP-R3U-FN-01 — non-exported supplemented Hamiltonian facade

Lifecycle: implemented. Permission: source maintenance.

Owner/canonical: `src/cartesian_base_hamiltonian.jl`;
[R3 usability supplemented workflow](r3_usability_supplemented_workflow.md).

Source: `src/cartesian_base_hamiltonian.jl`.

Validation: the facade section of
`test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl`.

Dependencies: `HP-R3U-ZDI-FN-01` for diatomic scope, shared base input
normalization, and the R3/RG same-construction path.

Permission: maintain
`cartesian_residual_gto_mwg_hamiltonian(system; basis, supplement, hamfile)`
as a non-exported direct-`CartesianIDAHamiltonian{Float64}` facade.

Non-goals: public API, wrapper/result objects, general molecular geometry,
heteronuclear/ECP inputs, solver work, or Cr2-specific behavior.

### HP-R3U-WIRE-01 — base-to-RG same-construction workflow

Lifecycle: implemented. Permission: source maintenance.

Owner/canonical: base-Hamiltonian composition;
[R3 usability supplemented workflow](r3_usability_supplemented_workflow.md).

Source: `src/cartesian_base_hamiltonian.jl` and
`src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl`.

Validation: `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl`.

Dependencies: one base construction, named supplement representation, RG
domain functions, producer-owned Coulomb expansion, and existing artifact
writer.

Permission: maintain one same-construction path from validated input through
the direct in-memory Hamiltonian and optional existing artifact.

Non-goals: post-hoc opaque-Hamiltonian augmentation, exposed stages, persistent
caches, payloads, or duplicate numerical algorithms.

### HP-R3U-TEST-01 — supplemented facade endpoint

Lifecycle: implemented validation contract. Permission: validation
maintenance.

Owner/canonical: R3 usability family;
[R3 usability supplemented workflow](r3_usability_supplemented_workflow.md).

Source: no production source permission.

Validation: the facade/artifact section of
`test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl`.

Dependencies: `HP-R3-TEST-01`, implemented facade, existing writer/readback,
and compact provenance contracts.

Permission: maintain H2 facade type, endpoint, malformed-input,
artifact/readback, and provenance checks.

Non-goals: normal test-suite wiring, new committed files, Be2/Cr2 committed
gates, private stage assertions, or source expansion.

## Implemented Residual Gaussian Domain

### HP-RG-FILE-01 — Residual Gaussian module files

Lifecycle: implemented. Permission: source maintenance.

Owner/canonical: `CartesianResidualGaussians`;
[RG domain](residual_gaussian_domain_module.md).

Source:

```text
src/cartesian_residual_gaussians/CartesianResidualGaussians.jl
src/cartesian_residual_gaussians/residual_basis.jl
src/cartesian_residual_gaussians/augmented_operators.jl
src/cartesian_residual_gaussians/mwg_interaction.jl
```

Validation: `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl`.
Dependencies: terminal basis, neutral raw blocks, and producer Coulomb policy.
Permission: maintain the internal module and owner files. Non-goals: public
export, facade parsing, artifacts, routes, protected policy, or solver work.

### HP-RG-OBJ-01 — residual Gaussian basis object

Lifecycle: implemented. Permission: source maintenance.

Owner/canonical: `CartesianResidualGaussianBasis`;
[RG domain](residual_gaussian_domain_module.md).
Source: `src/cartesian_residual_gaussians/residual_basis.jl`.
Validation: existing H2 endpoint. Dependencies: `HP-RG-FILE-01` and the
orthogonality/cutoff policy. Permission: maintain the exact live numerical
fields and invariants. Non-goals: status/report, route, artifact, or public
input fields.

### HP-RG-FN-01 — residual Gaussian basis construction

Lifecycle: implemented. Permission: source maintenance.

Owner/canonical: `build_residual_gaussian_basis(...)`;
[RG domain](residual_gaussian_domain_module.md).
Source: `src/cartesian_residual_gaussians/residual_basis.jl`.
Callers: terminal compatibility and protected-ladder composition.
Tests: H2 endpoint; `test/misc/runtests.jl` only for the separate explicit
numerical-complete policy. Dependencies: physical owner metadata and validated
`X/S_AA`. Permission: maintain owner-local selection and one final merge.
Non-goals: global selection, flooring, localization, or default changes.

### HP-RG-FN-02 — exact augmented operator transformation

Lifecycle: implemented. Permission: source maintenance.

Owner/canonical: `transform_augmented_operator(...)`;
[RG domain](residual_gaussian_domain_module.md).
Source: `src/cartesian_residual_gaussians/augmented_operators.jl`.
Callers/tests: terminal and protected-ladder composition; H2 endpoint.
Dependencies: exact `GG/GA/AA` blocks and a validated residual object.
Permission: maintain exact native-order kinetic, unit-nuclear, and moment
transforms. Non-goals: MWG, interaction rotation, artifact, or solver policy.

### HP-RG-FN-03 — moment-matched Gaussian descriptors

Lifecycle: implemented. Permission: source maintenance.

Owner/canonical: `moment_matched_gaussians(...)`;
[RG domain](residual_gaussian_domain_module.md).
Source: `src/cartesian_residual_gaussians/mwg_interaction.jl`.
Validation: H2 endpoint. Dependencies: final merged residual basis and exact
position/second moments. Permission: maintain finite positive final-residual
center/width descriptors. Non-goals: exact four-index integrals, basis
replacement, pre-merge descriptors, or localization policy.

### HP-RG-FN-04 — residual IDA interaction assembly

Lifecycle: implemented. Permission: source maintenance.

Owner/canonical: `assemble_residual_ida_interaction(...)`;
[RG domain](residual_gaussian_domain_module.md).
Source: `src/cartesian_residual_gaussians/mwg_interaction.jl`.
Validation: H2 endpoint and independent weight-aware `V_GM` reconstruction.
Dependencies: final descriptors, terminal supports/weights, and the one
producer Coulomb expansion. Permission: maintain unchanged `V_GG`,
weight-aware `V_GM`, and symmetric `V_MM`. Non-goals: `C' V C`, exact Vee,
artifacts, or solver policy.

### HP-RG-WIRE-01 — terminal and facade compatibility wiring

Lifecycle: implemented. Permission: compatibility/caller maintenance.

Owner/canonical: RG delegation boundary;
[RG domain](residual_gaussian_domain_module.md).
Source/callers:
`src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl`,
`src/cartesian_base_hamiltonian.jl`, and `bin/cartesian_ham_builder.jl`.
Validation: H2 facade/endpoint. Dependencies: `HP-RG-FN-01` through
`HP-RG-FN-04`. Permission: retain narrow live composition, validation, and
writer hooks. Non-goals: moving facade/artifact ownership into RG or restoring
migrated R3 physics helpers.

### HP-RG-TEST-01 — Residual Gaussian endpoint validation

Lifecycle: completed validation; active maintenance. Permission: test
maintenance, not source authority.

Owner/canonical: RG endpoint; [RG domain](residual_gaussian_domain_module.md).
Test: `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl`.
Dependencies: implemented RG functions and terminal composition.
Permission: maintain dimensions, exact operators, cutoff/provenance,
weight-aware `V_GM`, self-Coulomb, facade, and readback checks. Non-goals: new
fixture files, committed Be2/Cr2 gates, or private migration vocabulary.

### HP-RG-ORTHO-FN-01 — residual final-orthogonality robustness

Lifecycle: implemented current robustness. Permission: source maintenance.

Owner/canonical: final RG merge and identity validation;
[orthogonality/cutoff policy](residual_gaussian_orthogonality_robustness.md).
Source:
`src/cartesian_residual_gaussians/residual_basis.jl`, with compatibility
keyword plumbing in
`src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl`.
Validation: H2 endpoint and manager-log Passes 131-132.
Dependencies: positive owner metrics and one final merge.
Permission: maintain symmetric inverse-square-root normalization, hard
near-singular failure, and scale-aware final identity validation. Non-goals:
cutoff/default changes, flooring, selection changes, or broader workflow.

### HP-RG-ORTHO-TEST-01 — residual final-orthogonality validation

Lifecycle: completed evidence; active maintenance. Permission: test
maintenance, not source authority.

Owner/canonical: final RG numerical validation;
[orthogonality/cutoff policy](residual_gaussian_orthogonality_robustness.md).
Test: `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl`.
Evidence: manager-log Passes 131-132 retain the ignored N2 audit.
Dependencies: `HP-RG-ORTHO-FN-01` and the H2 endpoint.
Permission: maintain final `G-R`, `R-R`, merge, and endpoint checks.
Non-goals: new N2/Cr2 committed fixtures, source work, artifact schema,
driver, solver, ECP, or EGOI work.

### HP-RG-IDTOL-FN-01 — residual final-identity tolerance default

Lifecycle: implemented historical; superseded. Permission: none.

Owner/canonical: historical `identity_atol = 1e-8` transition;
[orthogonality/cutoff policy](residual_gaussian_orthogonality_robustness.md).
Historical source: `src/cartesian_residual_gaussians/residual_basis.jl` and
terminal keyword plumbing. Implementation: `47e56593c`; manager Passes
155-156. Dependencies: `HP-RG-ORTHO-FN-01` and the Be identity audit.
Disposition: `HP-RG-CUTOFF-FN-01` superseded this default. No source or test
work may use this ID; current maintenance uses active core/ORTHO/CUTOFF-02
authority. Non-goals: compatibility adapters or revival of `1e-8` as a
production value.

### HP-RG-IDTOL-TEST-01 — residual final-identity tolerance validation

Lifecycle: completed historical evidence; superseded. Permission: none.

Owner/canonical: Be identity-tolerance evidence;
[orthogonality/cutoff policy](residual_gaussian_orthogonality_robustness.md).
Evidence: manager Passes 155-156; no dedicated committed test file.
Dependencies: historical `HP-RG-IDTOL-FN-01`.
Disposition: preserve only for interpreting the Be cc-pV5Z/cVDZ audit and Git
history. Current tests use active RG/ORTHO/CUTOFF-02 authority. Non-goals: new
fixtures, source work, or workflow expansion.

### HP-RG-CUTOFF-FN-01 — residual occupation cutoff and identity tolerance defaults

Lifecycle: implemented historical; production cutoff superseded. Permission:
none.

Owner/canonical: historical `5e-8` cutoff and current identity-value origin;
[orthogonality/cutoff policy](residual_gaussian_orthogonality_robustness.md).
Historical source: `src/cartesian_residual_gaussians/residual_basis.jl` and
terminal keyword plumbing. Implementation: `f0b662dca`; manager Passes
171-173. Dependencies: owner-local selection and the Cr marginal-direction
audit. Disposition: `HP-RG-CUTOFF-FN-02` superseded the occupation cutoff;
`identity_atol = 5e-8` remains the current value explicitly preserved by
CUTOFF-02. No source work may use this ID. Non-goals: revival of the `5e-8`
production cutoff or compatibility machinery for old defaults.

### HP-RG-CUTOFF-TEST-01 — residual cutoff/tolerance validation

Lifecycle: completed historical evidence; superseded. Permission: none.

Owner/canonical: historical Cr/Be/H2 `5e-8` validation;
[orthogonality/cutoff policy](residual_gaussian_orthogonality_robustness.md).
Historical test: `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl`;
current assertions no longer use this value. Evidence: manager Passes 171-173.
Dependencies: historical `HP-RG-CUTOFF-FN-01`.
Disposition: preserve for interpreting old results only. Current maintenance
uses `HP-RG-CUTOFF-TEST-02`. Non-goals: restoring stale assertions, adding
fixtures, or granting source/workflow authority.

### HP-RG-CUTOFF-FN-02 — production residual cutoff tightening

Lifecycle: implemented current production policy. Permission: source
maintenance.

Owner/canonical: ordinary production residual cutoff;
[orthogonality/cutoff policy](residual_gaussian_orthogonality_robustness.md).
Source:
`src/cartesian_residual_gaussians/residual_basis.jl`, with compatibility
default plumbing in
`src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl`.
Implementation: `1f7f04e56`; manager Passes 183-184.
Dependencies: owner-local residual selection and current identity policy.
Permission: maintain `residual_occupation_cutoff = 1e-6` and explicit caller
overrides while preserving `identity_atol = 5e-8`, negative-metric checks,
one final merge, and strict orthogonality. Non-goals: cutoff/tolerance changes,
spectral guards, width defaults, flooring, MWG/artifact/driver/solver work, or
Cr2 production behavior.

### HP-RG-CUTOFF-TEST-02 — production residual cutoff validation

Lifecycle: completed evidence; active maintenance. Permission: test
maintenance, not source authority.

Owner/canonical: current production cutoff validation;
[orthogonality/cutoff policy](residual_gaussian_orthogonality_robustness.md).
Test: `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl`.
Evidence: manager Passes 183-185 retain the Be and residual-only Cr2 audits.
Dependencies: `HP-RG-CUTOFF-FN-02` and the existing H2 endpoint.
Permission: maintain in-memory and provenance `1e-6` assertions and current
bounded endpoint checks. Non-goals: new fixtures, source authority, spectral
guards, full HF, Cr2 artifacts/workflow, driver, solver, ECP, or EGOI work.

### HP-RG-NUMCOMP-FN-01 — numerical-complete residual basis and additive consumer

Status: implemented internal opt-in facility at `b2da7070c`.

Permission: source maintenance for the fixed `eta_num = 1e-10` numerical-
complete composition and existing private additive consumer only.

Owner: `CartesianResidualGaussians`, with narrow private additive-reference
composition near the protected ladder owner.

Canonical contract:
[Numerical-complete residual Gaussian basis](numerical_complete_residual_basis.md).

Approved source surfaces:

- `src/cartesian_residual_gaussians/residual_basis.jl` for narrow reuse of the
  existing builder only;
- `src/cartesian_residual_gaussians/augmented_operators.jl` for native
  `[G,R_num]` packet representation;
- `src/cartesian_protected_ladder_bundle.jl` for private in-memory composition;
- `src/cartesian_reference_density/atomic_hf_reference_packets.jl` and
  `src/cartesian_reference_density/screened_hartree_correction.jl` only for
  narrow reuse without contract or result-shape changes.

Permission summary: explicitly call the existing owner-local residual builder
with `eta_num = 1e-10`, injection disabled, and no compactness prefilter; build
exact augmented operators and MWG in native `M=[G,R_num]` order; validate
packet occupied capture after construction; and return the existing in-memory
screened-Hartree correction separately.

Dependencies: implemented ordinary residual basis, exact augmented operators,
MWG, atomic packets, reference-Hartree numerics, and screened-Hartree
correction assembly.

Non-goals: production/default cutoff changes, replacement/injection,
localization, public/driver/artifact/solver work, interaction rotation, EGOI,
Gaussian-array enrichment, or Cr2-specific behavior and claims.

### HP-RG-NUMCOMP-TEST-01 — numerical-complete residual validation

Status: completed bounded validation contract at `b2da7070c`.

Permission: validation maintenance for the existing compact rank/malformed-
metric and bounded H2 surfaces; ignored H2/Be2 and one gated Cr2 comparison
remain measurement-only.

Owner/canonical: numerical-complete residual validation;
[numerical-complete residual Gaussian basis](numerical_complete_residual_basis.md).

Approved surfaces:

- `test/misc/runtests.jl` for compact numerical-rank and malformed-metric
  coverage;
- `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl` for a bounded
  explicit opt-in path without changing ordinary defaults;
- ignored padded H2/Be2 probes and, only after those pass, one ignored Cr2
  fixed-density comparison.

Required gates: owner and final metric spectra, `G-R`/`R-R` identities, packet
capture/trace, broad occupation and `P_GR`, low modes, MWG finiteness/symmetry,
bounded endpoints, and inspected terminal due diligence. Near-threshold metric
failure, occupied capture failure, or a bad H2/Be2 mode stops the lane; do not
floor eigenvalues, inject occupied directions, or weaken thresholds.

Accepted evidence: package load; misc `59/59`; nested H2 augmented `67/67`;
facade `69/69`; reviewed padded Be2 construction, endpoint, and terminal due
diligence. One ignored Cr2 fixed imported-density comparison remains allowed;
it is not a production claim or further source authority.

### HP-RG-SPECTRAL-AUDIT-01 — residual-sector spectral audit

Status: approved measurement-only authority. This is not production source
authority.

Permission: measurement-only through ignored probes and external text/TSV
outputs; no tracked source, test, artifact, or workflow changes.

Owner/canonical: residual-sector numerical diagnostics;
[orthogonality/cutoff policy](residual_gaussian_orthogonality_robustness.md).

Evidence after `HP-RG-CUTOFF-FN-02`: the tightened
`residual_occupation_cutoff = 1.0e-6` performs the intended first cleanup,
dropping Cr2 retained residuals from `68 + 68` to `62 + 62`, but residual-only
spectra still show a low two-owner residual mode:

```text
min eig(K_RR)  =  0.3700413519
min eig(H1_RR) = -7.1647854052
owner weights  = approximately 0.5 / 0.5
```

Approved behavior:

- measurement-only residual-sector audit;
- compute retained residual count by owner;
- compute low eigenvalues of `K_RR`;
- compute low eigenvalues of
  `H1_RR = K_RR + sum_A Z_A U_A_RR`;
- report owner weights for low or flagged eigenvectors;
- report residual occupation composition of low or flagged eigenvectors;
- compare Cr2 residual spectra against one-center atom residual baselines when
  available;
- classify whether low modes are dominated by the smallest retained
  occupations or by otherwise healthy retained modes.

Approved surfaces:

- ignored `tmp/work/*.jl` probes only;
- durable text/TSV output under `/Users/srw/dmrgtmp/...` or CR2 run
  directories.

Forbidden:

- production source changes;
- committed tests or fixtures;
- artifact schema/provenance/reader/manifest changes;
- driver changes;
- MWG/IDA changes;
- dense Vee, full HF, or solver workflow;
- automatic residual pruning;
- kinetic or `H1_RR` spectral-guard implementation;
- cutoff, tolerance, owner grouping, residual-selection, or merge-policy
  changes.

Validation for later audit:

- `git diff --check`;
- package load;
- residual-only audit for one-center Cr atom baseline if available;
- residual-only audit for the current Cr2 fixture;
- no full HF and no new Hamiltonian artifact.

Failure rule: if the audit cannot reconstruct `K_RR`/`H1_RR` cheaply enough
from existing construction seams, stop and report the exact missing reusable
seam. Do not add source instrumentation as part of this lane.

### HP-RG-INJECT-AUDIT-01 — direct-G injection measurement audit

Status: completed historical measurement; no active permission.

Permission: none.

Canonical record:
[Default-off direct-G residual injection](residual_gaussian_injection_hybrid.md).

Evidence: ignored Cr/Cr2 probes and manager running-log history established
that direct-`G` replacement was well conditioned but did not remove the
tested low two-owner residual sector. The audit authorizes no source,
artifact, workflow, or renewed measurement work.

### HP-RG-INJECT-FN-01 — default-off direct-G injection compatibility

Status: implemented preservation-only internal compatibility facility.

Permission: source maintenance for existing default-off behavior only; no
feature expansion or new caller.

Owner:
`CartesianResidualGaussians`, with narrow terminal residual-GTO compatibility
plumbing.

Canonical contract:
[Default-off direct-G residual injection](residual_gaussian_injection_hybrid.md).

Approved source surfaces:

- `src/cartesian_residual_gaussians/residual_basis.jl`;
- `src/cartesian_residual_gaussians/augmented_operators.jl`;
- `src/cartesian_residual_gaussians/mwg_interaction.jl`;
- `src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl` only
  for existing private keyword/validation compatibility.

Permission summary: when explicitly enabled in memory, classify owner-local
principal supplement modes, globally clean the injected target, replace a
subspace of `G` with `F = [Y, G Q_perp]`, transform exact one-body
operators into `[F,R]`, and reserve residual MWG/IDA channels for true
residual directions. Disabled behavior must remain unchanged.

Dependencies: ordinary owner-local RG selection and exact augmented operator
blocks. This path is not a dependency of the current compact-main protected
builder.

Non-goals: default-on behavior, feature expansion, public/driver wiring,
injection-enabled artifacts, global residual selection, MWG descriptors for
injected directions, spectral pruning, solver/HF, or Cr2 production claims.
### HP-RG-PROTECT-INJECT-DESIGN-01 — protected-original compact-main injection design

Status: completed design rationale; not source authority by itself. The
implemented source lanes below realize the accepted contract.

Permission: design-only rationale maintenance; no source, test, artifact, or
workflow execution authority.

Owner and canonical contract:

- `CartesianResidualGaussians` protected geometry;
- [Protected-localized basis convention](protected_localized_basis.md).

Permission summary: define compact-first replacement over
`M = [G, R_compact]`, distinguish Gaussian Gram cleanup from representability,
preserve protected original spans, and treat unsupported broad directions as
basis-insufficiency diagnostics rather than MWG residual channels.

Dependencies: the ordered compact-first residual basis contract and exact
supplement overlap/mixed-overlap inputs.

Non-goals: source authority under this ID, direct historical `G` injection,
public input/defaults, artifacts, solver/HF, residual-policy changes, or Cr2
production claims.

### HP-RG-OCC-FIRST-INJECT-AUDIT-01 — occupied-first global injection measurement audit

Status: completed historical measurement record; not active source authority.

Permission: none.

Owner/canonical:
[occupied-first injection geometry](occupied_first_injection.md).

Purpose: established mandatory occupied protection and capture-spectrum
selection on Be/Ne before the source-backed helper was approved. Evidence is
retained in manager running-log Passes 323-324 and summarized by the canonical
contract below.

### HP-RG-OCC-FIRST-INJECT-FN-01 — occupied-first injection geometry

Lifecycle: implemented. Permission: source maintenance.

Owner/canonical: `CartesianResidualGaussians`;
[occupied-first injection geometry](occupied_first_injection.md).

Source: `src/cartesian_residual_gaussians/residual_basis.jl`, with read-only
packet/import coefficient inputs under their existing owners.

Permission: validate physical capture geometry, make supplied `Y_occ`
mandatory, keep pre-inclusion capture distinct from post-inclusion recovery,
and capture-select optional supplement directions. Weak rejected directions
never become MWG residual channels.

Boundary: the helper remains unwired into the protected builder and is not a
direct substitution for staged protected geometry over `M=[G,R_compact]`.

### HP-RG-OCC-FIRST-INJECT-TEST-01 — occupied-first validation

Lifecycle: completed validation contract. Permission: validation maintenance.

Owner/canonical: occupied-first geometry;
[occupied-first injection geometry](occupied_first_injection.md).

Tests: `test/misc/runtests.jl` and
`test/nested/cartesian_occupied_first_injection_runtests.jl`.

Permission: maintain synthetic pre/post and malformed-capture checks plus the
bounded real packet-driven Be/Ne PQS gate and terminal due diligence.

Non-goals: protected composition under these IDs, screened Hartree, EGOI,
shell-local injection, artifacts, public/solver workflow, exchange, or Cr2.

### HP-RG-PROTECT-ADDREF-FN-01 - protected additive atomic reference correction

Lifecycle: implemented by `0b778a676`. Permission: source maintenance.

Owner/canonical: protected member composition and existing reference owners;
[protected additive reference correction](protected_additive_reference_correction.md).

Source:

- `src/cartesian_residual_gaussians/residual_basis.jl`;
- `src/cartesian_residual_gaussians/augmented_operators.jl`;
- `src/cartesian_gaussian_raw_blocks/mixed_hartree_blocks.jl`;
- `src/cartesian_reference_density/atomic_hf_reference_packets.jl`;
- `src/cartesian_reference_density/screened_hartree_correction.jl`;
- private composition in `src/cartesian_protected_ladder_bundle.jl`.

Permission: build compact `R` once; use staged protected geometry with the
full-rank occupied union mandatory for basis protection; preserve original
per-packet occupied blocks for additive `P0`; build placed fitted-potential
`GG/GA/AA`; include all self and twice-cross `E0` terms; transform `J0`
through native protected/localized one-body operators; and return the existing
in-memory `ScreenedHartreeCorrection` plus reference diagnostics. The private
seam returns `(member, correction, reference)`; the no-reference path remains
unchanged.

Hard boundaries: packet self-integrity and structural mapping remain exact;
mapped overlap uses the existing `1e-10` numerical gate; mandatory occupied
capture failures stop; fitted-potential total/self/cross consistency is
reported rather than forced below `1e-8 Ha`; retired polished packets reject.

Non-goals: public input, corrected artifacts, protected atoms, counterpoise,
compact/high transfer, `Vee` transforms, `C' V C`, solvers, EGOI, exchange,
mapping/default changes, or Cr2 production claims.

### HP-RG-PROTECT-ADDREF-TEST-01 - additive-reference validation

Lifecycle: completed validation contract. Permission: validation maintenance.

Owner/canonical: protected additive composition;
[protected additive reference correction](protected_additive_reference_correction.md).

Tests: `test/misc/runtests.jl` and
`test/nested/cartesian_screened_hartree_correction_runtests.jl`.

Supporting atomic-packet evidence remains separately owned by
`HP-PQS-ATOMREF-PACKET-TEST-01`; this ID does not grant maintenance of that
test surface.

Permission: maintain mandatory-union recovery, packet embedding/failure,
per-packet trace, additive `P0/q0`, self/cross `E0`, placed raw-block,
protected/localized `J0`, correction-anchor, no-reference parity, and
ordinary-packet rejection/diagnostic checks. The accepted padded Be2 smoke is
structural evidence only; its retired polish-assisted energy value is not a
current endpoint gate.

### HP-RG-PROTECT-INJECT-FN-01 — staged protected-original geometry

Lifecycle: implemented internal/default-off. Permission: source maintenance.

Owner/canonical: `CartesianResidualGaussians`;
[protected-localized basis](protected_localized_basis.md).

Source: `src/cartesian_residual_gaussians/residual_basis.jl`.

Permission: consume the already-built compact residual object; build protected
and broad original subspaces; keep Gaussian Gram, representability, and
fake-RDM gates distinct; and return transform-ready `Z`, `B`, `Q_perp`, `F`,
and diagnostics. Rejected broad directions never become MWG channels.

### HP-RG-PROTECT-INJECT-TEST-01 — protected geometry validation

Lifecycle: completed validation contract. Permission: validation maintenance.

Owner/canonical: protected geometry;
[protected-localized basis](protected_localized_basis.md).

Evidence: source-backed staged-geometry probes and
`docs/src/developer/reports/cr2_staged_subspace_filter_870498b54/README.md`;
no dedicated committed unit-test file.

Non-goals: public/default workflow, artifacts, operator/interaction work under
these IDs, solvers, selection-policy changes, or Cr2 production claims.

### HP-RG-PROTECT-ONEBODY-AUDIT-01 — protected fixed-sector one-body audit

Status: completed historical measurement; not active source authority.

Permission: none.

Owner/canonical: protected exact one-body history;
[protected-localized basis](protected_localized_basis.md).

Evidence: manager running-log Passes 254-255 and
`docs/src/developer/reports/cr2_protected_onebody_audit_eaf05a38c/README.md`.
The audit established the dataflow later implemented under the source/test IDs
below.

### HP-RG-PROTECT-ONEBODY-FN-01 — exact protected one-body transform

Lifecycle: implemented internal/default-off. Permission: source maintenance.

Owner/canonical: `CartesianResidualGaussians`;
[protected-localized basis](protected_localized_basis.md).

Source: `src/cartesian_residual_gaussians/augmented_operators.jl`, with
transform-ready geometry from `residual_basis.jl`.

Permission: construct dense exact fixed-sector kinetic, per-center unit
nuclear, and assembled `H1_F` matrices through the actual protected/localized
one-body transform, with orthogonality and symmetry diagnostics.

### HP-RG-PROTECT-ONEBODY-TEST-01 — protected one-body validation

Lifecycle: completed validation contract. Permission: validation maintenance.

Owner/canonical: protected exact one-body operators;
[protected-localized basis](protected_localized_basis.md).

Evidence: source-backed dense replay and
`docs/src/developer/reports/cr2_protected_onebody_audit_eaf05a38c/README.md`;
no dedicated committed unit-test file.

Non-goals: matrix-action frameworks, public/default workflow, artifacts,
interaction rotation, solvers, residual policy, screened references, or Cr2
production claims.

### HP-RG-PROTECT-VEE-AUDIT-01 — protected interaction decision audit

Status: completed historical measurement; not source authority.

Permission: none.

Owner and canonical outcome:

- protected interaction measurement record;
- [Protected-localized basis convention](protected_localized_basis.md).

Validation/evidence surfaces:

- `tmp/work/cr2_protected_vee_algebra_debug.jl`;
- `tmp/work/protected_localized_injection_vee_probe.jl`;
- manager running-log Passes 269-270.

Recorded outcome: direct `C' V C` interaction rotation is invalid and must not
be reused. The viable convention is localized `L`, exact `H1_L`, and the
inherited pre-injection site-order `Vee_M` as `Vee_L`, with explicit
orthogonality, localization, finite/symmetry, and bounded physics checks.

Dependencies: the implemented protected geometry and exact one-body source
lanes above.

Non-goals: renewed measurement authority, tracked source authority under this
ID, alternative interaction rotations, public workflow, artifacts, solver,
screened-reference changes, or Cr2 production claims.

### HP-RG-PROTECT-ART-FN-01 — protected-localized Hamiltonian artifact variant

Status: implemented.

Owner: `src/cartesian_ida_hamiltonian.jl`, with artifact inputs supplied by
`src/cartesian_residual_gaussians/augmented_operators.jl`.

Canonical contract:
[Protected-localized artifact contract](protected_localized_artifact.md).

Source paths:

- `src/cartesian_ida_hamiltonian.jl`;
- `src/cartesian_residual_gaussians/augmented_operators.jl`.

Evidence: bounded protected write/readback and rejection smokes in
`docs/src/developer/pqs_manager_running_log.md` Pass 299 and implementation
commit `fd105b751`. The ordinary IDA artifact test is not protected-artifact
coverage.

Dependencies: the implemented
[protected-localized basis convention](protected_localized_basis.md) and its
already-computed native `H1_L`, inherited-site `Vee_L`, sectors, diagnostics,
and provenance.

Permission: write and read the recognized, versioned, opt-in
protected-localized artifact without changing its native matrix order.

Exclusions: ordinary artifact semantics, public/default driver workflow,
solver behavior, RG/injection selection, EGOI, rho0 or screened-Hartree,
alternative interaction transforms, and Cr2-specific production behavior.

### HP-RG-PROTECT-ART-TEST-01 — protected artifact validation

Status: implemented validation contract.

Owner and canonical contract:
`src/cartesian_ida_hamiltonian.jl` and
[Protected-localized artifact contract](protected_localized_artifact.md).

Evidence: manager Pass 299 and commit `fd105b751` record accepted bounded
write/readback and rejection smokes. No dedicated committed core
protected-artifact test exists. The external-GTO integration test constructs a
synthetic protected artifact for its separately owned sidecar identity checks;
it is not core protected-artifact coverage.

Dependencies: `HP-RG-PROTECT-ART-FN-01` and the ordinary reader's strict
artifact-kind boundary.

Permission: validate matrix/count/metadata roundtrip and reject wrong or
missing protected artifact identity and convention data.

Exclusions: solver, endpoint, public workflow, rho0, production-default, or
converged Cr2 claims.

### HP-RG-PROTECT-ARTLOC-FN-01 — protected artifact row-locality metadata

Status: implemented.

Owner: `protected_localized_row_locality` in
`src/cartesian_residual_gaussians/augmented_operators.jl`, with validation and
persistence in `src/cartesian_ida_hamiltonian.jl`.

Canonical contract:
[Protected-localized artifact contract](protected_localized_artifact.md).

Source paths:

- `src/cartesian_residual_gaussians/augmented_operators.jl`;
- `src/cartesian_ida_hamiltonian.jl`.

Evidence: bounded native-locality smokes in manager Pass 301 and implementation
commit `3fe2af697`. The ordinary IDA artifact test does not cover this metadata.

Dependencies: `HP-RG-PROTECT-ART-FN-01`, native `M` position operators, and
the implemented protected-localized transform.

Permission: attach validated native-order center, sector, inverse-permutation,
and optional all-or-none spread metadata to the protected artifact.

Exclusions: matrix or range reordering, label-derived centers, new
second-moment construction, artifact consumers, EGOI, additive-reference,
rho0 or screened-Hartree, solver/public-driver behavior, and Cr2 production
claims.

### HP-RG-PROTECT-ARTLOC-TEST-01 — row-locality validation

Status: implemented validation contract.

Owner and canonical contract:
`src/cartesian_residual_gaussians/augmented_operators.jl`,
`src/cartesian_ida_hamiltonian.jl`, and
[Protected-localized artifact contract](protected_localized_artifact.md).

Evidence: manager Pass 301 and commit `3fe2af697` record the accepted center,
inverse-permutation, sector, spread, and legacy-no-locality smokes. No
dedicated committed row-locality test file exists.

Dependencies: `HP-RG-PROTECT-ARTLOC-FN-01` and native-order artifact
roundtrip under `HP-RG-PROTECT-ART-TEST-01`.

Permission: validate native centers and sectors, deterministic inverse
permutations, optional all-or-none spreads, and legacy
`row_locality = nothing` readback.

Exclusions: z-sorted canonical matrices, solver/public workflow, rho0,
production defaults, or converged Cr2 claims.

### HP-RG-PROTECT-EGOI-AUDIT-01 - protected-localized EGOI measurement audit

Lifecycle: completed historical measurement; not active source, test,
artifact, or workflow authority.

Permission: none.

Owner and evidence:

- [Retained-GTO local-product EGOI](retained_gto_egoi.md);
- manager running-log Passes 302-308.

The completed audit established the retained-original `s1+s2`, local-product,
`M2` candidate. Machine-local probes and outputs remain historical evidence.

### HP-RG-PROTECT-EGOI-FN-01 - retained-GTO local-product EGOI helper

Lifecycle: approved pending internal source authority; not implemented in
committed source.

Owner:

- [Retained-GTO local-product EGOI](retained_gto_egoi.md).

Approved source surfaces:

- `src/hamiltonian_corrections.jl` as primary owner;
- `src/cartesian_residual_gaussians/residual_basis.jl` and
  `src/cartesian_residual_gaussians/augmented_operators.jl` only for retained
  source mapping or transform-ready `Qtarget`.

Dependencies:

- committed generic matrix-level EGOI routines;
- [protected-localized basis convention](protected_localized_basis.md);
- native protected-localized row geometry for the `M2` mask.

Permission: build the owner-balanced retained-original `s1+s2` target, native
`Qtarget`, local symmetric products, `AA-AA` / `BB-BB` / `AA-BB` acceptance
blocks, exactly local `M2` `DeltaV`, and compact diagnostics in memory.

Exclusions: public APIs/workflows, artifacts, solver/HF/MP2-NO integration,
selection changes, rho0/screened-Hartree, broad protected targets, AB overlap
products, `s3`/`p`/`d` promotion, and Cr2 production claims.

Uncommitted `src/hamiltonian_corrections.jl` additions do not satisfy or
change this lifecycle.

### HP-RG-PROTECT-EGOI-TEST-01 - retained-GTO EGOI validation

Lifecycle: approved pending validation contract; no committed protected
retained-GTO helper test exists.

Owner:

- [Retained-GTO local-product EGOI](retained_gto_egoi.md).

Permission: bounded H/Be/Be2 smokes plus an ignored Cr2 replay, with target
projection, rank, block residual, `DeltaV`, saturation, symmetry, and low-Fock
diagnostics. `max_disallowed_delta_v` must equal zero exactly.

Exclusions: production Cr2 HF, large committed fixtures, artifact/workflow
tests, and validation of uncommitted WIP.

### HP-RG-PROTECT-LADDER-XFER-AUDIT-01 - same-parent ladder transfer audit

Lifecycle: completed historical measurement; not active source, artifact, or
workflow authority.

Permission: none.

Owner and evidence:

- [Protected-localized ladder bundles](protected_localized_ladder.md);
- manager running-log ladder audit history.

The audit established exact cross-overlap-only transfer and target-member
Hamiltonian evaluation. Its reusable outcome is the implemented bundle
facility below.

### HP-RG-PROTECT-LADDER-BUNDLE-FN-01 - protected-localized ladder bundle facility

Lifecycle: implemented internal opt-in facility.

Owner:

- [Protected-localized ladder bundles](protected_localized_ladder.md).

Implemented source surfaces:

- `src/cartesian_protected_ladder_bundle.jl` as primary owner;
- `src/GaussletBases.jl` for its include;
- existing related owners `src/cartesian_representation_transfer.jl` and
  `src/cartesian_ida_hamiltonian.jl` without changing their contracts.

Dependencies:

- [protected-localized basis convention](protected_localized_basis.md);
- [protected-localized artifact contract](protected_localized_artifact.md);
- shared parent lattice, identical supplement, and one Coulomb expansion.

Permission: write a versioned directory manifest, protected member artifacts,
exact adjacent `S_BA` sidecars, optional native-order restart sidecars, and
bounded summaries; transfer only as `C_B = S_BA * C_A` and evaluate with
target `H1_L` / `Vee_L`.

Exclusions: generalized overlap, source-Hamiltonian or source-`Vee`
transforms, interaction rotation, ladder-owned representation sidecars,
protected member schema changes, solver/UHF, EGOI, screened-Hartree/rho0,
public workflow, and Cr2 production claims. The separately owned external-GTO
sidecar does not change this permission.

### HP-RG-PROTECT-LADDER-BUNDLE-TEST-01 - ladder bundle validation

Lifecycle: implemented validation contract.

Owner:

- [Protected-localized ladder bundles](protected_localized_ladder.md).

Validation evidence: source commit `3eaa812a9`, manager running-log Pass 313,
a bounded H2 bundle/readback smoke, and an ignored Cr2 `ns=7 -> ns=9` replay.

Permission: validate identity/layout, member readback, shared-parent proof,
same supplement/Coulomb policy, cross-overlap shape and singular values,
native restart ordering, transferred traces/orthogonality, fixed-density
target energy, provenance, and summary roundtrip.

Exclusions: committed large Cr2 fixtures, production HF, solver continuation,
and alternate artifact or transfer conventions.

### HP-RG-RHO0-GAL-AUDIT-01 - row-gauge rho0 audit

Status: completed and superseded historical measurement evidence.

Permission: none.

Owner and evidence:
[rho0 and reference-density correction history](rho0_reference_density_matrix.md).
No source or test authority remains.

### HP-RHO0-REFDENS-AUDIT-01 - fixed-P0 audit

Status: completed and superseded historical measurement evidence.

Permission: none.

Owner and evidence:
[rho0 and reference-density correction history](rho0_reference_density_matrix.md).

### HP-RHO0-REFDENS-FN-01 - historical candidate correction owner

Lifecycle: unapproved historical planning name. Permission: none.

Owner/history:
[rho0 and reference-density correction history](rho0_reference_density_matrix.md).

No source surface, caller, or implementation authority was approved. The
historical proposal is retained only in
[rho0 and reference-density correction history](rho0_reference_density_matrix.md).

### HP-RHO0-REFDENS-ERI-01 - historical candidate mixed-ERI owner

Lifecycle: unapproved historical planning name. Permission: none.

Owner/history:
[rho0 and reference-density correction history](rho0_reference_density_matrix.md).

No source surface, kernel, test, or implementation authority was approved.
Durable neutral mixed-Hartree numerics are governed by the implemented MIXH
families below, not by this candidate name.

### HP-RHO0-REFDENS-MIXH-AUDIT-01 - exact mixed-Hartree seam audit

Status: completed historical measurement evidence. Its durable result is
the implemented neutral MIXH/FEXACT family below.

Permission: none.

Owner/history:
[rho0 and reference-density correction history](rho0_reference_density_matrix.md).

### HP-RHO0-MIXH-GG-FN-01 - exact mixed-Hartree GG block

Status: implemented.

Permission: source maintenance for exact one-center finite symmetric-`P_A`
`GG` construction and compact diagnostics only.

Owner:
`CartesianGaussianRawBlocks`.

Canonical contract:
[reference Hartree numerics](reference_hartree_numerics.md).

Source:

- `src/cartesian_gaussian_raw_blocks/mixed_hartree_blocks.jl`;
- `src/cartesian_gaussian_raw_blocks/CartesianGaussianRawBlocks.jl` for
  internal module wiring;
- `src/gaussian_coulomb_reference.jl` only for the existing dense oracle or
  narrow pair-term algebra reuse.

Dependencies and permission: one-center finite symmetric `P_A`, reference
Gaussian supplement data, the existing Coulomb Gaussian expansion, and
terminal Gaussian-sum accumulation may produce exact finite symmetric
terminal/base `GG` plus compact diagnostics. No correction policy,
exchange, artifact, public workflow, solver, residual selection, or Cr/Cr2
authority is implied.

### HP-RHO0-MIXH-GG-TEST-01 - exact mixed-Hartree GG validation

Status: implemented validation authority.

Permission: validation maintenance for the stated compact gate and existing
consumer test only.

Owner/canonical: `CartesianGaussianRawBlocks`;
[reference Hartree numerics](reference_hartree_numerics.md).

Committed consumer test:
`test/nested/cartesian_screened_hartree_correction_runtests.jl`.

Compact gate: package load; bounded H/Be/Be2 finite/symmetric output;
angular and off-diagonal same-center pair coverage; dense-oracle spots; no
Cr/Cr2. Historical acceptance is manager Pass 284 and commit `efaee93f6`.

### HP-RHO0-MIXH-GAAA-FN-01 - exact mixed-Hartree GA/AA blocks

Status: implemented.

Permission: source maintenance for exact `GA`, symmetric `AA`, and compact
diagnostics only.

Owner:
`CartesianGaussianRawBlocks`.

Canonical contract:
[reference Hartree numerics](reference_hartree_numerics.md).

Source:

- `src/cartesian_gaussian_raw_blocks/mixed_hartree_blocks.jl`;
- `src/cartesian_gaussian_raw_blocks/CartesianGaussianRawBlocks.jl` for
  internal module wiring.

Dependencies and permission: reuse the validated one-center `P_A` pair
stream and Coulomb-expanded factors to form exact `GA` and symmetric `AA`,
including angular reference pairs and supplement rows. The wrapper remains
a neutral oracle seam. It grants no correction, exchange, artifact, public
workflow, solver, residual-selection, or Cr/Cr2 authority.

### HP-RHO0-MIXH-GAAA-TEST-01 - exact mixed-Hartree GA/AA validation

Status: implemented validation authority.

Permission: validation maintenance for the accepted bounded parity, angular,
dimension, and oracle gate only; no dedicated committed test file.

Owner/canonical: `CartesianGaussianRawBlocks`;
[reference Hartree numerics](reference_hartree_numerics.md).

Test ownership: no dedicated committed test file. The accepted compact
ignored-probe gate covered prior `GG` parity, finite `GA`, symmetric `AA`,
angular cases, and dense-oracle spots. Historical acceptance is manager
Pass 286 and commit `daac231d0`.

### HP-RHO0-MIXH-FEXACT-FN-01 - protected exact-Hartree transform

Status: implemented.

Permission: source maintenance for the exact protected fixed-sector transform,
localized `sym(W' * J0_F * W)`, and existing convenience composition only.

Owner:
`CartesianResidualGaussians`.

Canonical contract:
[reference Hartree numerics](reference_hartree_numerics.md).

Source:
`src/cartesian_residual_gaussians/augmented_operators.jl`.

Dependencies and permission: consume exact `GG/GA/AA`, existing protected
geometry `F = [Z, M Q_perp]`, and an existing localization `W` to form the
exact fixed-sector operator and `sym(W' * J0_F * W)`. It grants no geometry
selection, interaction transform, correction, artifact, public workflow,
solver, exchange, or Cr/Cr2 authority.

### HP-RHO0-MIXH-FEXACT-TEST-01 - protected exact-Hartree validation

Status: implemented validation authority.

Permission: validation maintenance for bounded H/Be/Be2 raw, fixed-sector,
localized, and oracle checks only; no dedicated committed test file.

Owner/canonical: `CartesianResidualGaussians`;
[reference Hartree numerics](reference_hartree_numerics.md).

Test ownership: no dedicated committed test file. The accepted compact
ignored-probe gate covered H/Be/Be2 raw-block replay, fixed/localized finite
symmetry, dimension and protected-geometry checks, and dense-oracle spots.
Historical acceptance is manager Pass 288 and commit `40a6f7e99`.

### HP-RHO0-FAPP-AUDIT-01 - approximate fixed-P0 Fock audit

Status: completed historical measurement evidence.

Permission: none.

Owner/history:
[rho0 and reference-density correction history](rho0_reference_density_matrix.md).

### HP-RHO0-FAPP-FN-01 - approximate IDA energy/Fock seam

Lifecycle: implemented but caller-free and dormant. Permission: none.

Owner/source: `src/cartesian_ida_hamiltonian.jl`.

Canonical lifecycle:
[rho0 and reference-density correction history](rho0_reference_density_matrix.md).
No new caller, source work, correction policy, or public surface is authorized.

### HP-RHO0-FAPP-TEST-01 - approximate IDA derivative validation

Lifecycle: completed historical validation. Permission: none.

Owner/history:
[rho0 and reference-density correction history](rho0_reference_density_matrix.md).

Evidence: ignored finite-difference gates recorded in the rho0 history; no
dedicated committed test or live caller remains.

### HP-RHO0-ANCHOR-FN-01 - old full-interaction anchor

Lifecycle: superseded. Permission: none.

Owner/history:
[rho0 and reference-density correction history](rho0_reference_density_matrix.md).

The old `Delta_F0_alpha/beta` interpretation is not a Hartree correction
contract and must not be revived.

### HP-RHO0-ANCHOR-TEST-01 - old anchor validation

Lifecycle: superseded historical evidence. Permission: none.

Owner/history:
[rho0 and reference-density correction history](rho0_reference_density_matrix.md).

### HP-RHO0-CORR-AUDIT-01 - corrected-Hamiltonian audit

Status: completed and superseded historical measurement evidence.

Permission: none.

Owner/history:
[rho0 and reference-density correction history](rho0_reference_density_matrix.md).

### HP-RHO0-JANCHOR-FN-01 - direct-Hartree anchor helper

Lifecycle: implemented but superseded in use and dormant. Permission: none.

Owner/history:
[rho0 and reference-density correction history](rho0_reference_density_matrix.md).

Source:

- `src/cartesian_ida_hamiltonian.jl`;
- `src/cartesian_residual_gaussians/augmented_operators.jl`.

Current replacement:
`src/cartesian_reference_density/screened_hartree_correction.jl`, governed
by [screened Hartree correction assembly](screened_hartree_correction_assembly.md)
and [screened Hartree residual density](screened_hartree_residual_density.md).

No new caller, source work, artifact, public workflow, solver, exchange, or
Cr/Cr2 work is authorized.

### HP-RHO0-JANCHOR-TEST-01 - direct-Hartree anchor validation

Lifecycle: completed historical validation. Permission: none.

Owner/history:
[rho0 and reference-density correction history](rho0_reference_density_matrix.md).

Evidence: ignored direct-anchor and finite-difference gates recorded in the
rho0 history. The live equivalent formula is validated under screened-Hartree
correction authority.

### HP-RHO0-XPAIR-AUDIT-01 - exchange/direct pairing question

Status: approved but deferred measurement-only question.

Permission: measurement-only after explicit design-manager reactivation, for
ignored H/Be/Be2 diagnostics only.

Owner and evidence:
[rho0 and reference-density correction history](rho0_reference_density_matrix.md).
It is not a current blocker or source lane and authorizes no tracked source,
test, artifact, public workflow, solver, exchange implementation, or Cr/Cr2
work.

## Implemented Cartesian Gaussian Raw-Block Nuclear Owner

### HP-CGRB-FILE-01 — neutral Cartesian Gaussian raw-block module files

Lifecycle: implemented. Permission: source maintenance.

Owner/canonical: `CartesianGaussianRawBlocks`;
[nuclear raw blocks](cartesian_gaussian_raw_blocks_nuclear.md).

Source:

```text
src/cartesian_gaussian_raw_blocks/CartesianGaussianRawBlocks.jl
src/cartesian_gaussian_raw_blocks/nuclear_blocks.jl
src/GaussletBases.jl
```

The `GaussletBases.jl` surface is internal include wiring only.

Validation: `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl` and
`test/core/runtests.jl`.

Dependencies: shared Gaussian analytic integrals, Qiu-White proxy primitives,
and the producer-resolved Coulomb expansion.

Permission: maintain the internal module and nuclear owner; no public export.

Non-goals: mixed-Hartree/non-nuclear expansion, routes, caches, artifacts,
drivers, or solver/public API.

### HP-CGRB-FN-01 — exact uncharged Gaussian nuclear raw blocks

Lifecycle: implemented. Permission: source maintenance.

Owner/canonical: neutral nuclear kernel;
[nuclear raw blocks](cartesian_gaussian_raw_blocks_nuclear.md).

Source:

```text
src/cartesian_gaussian_raw_blocks/nuclear_blocks.jl
```

Entry point: `gaussian_nuclear_raw_blocks_by_center(...)`.

Validation/evidence: H2/core tests above and manager-log Pass 078.

Dependencies: `HP-CGRB-FILE-01` and producer-wide Coulomb expansion parity.

Permission: maintain exact uncharged by-center `G-A`/`A-A` output and the
stable pairwise analytic factor formula.

Non-goals: physical charge application, terminal/residual transforms,
non-nuclear blocks, Hamiltonian assembly, or persistent providers.

### HP-CGRB-FN-02 — nuclear one-dimensional axis-family reuse

Lifecycle: implemented. Permission: source maintenance.

Owner/canonical: neutral nuclear kernel;
[nuclear raw blocks](cartesian_gaussian_raw_blocks_nuclear.md).

Source:

```text
src/cartesian_gaussian_raw_blocks/nuclear_blocks.jl
```

Validation/evidence: H2/core tests and manager-log Pass 082B.

Dependencies: `HP-CGRB-FN-01`, supplement axis-family inventory, and one
resolved Coulomb expansion.

Permission: maintain function-local family maps, unique center/family-pair
tables, orientation flags, term-first filling, and coupled primitive products.

Non-goals: independent x/y/z contraction, persistent caches, metadata,
non-nuclear work, route semantics, public/artifact/Cr2 workflow.

### HP-CGAI-FN-01 — optional Cartesian Gaussian axis helper

Lifecycle: unused optional authority; superseded as a performance endpoint.
Permission: dormant optional source support only.

Owner/canonical: nuclear axis-helper disposition;
[nuclear raw blocks](cartesian_gaussian_raw_blocks_nuclear.md).

Potential source surface:

```text
src/cartesian_gaussian_axis_integrals.jl
```

Implemented object/caller: none. The proposed in-place table helper never
landed; nuclear family reuse succeeded inside `nuclear_blocks.jl`.

Dependencies: only a future demonstrated need from `HP-CGRB-FN-02`.

Permission: no active work. A separately assigned pass may use this ID only
for the originally approved nonallocating nuclear factor/table support.

Non-goals: treating the later live allocating generic axis helpers as this
implementation, broad axis/raw-block rewrites, public API, caches, or route,
artifact, residual, Qiu-White, and Cr2 behavior.

### HP-CGRB-WIRE-01 — Residual Gaussian and Qiu-White rewiring

Lifecycle: implemented. Permission: caller maintenance.

Owner/canonical: neutral nuclear consumer boundary;
[nuclear raw blocks](cartesian_gaussian_raw_blocks_nuclear.md).

Live caller surfaces:

```text
src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl
src/ordinary_qw_raw_blocks.jl
src/ordinary_qw_operator_assembly.jl
```

Validation/evidence: H2 endpoint and manager-log Pass 078 Qiu-White/Be2/Cr2
parity.

Dependencies: `HP-CGRB-FN-01`; downstream residual/Qiu-White contracts remain
separate owners.

Permission: maintain the two direct neutral-kernel call sites and indirect
Qiu-White assembly compatibility.

Non-goals: residual selection/transforms, Qiu-White route objects, terminal or
parent construction, reports, artifacts, public API, or Cr2 workflow.

### HP-CGRB-TEST-01 — nuclear extraction validation

Lifecycle: validation completed; tracked coverage is indirect. Permission:
test maintenance.

Owner/canonical: neutral nuclear validation;
[nuclear raw blocks](cartesian_gaussian_raw_blocks_nuclear.md).

Test/evidence:

```text
test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl
test/core/runtests.jl
docs/src/developer/pqs_manager_running_log.md (Passes 078-082B)
```

Dependencies: implemented kernel/wiring and accepted ignored Qiu-White, Be2,
and Cr2 parity evidence.

Permission: maintain the H2 endpoint and stable-factor oracle. No dedicated
raw-block test file currently exists.

Non-goals: Cr2 workflow, artifact/public API, route/status/payload, or residual
internal-vocabulary tests.

## Implemented Cartesian Gaussian Raw-Block Non-Nuclear Owner

### HP-CGRB-NN-FILE-01 — non-nuclear raw-block file

Lifecycle: implemented. Permission: source maintenance.

Owner/canonical: `CartesianGaussianRawBlocks`;
[non-nuclear raw blocks](cartesian_gaussian_raw_blocks_non_nuclear.md).

Source:

```text
src/cartesian_gaussian_raw_blocks/non_nuclear_blocks.jl
src/cartesian_gaussian_raw_blocks/CartesianGaussianRawBlocks.jl
```

Validation: `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl`.

Dependencies: shared private axis-family and analytic axis-table helpers.

Permission: maintain the internal non-nuclear owner and module include.

Non-goals: root/public exports, nuclear/mixed-Hartree expansion, route/cache,
artifact, driver, or solver ownership.

### HP-CGRB-NN-FN-01 — exact non-nuclear Gaussian raw blocks

Lifecycle: implemented. Permission: source maintenance.

Owner/canonical: neutral non-nuclear kernel;
[non-nuclear raw blocks](cartesian_gaussian_raw_blocks_non_nuclear.md).

Source:

```text
src/cartesian_gaussian_raw_blocks/non_nuclear_blocks.jl
src/cartesian_gaussian_axis_integrals.jl
```

The axis-integral file supplies shared private analytic tables.

Entry points: `gaussian_non_nuclear_raw_blocks(...)` and
`gaussian_non_nuclear_overlap_blocks(...)`.

Validation/evidence: H2 endpoint and manager-log Passes 084B-086B.

Dependencies: `HP-CGRB-NN-FILE-01`, deterministic supplement order, and
shared axis-family/table conventions.

Permission: maintain exact `G-A`/`A-A` overlap, kinetic, x/y/z, and x2/y2/z2
blocks, symmetry, family reuse, and the overlap-only path.

Non-goals: status/provider/cache objects, nuclear/mixed-Hartree blocks,
terminal `G-G`, residual transforms, routes, artifacts, or public API.

### HP-CGRB-NN-WIRE-01 — Residual Gaussian and Qiu-White non-nuclear rewiring

Lifecycle: implemented. Permission: caller maintenance.

Owner/canonical: neutral non-nuclear kernel;
[non-nuclear raw blocks](cartesian_gaussian_raw_blocks_non_nuclear.md).

Source/callers:

```text
src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl
src/ordinary_qw_raw_blocks.jl
src/ordinary_qw_operator_assembly.jl
```

Validation: `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl`
and manager-log Passes 084B-086B.

Dependencies: `HP-CGRB-NN-FN-01`, Residual Gaussian exact operators, and the
live Qiu-White operator-assembly path.

Permission: maintain the Residual Gaussian and main diatomic Qiu-White callers
of the neutral blocks.

Non-goals: rewiring or deleting the retained atomic-reference, factor-term,
hybrid-sidecar, dense-parent-probe, or CPB/provider Qiu-White helpers; changing
residual transforms, terminal `G-G`, routes, caches, artifacts, or public API.

### HP-CGRB-NN-TEST-01 — non-nuclear extraction validation

Lifecycle: completed validation. Permission: test maintenance.

Owner/canonical: neutral non-nuclear kernel;
[non-nuclear raw blocks](cartesian_gaussian_raw_blocks_non_nuclear.md).

Test/evidence:

```text
test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl
docs/src/developer/pqs_manager_running_log.md (Passes 084B-086B)
```

Dependencies: implemented kernel/wiring and accepted ignored Qiu-White, Be2,
and Cr2 parity evidence.

Permission: maintain the H2 endpoint and existing indirect parity coverage. No
dedicated raw-block test file currently exists.

Non-goals: Cr2 workflow, artifacts, public API, route/status/payload, or
residual internal-vocabulary tests.

## Implemented R3 Exact-Operator Optimizations

### HP-R3GG-FN-01 — terminal G-G product matrices

Lifecycle: implemented. Permission: source maintenance.

Owner/canonical: `CartesianFinalBasisRealization`;
[terminal G-G products](r3_terminal_gg_product_matrices.md).

Source:

```text
src/cartesian_final_basis_realization/pqs_terminal_one_body.jl
src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl
```

Dependencies: realized terminal blocks and Residual Gaussian exact transforms.

Permission: maintain exact kinetic and first/second moment `G-G` assembly with
function-local scratch reuse.

Non-goals: Gaussian `G-A`/`A-A` raw blocks, unit nuclear, residual/MWG policy,
persistent caches, metadata, artifacts, public API, or solver workflow.

### HP-R3GG-TEST-01 — terminal G-G validation

Lifecycle: completed validation. Permission: test maintenance.

Owner/canonical: `CartesianFinalBasisRealization`;
[terminal G-G products](r3_terminal_gg_product_matrices.md).

Test/evidence:

```text
test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl
docs/src/developer/pqs_manager_running_log.md (Passes 087-087B)
```

Dependencies: `HP-R3GG-FN-01` and accepted Be2/Cr2 ignored parity evidence.

Permission: maintain the focused endpoint, parity, finiteness, and symmetry
coverage. No separate product-framework fixture is authorized.

Non-goals: route, artifact, public workflow, or Cr2 production tests.

### HP-R3UN-FN-01 — terminal unit-nuclear U_GG Gaussian sum

Lifecycle: implemented. Permission: source maintenance.

Owner/canonical: `CartesianFinalBasisRealization`;
[unit-nuclear Gaussian sum](r3_unit_nuclear_ugg_gaussian_sum.md).

Source:

```text
src/cartesian_final_basis_realization/pqs_terminal_one_body.jl
src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl
```

Dependencies: realized terminal blocks and the producer-wide Coulomb expansion.

Permission: maintain exact uncharged by-center `U_GG`, term-first assembly,
and function-local scratch reuse.

Non-goals: physical charge application, neutral `G-A`/`A-A` raw blocks,
residual/MWG policy, persistent caches, metadata, artifacts, or public API.

### HP-R3UN-TEST-01 — terminal unit-nuclear validation

Lifecycle: completed validation. Permission: test maintenance.

Owner/canonical: `CartesianFinalBasisRealization`;
[unit-nuclear Gaussian sum](r3_unit_nuclear_ugg_gaussian_sum.md).

Test/evidence:

```text
test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl
docs/src/developer/pqs_manager_running_log.md (Passes 089-089B)
```

Dependencies: `HP-R3UN-FN-01` and accepted Be2/Cr2 ignored parity evidence.

Permission: maintain the focused endpoint, fallback, parity, finiteness, and
symmetry coverage. No separate Gaussian-sum fixture is authorized.

Non-goals: route, artifact, public workflow, or Cr2 production tests.

### HP-R3BASE-FN-01 — same-construction base K/U reuse

Lifecycle: implemented. Permission: source and approved-caller maintenance.

Owner/canonical: `CartesianFinalBasisRealization` plus producer composition;
[same-construction base reuse](r3_same_construction_base_reuse.md).

Source:

```text
src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl
src/cartesian_base_hamiltonian.jl
```

Additional live consumer: `src/cartesian_protected_ladder_bundle.jl`, under
its separately owned protected-ladder authority.

Dependencies: `HP-R3GG-FN-01`, `HP-R3UN-FN-01`, one shared base construction,
and exact augmented transforms.

Permission: maintain trusted `base_kinetic` and `base_unit_nuclear` handoff,
dimension/center checks, and live exact recomputation fallbacks.

Non-goals: public inputs, persisted trust proofs, caches, provenance or artifact
fields, operator rewrites, residual/MWG policy, or solver behavior.

### HP-R3BASE-TEST-01 — same-construction reuse validation

Lifecycle: completed validation. Permission: test maintenance.

Owner/canonical: producer composition;
[same-construction base reuse](r3_same_construction_base_reuse.md).

Test/evidence:

```text
test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl
docs/src/developer/pqs_manager_running_log.md (Passes 109-110)
```

Dependencies: implemented reuse and fallback paths plus accepted ignored
allocation/parity evidence.

Permission: maintain focused fallback/reuse parity, exact-operator, endpoint,
and artifact-readback coverage.

Non-goals: dedicated cache/provenance tests, Cr2 workflow, or public API tests.

### HP-R3BASE-DRV-WIRE-01 — canonical-driver K/U reuse wiring

Lifecycle: implemented. Permission: caller maintenance.

Owner/canonical: canonical Cartesian driver;
[same-construction base reuse](r3_same_construction_base_reuse.md).

Source:

```text
bin/cartesian_ham_builder.jl
```

Dependencies: `HP-R3BASE-FN-01` and the supported supplemented driver path.

Permission: maintain passing the already-built base kinetic and by-center unit
nuclear matrices to supplemented exact-operator construction.

Non-goals: new driver inputs, hooks, timing labels, stages, artifacts, kernels,
or base-only behavior changes.

### HP-R3BASE-DRV-TEST-01 — canonical-driver reuse validation

Lifecycle: completed validation. Permission: evidence maintenance.

Owner/canonical: canonical Cartesian driver;
[same-construction base reuse](r3_same_construction_base_reuse.md).

Test/evidence:

```text
test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl
docs/src/developer/pqs_manager_running_log.md (Pass 110)
```

No dedicated committed driver test exists for these keyword handoffs.

Dependencies: `HP-R3BASE-DRV-WIRE-01` and supported H2/Be2 supplemented
driver evidence.

Permission: maintain the accepted call-site validation record.

Non-goals: new fixtures, hooks, public inputs, artifacts, or Cr2 workflow.

## Implemented Canonical Cartesian Driver

Canonical contract:
[Cartesian driver usability workflow](cartesian_driver_usability_workflow.md).

### HP-DRV-FILE-01 — canonical driver file

Lifecycle: implemented. Permission: source maintenance.

Canonical contract:
[Cartesian driver usability workflow](cartesian_driver_usability_workflow.md).

Owner/source: `bin/cartesian_ham_builder.jl`.

Permission: maintain the canonical trusted local scientific driver. No other
`bin`, tool, source, test, or committed input fixture is authorized by this ID.

### HP-DRV-FN-01 — compact functional driver workflow

Lifecycle: implemented. Permission: source maintenance.

Canonical contract:
[Cartesian driver usability workflow](cartesian_driver_usability_workflow.md).

Owner/source: `bin/cartesian_ham_builder.jl`, with separately authorized
non-exported producer stages.

Permission: maintain trusted input-file and `key=value` overrides, visible
system/basis/supplement construction, `basisname === nothing` base selection,
coarse physics timing, terminal inventory/due diligence, artifact write, and
optional readback. Exact live inputs and defaults are canonical in the linked
driver contract.

Boundary: the driver does not currently expose `coulomb_accuracy`; omission
uses the producer's compact default. It is not a parser framework, second
public API, route diagnostic, solver, artifact-schema editor, benchmark
harness, or Cr2-specific workflow.

### HP-DRV-NEST-FN-01 — construction-family driver input

Lifecycle: implemented. Permission: source maintenance.

Owner/canonical: `bin/cartesian_ham_builder.jl`;
[Cartesian driver usability workflow](cartesian_driver_usability_workflow.md).

Permission: maintain visible `nesting = :pqs | :wl`, default `:pqs`, in driver
contract/summary/readback facts. It is a construction-family choice, not a
route diagnostic.

### HP-DRV-NEST-WIRE-01 — construction-family route mapping

Lifecycle: implemented. Permission: source maintenance.

Owner/canonical: driver and base normalization;
[Cartesian driver usability workflow](cartesian_driver_usability_workflow.md).

Source: `bin/cartesian_ham_builder.jl` and
`src/cartesian_base_hamiltonian.jl`.

Permission: map `:pqs` to `:pqs_source_box` and `:wl` to
`:white_lindsey_low_order`, preserve public stage/artifact behavior, and reject
unsupported combinations without exposing internal route vocabulary.

### HP-DRV-NEST-TEST-01 — construction-family validation

Lifecycle: completed validation contract. Permission: validation maintenance.

Owner/canonical: driver construction-family input;
[Cartesian driver usability workflow](cartesian_driver_usability_workflow.md).

Evidence: default PQS, supported WL base/supplemented artifact/readback, and
unsupported-combination rejection smokes.

Family-wide non-goals: route algorithms, shell/lowering policy, old WL
materialization, raw-block/RG/MWG/IDA changes, artifacts, public API, solvers,
fixtures, or Cr2-specific workflow.

### HP-DRV-STAGE-FN-01 — visible physics-stage producer surface

Lifecycle: implemented. Permission: source maintenance.

Owner/canonical: staged base/supplement producer;
[Cartesian driver usability workflow](cartesian_driver_usability_workflow.md).

Source: `src/cartesian_base_hamiltonian.jl` plus behavior-preserving operator
factoring in `pqs_source_box_low_order_materialization.jl`,
`pqs_terminal_one_body.jl`, and `pqs_terminal_residual_gto.jl`.

Permission: maintain separate non-exported stages for working basis,
product/moment, unit-nuclear, IDA/MWG interaction, residual augmentation, and
Hamiltonian assembly so the driver can bind and time physical objects. Facades
remain wrappers over the same construction.

Non-goals: public API, route-stage payloads, diagnostics/status clouds,
persistent caches, raw-block switches, kernel instrumentation, solvers, or
artifact changes.

### HP-DRV-STAGE-WIRE-01 — canonical driver staged wiring

Lifecycle: implemented. Permission: source maintenance.

Owner/canonical: `bin/cartesian_ham_builder.jl`;
[Cartesian driver usability workflow](cartesian_driver_usability_workflow.md).

Permission: call the named producer stages directly and print coarse timings
for product/moment, unit-nuclear, interaction, and assembly work. Do not replace
them with an opaque wrapper or expose underscored route stages, stop controls,
providers, allocation probes, or solver controls.

### HP-DRV-STAGE-TEST-01 — staged driver validation

Lifecycle: completed validation contract. Permission: validation maintenance.

Owner/canonical: staged driver workflow;
[Cartesian driver usability workflow](cartesian_driver_usability_workflow.md).

Evidence: bounded base and supplemented driver timing/summary,
artifact/readback, expected-dimension, and contract-printing smokes.

Non-goals: committed fixtures/tests, route diagnostics, solver runs, or
Cr2-specific driver behavior.

### HP-DRV-INV-FN-01 — canonical driver terminal-region inventory

Lifecycle: implemented. Permission: source maintenance.

Owner/canonical: base producer and canonical driver;
[Cartesian driver usability workflow](cartesian_driver_usability_workflow.md).

Source: `src/cartesian_base_hamiltonian.jl` and
`bin/cartesian_ham_builder.jl`, with compact native accessors only where
already owned.

Permission: maintain the bounded human-facing terminal-region inventory with
region/lowering/realization kind, shell index, support/retained counts,
compression, identity/product status, index/physical bounds, slab facts, and
base/supplemented dimensions.

Non-goals: route/source/pair/raw-block dumps, new controls, broad payloads,
artifacts/readers, numerical or shell policy, RG/MWG/IDA, solvers, or Cr2.

### HP-DRV-INV-TEST-01 — terminal-region inventory validation

Lifecycle: completed validation contract. Permission: validation maintenance.

Owner/canonical: terminal-region inventory;
[Cartesian driver usability workflow](cartesian_driver_usability_workflow.md).

Evidence: bounded PQS/WL base and supplemented inventory output,
artifact/readback parity, required geometry/slab columns, and bounded-output
checks. No dedicated committed fixture is owned by this ID.

### HP-DRV-SHELLDD-FN-01 — terminal shellification due-diligence report

Lifecycle: implemented. Permission: source maintenance.

Owner/canonical: base producer report and canonical driver presentation;
[terminal shellification due diligence](terminal_shellification_due_diligence.md).

Source: `src/cartesian_base_hamiltonian.jl` and
`bin/cartesian_ham_builder.jl`, with compact native accessors only where
already owned.

Permission: maintain the bounded in-memory/report table joining terminal
inventory with retained/support facts; system/geometry, axis/center/weight,
dimension, and shell-row diagnostics; actual and expected source shapes;
retained/final ranges; slab metadata; and advisory warning flags.

Contract: producer/driver workflows expose the report, consumers inspect it
before interpreting energies or basis behavior, and warnings remain advisory
unless separate policy makes them fatal. Gausslet/IDA weight summaries are
diagnostic and are not MWG/residual weights or proof of quadrature quality.

Non-goals: new public inputs, artifacts/readers, shell/source/retained policy,
new longitudinal-resolution rules, broad payloads or dumps, RG/MWG/IDA,
Hamiltonian or solver changes, and Cr2 workflow.

### HP-DRV-SHELLDD-TEST-01 — terminal shellification due-diligence validation

Lifecycle: completed validation contract. Permission: validation maintenance.

Owner/canonical: terminal due diligence;
[terminal shellification due diligence](terminal_shellification_due_diligence.md).

Evidence: bounded H2/H2+ producer/driver reports, rectangular-shell warning,
system/axis/weight/dimension/shell-row inspection, bounded-output checks, and
artifact/readback parity. No dedicated committed fixture is owned by this ID.

### HP-PQS-ASPECTSHELL-FN-01 — PQS aspect-aware source modes

Lifecycle: implemented. Permission: source maintenance.

Owner/canonical: PQS terminal low-order route and multilayer shell realization;
[PQS complete-shell aspect-aware source modes](pqs_complete_shell_aspect_source_modes.md).

Source: `src/pqs_source_box_route_driver_helpers.jl`,
`src/pqs_multilayer_shell_region_plan.jl`,
`src/pqs_multilayer_shell_source_plan.jl`, and due-diligence shape reporting in
`src/cartesian_base_hamiltonian.jl`.

Permission: preserve the existing post-shellification angular-band `L`
selection and one authoritative `(q,q,L)` shape through lowering, retention,
realization, and reporting.

### HP-PQS-ASPECTSHELL-TEST-01 — aspect-source validation

Lifecycle: completed validation contract. Permission: validation maintenance.

Owner/canonical: PQS aspect-aware shell source policy;
[PQS complete-shell aspect-aware source modes](pqs_complete_shell_aspect_source_modes.md).

Evidence: bounded H2/H2+ noncubic retained-count, finite/symmetric
artifact/readback, and due-diligence checks recorded in manager Pass 247.

Non-goals: a new longitudinal-resolution policy, public inputs, WL policy,
shell ownership, thin-slab/direct-core changes, artifacts, RG/MWG/IDA,
injection, solvers, or Cr2 claims.
### HP-DRV-TEST-01 — driver workflow validation

Lifecycle: completed validation contract. Permission: validation maintenance.

Owner/canonical: canonical driver;
[Cartesian driver usability workflow](cartesian_driver_usability_workflow.md).

Evidence: bounded base and supplemented H2 contract/summary,
artifact/readback, expected-dimension, and optional Be2 usability smokes.
Inputs remain ignored local files; this ID owns no committed fixture or solver
gate.

## Implemented Canonical Driver Atom Workflow

Canonical contract:
[Cartesian driver atom workflow](cartesian_driver_atom_workflow.md).

### HP-DRV-ATOM-FN-01 — explicit base atom driver workflow

Lifecycle: implemented. Permission: source maintenance.

Canonical contract:
[Cartesian driver atom workflow](cartesian_driver_atom_workflow.md).

Owner/source: `bin/cartesian_ham_builder.jl`.

Permission: maintain `Natom=1`, `basisname === nothing` base-atom selection;
explicit origin, charge, spin-sector, neutral-count, `ns`, `core_spacing`,
`s_factor`, source-span/nesting, and radius-from-padding inputs; and clear
unsupported-input failures. There is no public `mode=:base` input.

### HP-DRV-ATOM-WIRE-01 — driver atom-to-base-facade wiring

Lifecycle: implemented. Permission: source maintenance.

Canonical contract:
[Cartesian driver atom workflow](cartesian_driver_atom_workflow.md).

Owner/source: `bin/cartesian_ham_builder.jl`.

Permission: pass the explicit atom contract through the same named base
producer stages, artifact writer, provenance, and readback as the base facade.
Supplemented atoms remain separately governed by `HP-COMP-SUPPATOM-*`.

Non-goals: package algorithm changes, atom-only Hamiltonian construction,
route reports/status payloads, artifact fields, or unsupported atom broadening.

### HP-DRV-ATOM-TEST-01 — base atom driver validation

Lifecycle: completed validation contract. Permission: validation maintenance.

Owner/canonical: driver atom workflow;
[Cartesian driver atom workflow](cartesian_driver_atom_workflow.md).

Evidence: origin-centered H base artifact/readback and malformed origin,
electron-count, size/spacing, and unsupported-input smokes. Supplemented atom
validation remains owned by `HP-COMP-SUPPATOM-TEST-01`.

### HP-DRV-ATOM-CLEAN-01 — remove hidden atom `d` driver residue

Lifecycle: implemented cleanup. Permission: source preservation.

Canonical contract:
[Cartesian driver atom workflow](cartesian_driver_atom_workflow.md).

Owner/source: `bin/cartesian_ham_builder.jl`.

Permission: preserve the absence of hidden atom `d`; visible atom basis uses
`ns`, `core_spacing`, `radius`, and current optional fields. Do not restore
the compatibility field or use this ID for new inputs, diagnostics, source
algorithms, tools, fixtures, or Cr2 workflow.

## Completed Complete-Core-Shell RHF Retirement

Status: completed and closed by `28e9b2c84`. The entries below are historical
deletion authority and no longer authorize source or test work.

This section records the deletion lane in
`complete_core_shell_rhf_retirement.md`. The old complete-core-shell RHF stack
was stale route-era workflow machinery, not a live Cartesian Hamiltonian
producer path.

### HP-RETIRE-CCS-RHF-FN-01 — remove stale RHF payload stack

Status: completed and closed by `28e9b2c84`.

Permission: none.

Owner/history:
[complete-core-shell RHF retirement](complete_core_shell_rhf_retirement.md).

Approved source files:

```text
src/GaussletBases.jl
src/pqs_multilayer_complete_core_shell_rhf.jl
```

Approved behavior:

- remove the `pqs_multilayer_complete_core_shell_rhf.jl` include from
  `src/GaussletBases.jl`;
- delete `src/pqs_multilayer_complete_core_shell_rhf.jl`;
- remove only docs/index references that describe this RHF stack as active
  current code, if encountered during the deletion pass;
- add no replacements, adapters, compatibility wrappers, status objects,
  payload objects, reports, or tests.

Pre-deletion evidence: focused search found no live `src`, `bin`, `test`, or
`tool` caller outside the file itself and the root include. The file carried
old payload and blocked-status vocabulary such as
`pqs_multilayer_complete_core_shell_rhf_input_contract`,
`pqs_multilayer_complete_core_shell_rhf_scf_payload`, and
`pqs_multilayer_complete_core_shell_rhf_one_step_payload`. The current
CR2-facing producer path is `bin/cartesian_ham_builder.jl`, staged producer
functions, and `CartesianIDAHamiltonian` artifacts.

This ID does not approve canonical driver changes, source changes outside the
approved file/include except minimal stale active-reference cleanup, changes to
`pqs_multilayer_complete_core_shell_h1.jl`,
`pqs_complete_core_shell_final_basis.jl`, or
`pqs_source_box_low_order_materialization.jl`, ordinary/Qiu-White donor-kernel
changes, artifact schema/provenance/reader changes, route/shellification/
terminal-lowering/raw-block/RG/MWG/IDA changes, Hamiltonian assembly changes,
committed tests/fixtures, or Cr2 workflow.

Historical failure rule: if any live `src`, `bin`, `test`, or `tool` caller
depends on the RHF stack, make no source commit and report the exact caller.
Do not preserve the path through an adapter.

### HP-RETIRE-CCS-RHF-TEST-01 — retirement validation

Status: completed validation evidence; closed to new test work.

Permission: none.

Owner/history:
[complete-core-shell RHF retirement](complete_core_shell_rhf_retirement.md).

Approved validation:

- `git diff --check`;
- package load;
- focused `rg` showing no remaining live references to
  `pqs_multilayer_complete_core_shell_rhf`;
- canonical small base artifact/readback smoke;
- canonical small supplemented artifact/readback smoke;
- H2 Residual Gaussian endpoint remains unchanged;
- no Cr2 run.

No committed test, fixture, replacement path, adapter, status/report/payload
object, artifact-schema validation, or Cr2 workflow is approved.

## Completed Route-Driver Materialization Workflow Retirement

Status: completed and closed by `e2e164e9b`; the ladder-runner follow-up was
completed by `77fa2700b`. The entries below are historical deletion authority
and no longer authorize source, test, tool, or docs work.

This section records the retirement/quarantine lane in
`route_driver_materialization_retirement.md`. The old route-driver
materialization/report/save wrapper workflow is not the canonical Cartesian
producer path; the current path is the staged human-facing driver plus
`CartesianIDAHamiltonian` artifacts.

### HP-RETIRE-DRV-MAT-FN-01 — remove old materialization/report/save wrappers

Status: completed and closed by `e2e164e9b`.

Permission: none.

Owner/history:
[route-driver materialization retirement](route_driver_materialization_retirement.md).

Approved source files:

```text
src/pqs_source_box_route_driver_helpers.jl
src/pqs_source_box_low_order_materialization.jl
src/pqs_source_box_route_driver_reporting.jl
src/GaussletBases.jl
```

`src/GaussletBases.jl` is approved only if the reporting include becomes unused
after the wrapper/report/save path is removed.

Approved behavior:

- remove `cartesian_materialization`, `cartesian_print_summary`,
  `cartesian_print_details`, and `cartesian_save`;
- remove `_pqs_source_box_route_driver_materialization`,
  `_pqs_source_box_route_driver_print_materialization`, and
  `_pqs_source_box_route_driver_save` if they are used only by those wrappers;
- remove a related old White-Lindsey atomic pure-gausslet materialization
  helper only if it becomes uncalled;
- do not add replacement wrappers, adapters, status fields, payload objects, or
  tests.

Pre-deletion evidence: focused search found no hits for the audited names in
`bin/cartesian_ham_builder.jl`. The CR2-facing artifact workflow used the
canonical staged producer and `CartesianIDAHamiltonian` artifacts. Remaining
hits were old wrapper definitions, tools/harnesses, stale docs-policy
assertions, and stale compact-doc references.

This ID does not approve canonical driver changes, current staged producer
function changes, artifact schema/provenance/reader/manifest changes, route,
shellification, terminal-lowering, raw-block, Residual Gaussian, MWG, IDA, or
Hamiltonian assembly changes, changes to
`pqs_multilayer_complete_core_shell_h1.jl`,
`pqs_complete_core_shell_final_basis.jl`, broad ordinary/Qiu-White donor-kernel
retirement, replacement wrappers, adapters, status fields, payloads, new tests,
or Cr2 workflow.

Historical failure rule: if any current canonical producer path or public
artifact workflow depends on these wrappers, make no source commit and report
the exact dependency. Do not preserve the wrapper workflow through compatibility
adapters.

### HP-RETIRE-DRV-MAT-TOOL-01 — old wrapper-tool quarantine

Status: completed and closed by `e2e164e9b`.

Permission: none.

Owner/history:
[route-driver materialization retirement](route_driver_materialization_retirement.md).

Approved tool files:

```text
tools/cartesian_driver_harness.jl
tools/cr2_cartesian_ida_stage_probe.jl
tools/cartesian_driver_ladder_lib.jl
tools/h2_pqs_base_hamiltonian_smoke.jl
```

Approved behavior:

- delete or quarantine only old tools that exist to drive the retired wrapper
  workflow;
- do not move route diagnostics, ladder probing, stage stops, or wrapper
  behavior into `bin/cartesian_ham_builder.jl`;
- do not create replacement tool frameworks.

### HP-RETIRE-DRV-MAT-DOC-01 — active docs cleanup

Status: completed and closed by `e2e164e9b`.

Permission: none.

Owner/history:
[route-driver materialization retirement](route_driver_materialization_retirement.md).

Approved docs files:

```text
docs/src/developer/algorithm_implementation_index.md
docs/src/developer/designs/cartesian_hamiltonian_producer/implementation_slices.md
docs/src/developer/designs/cartesian_hamiltonian_producer/registry.md
docs/src/developer/designs/cartesian_hamiltonian_producer/r1_public_base_producer.md
docs/src/developer/pqs_manager_running_log.md
```

Approved behavior:

- stop describing the old wrapper workflow as canonical or active current
  authority;
- keep historical references historical;
- keep the canonical staged driver and current producer/artifact path
  unchanged.

### HP-RETIRE-DRV-MAT-TEST-01 — retirement validation

Status: completed validation evidence; closed to new test work.

Permission: none.

Owner/history:
[route-driver materialization retirement](route_driver_materialization_retirement.md).

Approved test file:

```text
test/docs/cartesian_ham_builder_policy_runtests.jl
```

Approved behavior:

- remove or update only old route-stage wrapper assertions that assume the
  canonical driver should avoid calling these now-retired wrappers.

Approved validation:

- `git diff --check`;
- package load;
- focused `rg` over `src`, `bin`, `test`, and `tools` showing no remaining
  live references to the retired wrapper names;
- canonical small base artifact/readback smoke;
- canonical small supplemented artifact/readback smoke;
- H2 Residual Gaussian endpoint remains unchanged;
- docs-policy test either passes after update or the stale assertion is removed
  under this authority;
- no Cr2 run.

No new committed test or fixture is approved.

### HP-RETIRE-LADDER-RUNNERS-FN-01 — delete dangling ladder runners

Status: completed and closed by `77fa2700b`.

Permission: none.

Owner/history:
[route-driver materialization retirement](route_driver_materialization_retirement.md).

Approved tool files:

```text
tools/run_cartesian_driver_ladder.jl
tools/run_cartesian_line_ladder.jl
```

Approved behavior:

- delete both runner scripts;
- do not add replacements;
- do not modify `bin/cartesian_ham_builder.jl`;
- do not modify `tools/cartesian_driver_ladder_lib.jl` unless a later
  docs-only amendment explicitly approves deleting that quarantined library.

These scripts were only entrypoints into the retired route-driver ladder
workflow after `HP-RETIRE-DRV-MAT-*`. This ID does not approve canonical driver
changes, source changes, test changes except validation scans, artifact/
provenance/reader changes, route/shellification/terminal-lowering/raw-block/
RG/MWG/IDA/Hamiltonian assembly changes, new wrappers, adapters, status fields,
payloads, reports, tools, tests, or Cr2 workflow.

Historical failure rule: if any live source, canonical workflow, or approved
tool still depends on these runner scripts, make no commit and report the
exact dependency. Do not preserve them through an adapter.

### HP-RETIRE-LADDER-RUNNERS-TEST-01 — ladder runner deletion validation

Status: completed validation evidence; closed to new test work.

Permission: none.

Owner/history:
[route-driver materialization retirement](route_driver_materialization_retirement.md).

Approved validation:

- `git diff --check`;
- package load;
- focused `rg` over `src`, `bin`, `test`, and `tools` for
  `run_cartesian_driver_ladder`, `run_cartesian_line_ladder`, and
  `cartesian_driver_ladder_lib`;
- canonical small base artifact/readback smoke;
- no Cr2 run.

No committed test, replacement tool, adapter, or Cr2 workflow is approved.
After this deletion pass, pause the cleanup lane unless a later amendment names
another stale surface.

## Implemented Homonuclear Z-Axis Diatomic Supplemented Workflow

### HP-R3U-ZDI-FN-01 — homonuclear z-axis diatomic supplemented facade

Lifecycle: implemented. Permission: source maintenance.

Owner/canonical: base-Hamiltonian input normalization;
[R3 homonuclear supplemented workflow](r3_homonuclear_diatomic_supplemented_workflow.md).

Source: `src/cartesian_base_hamiltonian.jl`.

Validation: H2 facade gate plus accepted H2/Be2 and optional ignored Cr2
evidence in manager Pass 093; implementation commit `c57e709e7`.

Dependencies: `HP-R3U-FN-01`, explicit system/basis/supplement inputs, and the
shared base producer.

Permission: maintain explicit neutral all-electron homonuclear two-center
z-axis validation and optional trusted supplement `basisfile`.

Non-goals: heteronuclear/general geometry, ECP/charged systems, solver work,
public API, artifact changes, or Cr2-specific branches/defaults/fixtures.

### HP-R3U-ZDI-WIRE-01 — canonical driver supplemented-mode wiring

Lifecycle: implemented. Permission: source maintenance.

Owner/canonical: canonical driver composition;
[R3 homonuclear supplemented workflow](r3_homonuclear_diatomic_supplemented_workflow.md).

Source: `bin/cartesian_ham_builder.jl`.

Validation: accepted H2/Be2 driver artifact/readback evidence in manager Pass
095; original implementation commit `3a4933812`.

Dependencies: `HP-DRV-*`, `HP-DRV-STAGE-*`, and the implemented supplemented
producer stages.

Permission: maintain generic supplemented-mode wiring for explicit homonuclear
z-axis inputs through the same producer construction.

Non-goals: Cr2-specific workflow, committed fixtures, route diagnostics,
artifact changes, public exports, solver behavior, or broad driver growth.

### HP-R3U-ZDI-TEST-01 — homonuclear diatomic validation

Lifecycle: implemented validation contract. Permission: validation
maintenance.

Owner/canonical: homonuclear supplemented workflow;
[R3 homonuclear supplemented workflow](r3_homonuclear_diatomic_supplemented_workflow.md).

Source: no production source permission.

Validation: standalone H2 facade gate, accepted H2/Be2 driver artifact
readback, and optional ignored/user-run Cr2 stress after those gates pass.

Dependencies: `HP-R3U-TEST-01`, `HP-DRV-TEST-01`, and generic homonuclear
input validation.

Permission: maintain bounded H2/Be2 correctness checks and optional ignored
generic stress evidence.

Non-goals: new committed tests or Cr2 fixtures, heteronuclear/general-geometry
gates, solver runs, or artifact schema expansion.

## Measurement Records

These entries record bounded measurement authority and its disposition. A
completed record grants no continuing work unless its own lifecycle says
otherwise.

### HP-COMP-ANGBOX-AUDIT-01 — angular-balanced shellification geometry audit

Lifecycle: completed historical measurement. Permission: none.

Owner/canonical: common shell geometry;
[common terminal shell decomposition](common_terminal_shell_decomposition.md).

Outcome: the ignored angular-balance inventory established the native
angular-z-extension policy later implemented under `HP-COMP-ANGBOX-*`.
It grants no source, test, artifact, driver, solver, or Cr2 authority.

### HP-R3REM-AUDIT-01 — remaining exact-operator allocation audit

Lifecycle: completed historical measurement. Permission: none.

Owner/canonical: historical exact-operator attribution;
[remaining allocation audit](r3_remaining_exact_operator_allocation_audit.md).

Evidence:

```text
docs/src/developer/pqs_manager_running_log.md (Passes 088-089B)
```

Outcome: attributed the largest in-wrapper remainder to terminal unit-nuclear
`U_GG`; the separately authorized implementation is now canonical in
[unit-nuclear Gaussian sums](r3_unit_nuclear_ugg_gaussian_sum.md).

Non-goals: this completed audit grants no source, test, tool, driver, route,
kernel, artifact, public API, solver, or Cr2 workflow authority.

## Rejected Or Deferred

### HP-RES-01 — terminal basis build result — rejected

Lifecycle: rejected. Permission: none.

Canonical boundary:
[terminal basis and base assembly](terminal_basis_and_base_assembly.md).

Do not introduce a persistent terminal-basis result wrapper. The realizer
returns `CartesianTerminalBasisRealization` on success.

### HP-CHANGE-01 — return shell overlap from existing shell plan — rejected/deferred

Lifecycle: rejected as standalone authority. Permission: none.

Canonical boundary:
[terminal basis and base assembly](terminal_basis_and_base_assembly.md).

This can be a helper detail under HP-FN-00, but it is not standalone authority.

### HP-OBJ-03 — generic build-result wrapper — rejected

Lifecycle: rejected. Permission: none.

Canonical boundary:
[terminal basis and base assembly](terminal_basis_and_base_assembly.md).

Do not introduce `CartesianHamiltonianBuildResult`, another payload, or a broad
status wrapper around `CartesianIDAHamiltonian`.

### HP-TEST-01 — new committed terminal smoke — rejected

Lifecycle: rejected. Permission: none.

Canonical boundary:
[terminal basis and base assembly](terminal_basis_and_base_assembly.md).

No new committed terminal smoke/probe is approved. Use existing smokes or
ignored `tmp/work` validation unless a later design explicitly approves a test.

### R3 Be2/Cr2 readiness guardrail — deferred

R3-A/B/C are implemented for the narrow H2 residual-GTO/MWG endpoint and
compact supplemented artifact. The following items are closed for that narrow
path and should not be listed as future blockers by default:

- same-construction internal path for the accepted H2 R3 construction;
- deterministic rank-deficient handling in the legacy global-selection
  implementation, now superseded by the approved owner-local selection source
  correction;
- one-shot parent-by-supplement analytic exact-block organization for R3-A
  mixed/self blocks, avoiding repeated CPB-per-terminal-block construction on
  the Be2 proxy;
- independent weight-aware final-basis `V_GM` validation for R3-B;
- compact `supplement_provenance/` artifact group for R3-C.

Do not present Cr2 or broader residual-GTO/MWG supplement support as approved
until a later docs-only amendment chooses and closes the next lane. The
non-exported H2/Be2 R3 usability facade is approved separately by
`HP-R3U-FILE-01`, `HP-R3U-FN-01`, `HP-R3U-WIRE-01`, and `HP-R3U-TEST-01`.
Remaining deferred lanes are:

- implementation and validation of the approved owner-local source correction,
  including the updated H2 scalar and no full Cr2 Hamiltonian;
- public export, public examples, or driver workflow beyond the non-exported
  R3 usability facade;
- a Cr2-readiness lane for measurement-only candidate/rank/memory forecasting,
  with no full Cr2 Hamiltonian yet;
- a basis/supplement-realism lane for validated supplement choices, basis
  labels, and filtering policy beyond the first H2 fixture;
- bounded or streamed residual MWG term storage if higher residual rank makes
  the current dense residual term storage costly;
- allocation-free or bounded-allocation validation reductions for large dense
  matrices.

### Composition family index

The bounded atom/diatomic, PQS/WL, base/supplemented matrix is implemented.
Its individually addressable source/test entries and evidence appear above;
the live numerical contract is
[Nesting/supplement composition](nesting_supplement_composition_plan.md).
WL compact retained-basis and boundary-count policies remain separate
registered families. General geometry, public export, solver/ECP, and
Cr2-specific workflow remain outside this composition authority.
