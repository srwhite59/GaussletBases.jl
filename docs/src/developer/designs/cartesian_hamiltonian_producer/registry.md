# Registry

Only entries marked approved/implemented authorize work on the exact surface
they describe. Measurement-only entries do not authorize production source
edits. Candidate or rejected entries do not authorize implementation.

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
belong in the linked canonical subsystem document. A heading
that currently lists a source/test pair is one shared lane record and applies
to both IDs; a later machine-registry pass must normalize the IDs without
dropping either one.

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

## Approved For White-Lindsey Terminal Basis Implementation

This section approves only the narrow terminal-basis seam recorded in
`white_lindsey_terminal_basis_realization.md`.

### HP-WLTERM-FILE-01 — optional WL terminal realization file

Approved source files:

```text
src/cartesian_final_basis_realization/white_lindsey_terminal_basis_realization.jl
src/cartesian_final_basis_realization/CartesianFinalBasisRealization.jl
```

This ID is optional. It approves creating a small WL-specific terminal
realization sibling and adding its include to `CartesianFinalBasisRealization`
only if extending `pqs_terminal_basis_realization.jl` directly would obscure
the distinct WL boundary-stratum construction. No public export, root include,
new module, new basis object, artifact, report, or status/payload object is
approved.

### HP-WLTERM-FN-01 — WL low-order terminal basis realization

Approved source files:

```text
src/cartesian_final_basis_realization/pqs_terminal_basis_realization.jl
src/cartesian_final_basis_realization/white_lindsey_terminal_basis_realization.jl
```

Approved target: realize terminal final-basis blocks for the existing
`:white_lindsey_low_order` route family and return the existing
`CartesianTerminalBasisRealization`.

Allowed behavior:

- support direct identity terminal blocks on their authoritative owned rows;
- realize WL boundary-stratum/product terminal blocks from existing native
  terminal support, retained-rule, and transform records;
- use only support-local coefficients on `support.support_indices` /
  `support.support_states`;
- preserve deterministic terminal support, lowering, retained-record, and
  transform-contract order;
- validate disjoint owned supports and block-local identity overlaps under the
  same structural terminal-basis policy as PQS.

This ID does not approve old WL H1/H1+J materialization adaptation, new
Hamiltonian objects, new route-stage objects, route reports, status/result
payloads, diagnostics, artifact/schema changes, public API/export changes,
raw-block changes, Residual Gaussian/MWG/IDA changes, terminal
shellification-policy changes, retained-selection-policy changes, route
skeleton construction changes, or source files outside the approved surfaces.

Failure rule: if WL boundary-stratum final basis cannot be materialized from
existing terminal lowering, retained-unit, and transform records without
broader route redesign, make no source commit and report the exact missing
native fact.

### HP-WLTERM-WIRE-01 — WL route helper terminal-basis wiring

Approved source file:

```text
src/pqs_source_box_route_driver_helpers.jl
```

Approved behavior:

- remove or narrow the PQS-only terminal-basis guard so
  `route_family = :white_lindsey_low_order` can call the approved WL terminal
  realizer when native terminal records are available;
- keep `route_family = :pqs_source_box` behavior unchanged;
- keep route skeleton construction semantics, route recipe semantics,
  shellification behavior, terminal lowering order, retained-rule order,
  public driver contract, and artifact schema unchanged;
- reject or return a clear unsupported route error when the WL route lacks the
  native facts required for terminal-basis realization.

This ID does not approve adapting old WL materialization, broad route redesign,
supplemented WL behavior, driver input changes, diagnostics/status/report
expansion, raw-block switches, stop-after controls, or source files outside the
approved route helper file.

### HP-WLTERM-TEST-01 — WL terminal-basis validation

Approved validation:

- `git diff --check`;
- package load;
- current default `nesting = :pqs` atom or H2 base artifact/readback remains
  unchanged;
- `nesting = :wl` base atom artifact/readback;
- `nesting = :wl` base H2 artifact/readback if the existing WL diatomic route
  has sufficient native terminal records;
- H2 residual-GTO/MWG PQS endpoint remains unchanged if terminal realization
  code is touched;
- clear unsupported-input/blocker report if WL H2 cannot be realized from
  current native records.

No Cr2 run, supplemented WL run, committed fixture, committed test file,
solver/RHF/ECP/EGOI workflow, artifact schema validation, or broad WL workflow
validation is approved.

## Approved First Composition Lane: WL Z-Axis Diatomic Base

This section promotes the first
`nesting_supplement_composition_plan.md` placeholder. It approves only the
`geometry = z-axis diatomic`, `nesting = :wl`, `supplement = off` base path.

### HP-COMP-WLDIAT-FN-01 — WL z-axis diatomic base terminal records

Approved source files:

```text
src/pqs_source_box_diatomic_complete_core_shell.jl
src/cartesian_terminal_shellification_geometry.jl
src/cartesian_terminal_lowering/selection.jl
src/cartesian_terminal_lowering/region_contracts.jl
src/pqs_source_box_route_driver_helpers.jl
src/cartesian_final_basis_realization/pqs_terminal_basis_realization.jl
src/cartesian_final_basis_realization/white_lindsey_terminal_basis_realization.jl
src/cartesian_final_basis_realization/CartesianFinalBasisRealization.jl
src/cartesian_base_hamiltonian.jl
```

Approved goal:

```text
Natom = 2
nesting = :wl
basisname = nothing
```

must produce an existing `CartesianIDAHamiltonian{Float64}` artifact/readback
through the same `CartesianTerminalBasisRealization` and staged base
Hamiltonian path used by the PQS producer.

Allowed behavior:

- produce native White-Lindsey z-axis diatomic terminal support,
  shellification, lowering, retained-rule, and transform records needed by the
  existing WL terminal realizer;
- preserve the existing `:white_lindsey_low_order` construction-family route
  and deterministic support/lowering/retained/transform order;
- realize WL boundary-stratum/product terminal blocks as owned-support
  terminal blocks in the existing `CartesianTerminalBasisRealization`;
- reuse the existing staged product, unit-nuclear, IDA, Hamiltonian
  construction, writer, and reader path;
- in `src/cartesian_base_hamiltonian.jl`, add only narrow staged/facade wiring
  required by WL z-axis diatomic base construction and the truthful route
  provenance value `:z_axis_diatomic_wl_base`.

The `:z_axis_diatomic_wl_base` route value is approved as a value under the
existing `producer_provenance/route` and `recipe_provenance/route` keys. It is
not an artifact schema change.

Forbidden:

- driver public input changes or driver special cases;
- old WL H1/H1+J materialization revival or adaptation;
- artifact schema changes, matrix-key changes, reader behavior changes,
  manifest shape changes, public API/export changes, or new Hamiltonian
  wrapper/result objects;
- Residual Gaussian, MWG/IDA, supplement, ECP, solver/RHF, or Cr2 workflow
  work;
- route diagnostics, stop-after controls, report/status/payload fields,
  raw-block switches, route-stage labels, or broad route skeleton redesign;
- committed tests, committed fixtures, committed driver input files, or source
  files outside the approved list.

Failure rule: if WL z-axis diatomic base construction requires adapting the old
WL H1/H1+J materialization path, changing artifact schema/reader behavior,
adding driver special cases, or creating a parallel Hamiltonian builder, make
no source commit and report the blocker.

Line budget: at most `250` added `src` lines, with deletion or simplification
of obsolete blocker-only WL diatomic guards expected where practical. Stop for
a new amendment if the pass needs broader route skeleton redesign, source
files outside the approved list, or persistent payload/cache objects.

### HP-COMP-WLDIAT-TEST-01 — WL z-axis diatomic base validation

Approved validation:

- `git diff --check`;
- package load;
- current default `nesting = :pqs` H2 base artifact/readback remains
  unchanged;
- existing `nesting = :wl` one-center atom artifact/readback remains
  unchanged;
- `nesting = :wl` z-axis H2 base artifact/readback succeeds through the staged
  base Hamiltonian path;
- direct provenance inspection confirms `nesting = :wl` and
  `route = :z_axis_diatomic_wl_base` for the WL H2 artifact;
- H2 residual-GTO/MWG PQS endpoint remains unchanged if terminal realization
  code is touched;
- no Cr2 run.

No supplemented WL run, committed test file, committed fixture, driver
contract test, solver/RHF/ECP/EGOI validation, route-diagnostic validation, or
Cr2 fixture is approved.

## Approved Correction Lane: WL Z-Axis Diatomic Compact Retained Basis

This section records the follow-up design decision after the WL diatomic
terminal-record endpoint exposed the remaining placeholder-like retained-basis
shape. The current `nesting = :wl`, `Natom = 2` path can be mechanically
realized, but it still follows:

```text
elongated shared complete shell
-> boundary CPB strata
-> one retained identity unit per stratum
-> identity terminal blocks
```

That is not the intended compact White-Lindsey retained basis and should not be
used as the production PQS/WL comparison story. The observed audit evidence was
that bounded WL diatomic `ns = 4/5` examples built an elongated shell with
support-size scale `(5,5,9) - (3,3,7) = 162`, rather than the cubic shell-size
scales `4^3 - 2^3 = 56` and `5^3 - 3^3 = 98`; WL retained-unit lowering then
split that shell into 26 boundary-stratum units that the terminal realizer
appended as full-support identity blocks. For `ns = 6`, contact-core geometry
can consume the bounded support and collapse the terminal basis to one direct
identity block.

### HP-WLDIAT-COMPACT-FN-01 — WL diatomic compact retained basis

Approved source files:

```text
src/cartesian_shellification/terminal_geometry.jl
src/cartesian_terminal_lowering/region_contracts.jl
src/cartesian_retained_units/lower_contract_units.jl
src/cartesian_retained_unit_transform_contracts/unit_contracts.jl
src/cartesian_final_basis_realization/white_lindsey_terminal_basis_realization.jl
src/pqs_source_box_route_driver_helpers.jl
```

`src/pqs_source_box_route_driver_helpers.jl` is approved only for narrow route
wiring if the compact WL retained-unit facts must be passed to the existing WL
terminal-basis seam.

Approved behavior:

- preserve the WL unit-based implementation model: faces, edges, corners, and
  small boundary units after shellification;
- do not force a persistent shell object after retained-unit splitting;
- make each WL unit carry or realize the intended compact retained basis from
  products of one-dimensional contractions on the authoritative owned unit
  support;
- treat identity realization as valid only for true direct/core identity
  units, not for WL boundary-stratum retained units;
- use deleted WL coefficient helpers only as historical donor/reference
  material for the compact CPB-local product-of-1D coefficient primitive;
- preserve deterministic geometry, lowering, retained-unit, transform-contract,
  and terminal-block ordering;
- keep the same public `ns` as the fair starting input for PQS/WL comparison,
  while allowing WL-specific geometry/contact cases and not promising equal
  final dimensions.

Forbidden:

- driver changes;
- artifact schema, provenance, matrix-key, reader, or manifest changes;
- PQS behavior changes;
- Hamiltonian assembly changes;
- raw-block, Residual Gaussian, MWG/IDA, Qiu-White, supplement, solver/ECP,
  or Cr2 workflow changes;
- old WL route-global stack, reports, adapters, or H1/H1+J materialization
  revival or adaptation;
- broad route diagnostics, report/status/payload fields, raw-block switches,
  retained-rule dumps, or route-stage labels;
- fake compactness by dropping support rows, relabeling full-support identity
  units, or changing the driver comparison;
- committed tests, committed fixtures, or committed driver input files.

Failure rule: if compact WL retained units require construction-native facts
that are not currently available, make no source commit and report the exact
missing fact. Do not fake compactness by deleting rows, changing public input
semantics, or rerouting through old WL materialization. If an essential
primitive exists only in deleted WL files, restore or re-express only that
primitive behind the current terminal-basis boundary; do not restore the old
route-global framework around it.

Line budget: at most `250` added `src` lines unless a later source blurb
narrows or revises the budget after auditing the exact live callers.

### HP-WLDIAT-COMPACT-TEST-01 — WL compact-basis validation

Approved validation:

- `git diff --check`;
- package load;
- small H2 or Be2 WL base artifact/readback;
- small H2 or Be2 WL supplemented artifact/readback if the compact base path
  works through the existing supplemented boundary;
- PQS base and supplemented smokes remain unchanged;
- WL retained dimension is compared against expected shell-size scale for
  bounded `ns = 4/5` examples;
- finite/symmetric `K` and `V` checks;
- no Cr2 run.

No committed test file, committed fixture, driver contract test,
solver/RHF/ECP/EGOI validation, route-diagnostic validation, or Cr2 fixture is
approved.

## Approved Correction Lane: WL Boundary-Stratum Retained-Count Parity

This section records the follow-up policy correction after the compact WL
diatomic terminal-basis source pass. The full-support identity bug is fixed,
but `nesting = :wl`, `ns = 4` still follows the inherited symmetric-odd donor
rule and produces 26 boundary columns rather than the nominal shell count
`4^3 - 2^3 = 56`.

That remaining behavior is not an acceptable WL policy. The odd-side rule is a
direct core-block centering requirement: a nucleus-centered direct core should
have odd side length so the nucleus is centered. Boundary shells and boundary
strata outside that core do not require odd side counts and must retain the
requested shell contraction count.

### HP-WLDIAT-PARITY-FN-01 — WL boundary retained-count parity

Approved source file:

```text
src/cartesian_final_basis_realization/white_lindsey_terminal_basis_realization.jl
```

Approved behavior:

- preserve odd-side enforcement for true direct nucleus-centered core blocks;
- for WL boundary shell strata, use the requested boundary retained count
  without symmetric-odd coercion;
- for public `nesting = :wl`, `ns = 4`, route-local `q = 2` must retain the
  shell count `4^3 - 2^3 = 56`;
- for public `nesting = :wl`, `ns = 5`, retain `5^3 - 3^3 = 98`;
- keep the compact WL product-of-1D coefficient construction and existing
  terminal-basis boundary.

Forbidden:

- driver changes;
- public `ns` normalization or route-local `q` rule changes;
- route skeleton, shellification, terminal lowering, retained-unit metadata
  shape, or contract-plan changes;
- direct/core identity behavior changes;
- artifact schema/provenance, matrix-key, reader, or manifest changes;
- PQS behavior changes;
- Hamiltonian assembly changes;
- raw-block, Residual Gaussian, MWG/IDA, Qiu-White, supplement, solver/ECP,
  or Cr2 workflow changes;
- old WL route-global stack, reports, adapters, or H1/H1+J materialization
  revival or adaptation;
- broad route diagnostics, report/status/payload fields, raw-block switches,
  retained-rule dumps, or route-stage labels;
- committed tests, committed fixtures, or committed driver input files.

Failure rule: if fixing boundary parity requires changing source files outside
`src/cartesian_final_basis_realization/white_lindsey_terminal_basis_realization.jl`,
public `ns` semantics, route/shellification/terminal-lowering contracts,
artifact schema, or old WL materialization, make no source commit and report
the exact blocker.

Line budget: target under `30` added `src` lines, with no new persistent shape.

### HP-WLDIAT-PARITY-TEST-01 — WL parity validation

Approved validation:

- `git diff --check`;
- package load;
- WL H2 or Be2 z-axis diatomic `ns = 4` retained boundary count / dimension
  demonstrates 56 boundary columns rather than 26;
- WL H2 or Be2 z-axis diatomic `ns = 5` retained boundary count demonstrates
  98 boundary columns;
- small WL base artifact/readback smoke;
- small WL supplemented artifact/readback smoke if bounded by the existing
  supplemented boundary;
- PQS H2 residual-GTO/MWG endpoint remains unchanged;
- finite/symmetric `K` and `V` checks for the WL smoke;
- no Cr2 run.

No committed test file, committed fixture, driver contract test,
solver/RHF/ECP/EGOI validation, route-diagnostic validation, or Cr2 fixture is
approved.

## Composition Input: Public ns Direct-Core Side

Canonical contract:
[Public ns direct-core side parity](public_ns_core_side_parity.md).

### HP-COMP-NSCORE-FN-01 — public ns direct-core side

Lifecycle: implemented. Permission: source maintenance.

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

Owner/source:
`src/cartesian_shellification/terminal_geometry.jl`, with narrow caller
plumbing in `src/pqs_source_box_route_driver_helpers.jl`.

Permission: maintain route-family-free direct-core/shell regions, ordering,
coverage, and owned support before PQS/WL retained construction diverges.

### HP-COMP-SHELLGEOM-TEST-01 — common shell validation

Status: completed validation evidence.

Evidence: bounded same-input PQS/WL shell/support comparisons and downstream
artifact/endpoint smokes recorded in the manager log. No committed fixture is
owned by this ID.

### HP-COMP-SHELLGEOM-DIAT-FN-01 — diatomic common shellifier entry

Status: implemented.

Owner/source:
`src/cartesian_shellification/terminal_geometry.jl` and narrow
`src/pqs_source_box_route_driver_helpers.jl` caller plumbing.

Permission: feed the same public `ns`, direct-core side, centers, bond axis,
and parent facts into common z-axis diatomic shellification before family
lowering. Central-gap/contact redesign is not approved.

### HP-COMP-SHELLGEOM-DIAT-TEST-01 — diatomic shellifier-entry validation

Status: completed validation evidence.

Evidence: bounded PQS/WL same-function/same-argument and shell-region parity;
no Cr2 committed gate.

### HP-COMP-OUTERMM-FN-01 — outer-mismatch-only correction

Status: superseded; no permission.

Superseded by `HP-COMP-THINSLAB-FN-01`. Do not restore a separate
outer-mismatch path.

### HP-COMP-OUTERMM-TEST-01 — outer-mismatch-only validation

Status: superseded; no permission.

Historical evidence is subsumed by common thin-slab validation.

### HP-COMP-ANGBOX-FN-01 — angular-balanced diatomic shellification

Status: implemented.

Owner/source:
`src/cartesian_shellification/terminal_geometry.jl`.

Permission: emit native ordered `:angular_z_extension_slab` stacks so the
ordinary shell body plus planned axial extensions realizes the physical
outer-nucleus angular target. It does not change real-shell retained policy or
central-gap/contact ownership.

### HP-COMP-ANGBOX-TEST-01 — angular shellification validation

Status: completed validation evidence.

Evidence: ignored H2/Be2/Cr2-style geometry inventories and accepted
shellification Passes 179/186; no production Cr2 claim.

### HP-COMP-FACEPROD-FN-01 — neutral terminal face-product helper

Status: implemented.

Owner/source:

- `src/cartesian_final_basis_realization/terminal_face_product_blocks.jl`;
- module include in
  `src/cartesian_final_basis_realization/CartesianFinalBasisRealization.jl`;
- consumers in PQS and White-Lindsey terminal realization.

Permission: one route-neutral face/face-stack coefficient assembly over fixed
normal-axis indices. It is not a new terminal-basis policy.

### HP-COMP-FACEPROD-TEST-01 — face-product validation

Status: completed validation evidence.

Evidence: exact White-Lindsey facet coefficient parity and compact slab reuse;
no committed fixture.

### HP-COMP-THINSLAB-FN-01 — common compact thin-slab lowering

Status: implemented.

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

Evidence: PQS/WL H2/Be2 artifact/readback, compact retained counts, exact WL
facet parity, and H2 residual-GTO endpoint replay; no committed Cr2 gate.

### HP-COMP-THINSLAB-META-FN-01 — thin-slab inventory metadata

Status: implemented.

Owner/source:
`src/cartesian_terminal_shellification_geometry.jl`.

Permission: describe native compact slab kinds consistently in internal
inventory/scaffold summaries. It does not materialize coefficients or create
artifact/report payloads.

### HP-COMP-THINSLAB-META-TEST-01 — thin-slab inventory validation

Status: completed validation evidence.

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

Source owners:
`src/cartesian_nested_faces.jl` and narrow existing PQS source-axis plumbing.
No new production file was created under `CartesianRawProductSources`.

### HP-MCOMX-OBJ-01 — mapped source specification

Status: implemented internal construction specification.

Permission: fixed protected-P2, mapped-Chebyshev, lambda/no-sqrt-J, and
physical-localization facts. No public export or general tuning object.

### HP-MCOMX-FN-01 — mapped source-span construction

Status: implemented.

Owner/source:
`src/cartesian_nested_faces.jl`.

Permission: construct normalized-local-coordinate mapped enrichment before the
existing physical-coordinate COMX cleanup. Ordinary behavior remains default.

### HP-MCOMX-WIRE-01 — PQS axis-transform wiring

Status: implemented.

Owner/source:
`src/cartesian_pair_block_materialization/pqs_source_axis_transforms.jl`.

Permission: pass the internal source-span choice into the existing doside seam
and return ordinary carried `AxisSourceTransformFact` objects.

### HP-MCOMX-TEST-01 — mapped source validation

Status: completed validation evidence.

Evidence: `ns=5/6/7` rank/protected-span/coordinate/orthogonality checks and
bounded H/He/H2 comparisons. He `ns=5` blocks default promotion.

### HP-MCOMX-TERM-FN-01 — terminal shell-seed consumption

Status: implemented.

Owner/source:

- `src/cartesian_final_basis_realization/pqs_terminal_basis_realization.jl`;
- narrow module import/include support in
  `src/cartesian_final_basis_realization/CartesianFinalBasisRealization.jl`.

Permission: validate and consume materialized carried axis facts as the PQS
shell seed while preserving ordinary fallback, boundary selection, support
restriction, Lowdin, and canonicalization.

### HP-MCOMX-TERM-TEST-01 — terminal seam validation

Status: completed validation evidence.

Evidence: carried-fact coefficient parity/difference checks plus ordinary and
supplemented bounded endpoints; no committed Cr2 gate.

### HP-MCOMX-DRV-FN-01 — canonical source-span selector

Status: implemented.

Owner/source:

- `bin/cartesian_ham_builder.jl`;
- `src/cartesian_base_hamiltonian.jl`;
- narrow `src/pqs_source_box_route_driver_helpers.jl` propagation.

Permission: expose `source_span = :ordinary | :mapped_comx`, default ordinary,
and reject mapped-COMX with White-Lindsey. This is not a diagnostic route
switch.

### HP-MCOMX-DRV-TEST-01 — driver source-span validation

Status: completed validation evidence.

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

Evidence and test surfaces:

- `docs/src/developer/pqs_manager_running_log.md`, Pass 139;
- `test/driver_public/cartesian_base_hamiltonian_runtests.jl`, current
  downstream base-facade regression.

The accepted gate covered H2 and bounded Be2 PQS/WL artifact/readback plus
invalid geometry, charge, and electron-count rejection. No dedicated
ID-specific committed test file was added.

### HP-COMP-SUPPWL-FN-01 — supplemented WL z-axis diatomics

Lifecycle: implemented. Permission: source maintenance.

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

Evidence and test surfaces:

- `docs/src/developer/pqs_manager_running_log.md`, Pass 141;
- `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl`, current
  downstream supplemented-path regression.

The accepted gate covered H2 PQS/WL and bounded Be2 WL artifact/readback,
finite/symmetric operators, and invalid-system rejection. No dedicated
ID-specific committed test file was added.

### HP-COMP-SUPPATOM-FN-01 — supplemented one-center atoms

Lifecycle: implemented. Permission: source maintenance.

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

Evidence and test surfaces:

- `docs/src/developer/pqs_manager_running_log.md`, Pass 146;
- `test/driver_public/cartesian_base_hamiltonian_runtests.jl`, current
  downstream one-center base regression.

The accepted gate covered atom base/supplemented PQS/WL construction, bounded
Be, radius sensitivity, and H2/Be2 regression smokes. No dedicated
ID-specific committed test file was added.

### HP-COMP-NS-FN-01 — public ns and derived q

Lifecycle: implemented. Permission: source maintenance.

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

Evidence and test surfaces:

- `docs/src/developer/pqs_manager_running_log.md`, Pass 145;
- `test/driver_public/cartesian_base_hamiltonian_runtests.jl`, current
  downstream public-input/provenance regression.

The accepted gate covered atom/diatomic PQS/WL base construction, bounded
supplemented smokes, legacy-`q` compatibility, inconsistent input rejection,
and provenance. No dedicated ID-specific committed test file was added.

### HP-COMP-WLNS-FN-01 — WL diatomic ns guard

Lifecycle: implemented. Permission: source maintenance.

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

### HP-PQS-ATOMREF-PACKET-FN-01 / HP-PQS-ATOMREF-PACKET-TEST-01 — reusable atomic HF reference packets

Status: implemented internal facility.

Owner: `CartesianReferenceDensity` atomic packet subsystem.

Canonical contract:

- [Atomic HF reference packets](atomic_hf_reference_packets.md)

Approved source surface:

- `src/cartesian_reference_density/CartesianReferenceDensity.jl`;
- `src/cartesian_reference_density/atomic_hf_reference_packets.jl`;
- `src/GaussletBases.jl` only for include/qualified access wiring;
- `src/cartesian_reference_density/screened_hartree_correction.jl` only for
  packet-convergence rejection at consumption;
- optional narrow use of existing exact Hartree helpers in
  `src/cartesian_gaussian_raw_blocks/mixed_hartree_blocks.jl` for validation
  and bounded packet consumption without changing their contract;
- `data/legacy/BasisSets` header only and
  `docs/src/developer/legacy_basissets_provenance.md` for the bounded vendored-data
  provenance correction; the scientific body is not editable under this ID.

Approved test surface:

- `test/nested/cartesian_atomic_hf_reference_packet_runtests.jl`;
- `test/nested/cartesian_screened_hartree_correction_runtests.jl` only for
  unconverged packet-consumer rejection and protected owner-local embedding
  equivalence/failure coverage;
- `test/misc/runtests.jl` only for the vendored-basis identity/parser
  regression.

Permission summary: build, validate, write/read, and boundedly consume the
implemented Be-core and Ne-all-electron one-center packet references under the
canonical determinant/density/potential contract. Stored packet overlap
fingerprints remain exact self-integrity checks. A translated/reconstructed
owner-local supplement block must match all identity/order/owner/placement
facts exactly, but its overlap is accepted by the unchanged numerical
`norm(..., Inf) <= 1e-10` gate rather than raw-byte fingerprint equality. One
nested internal mapping summary may carry the three fingerprints, mapped exact-
match boolean, two overlap differences, and tolerance.

Packet fitting uses only the ordinary determinant -> density-fit -> radial-
potential-fit pipeline. Density fits own `E0`; the current 33-term potential
fit is an approximate `J0` evaluator built from the compact 45-term scaffold.
Its radial, tail, matrix, and `Tr(P0*J0_fit)-E0_fit` consistency errors are
reported. Existing packets containing retired moment-polish provenance are
rejected and require regeneration.

Amendment source limit: the embedding-equivalence follow-on may edit only
`src/cartesian_reference_density/atomic_hf_reference_packets.jl`, plus the
existing private additive-reference caller if the nested diagnostic must be
forwarded. This does not reopen the other packet or correction surfaces.

Retirement cleanup under these packet IDs must delete
`_POTENTIAL_MOMENT_DISTANCES`, `_determinant_potential_moments`,
`_polish_atomic_reference_potential`, packet fields, writer/readback support,
and focused tests; make packet construction consume the ordinary radial fit
directly; and reject retired provenance. Matching energy-consistency behavior
belongs to the active screened-Hartree and protected-additive IDs. No
compatibility adapter or replacement fit policy is approved; the source/test
change should be materially line-negative.

Non-goals: production corrected Hamiltonians, artifact integration beyond the
packet, public defaults, solver workflow, exchange, EGOI, row-gauge rho0/P0,
Cr/Cr2 claims, inferred occupancy, or fitted terms as protected orbitals.

### HP-PQS-ATOMREF-POTMOM-FN-01 / HP-PQS-ATOMREF-POTMOM-TEST-01 - determinant-moment fitted-potential polish

Lifecycle: retired. Permission: none.

This same-day false start was approved in commit `9739c22a6` and implemented
as the moment-polish portion of manager Pass 353. It adjusted fitted-potential
coefficients against determinant moments on a fixed separation grid to force
one padded Be2 energy-consistency value below `1e-8 Ha`. That is not a generic
atomic or molecular fitting principle and must not be retained through an
adapter.

Historical evidence remains in manager Passes 351 and 353. Durable packet,
screened-Hartree, and additive-reference behavior is now owned by their active
IDs and canonical contracts; the ordinary radial fit supersedes the retired
behavior. The retired IDs do not authorize source work, packet compatibility,
molecule-trained fitting, or consumption of polished packets.

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

### HP-PQS-SCREEN-HARTREE-CORR-FN-01 / HP-PQS-SCREEN-HARTREE-CORR-TEST-01 - internal screened-Hartree correction assembly

Status: implemented internal facility.

Owner and canonical contract:

- `CartesianReferenceDensity`;
- [Screened Hartree correction assembly](screened_hartree_correction_assembly.md).

Implemented source surfaces:

- `src/cartesian_reference_density/CartesianReferenceDensity.jl`;
- `src/cartesian_reference_density/screened_hartree_correction.jl`;
- `src/cartesian_reference_density/atomic_hf_reference_packets.jl` for narrow
  packet validation and field evaluation;
- `src/GaussletBases.jl` for module wiring.

Implemented test surface:

- `test/nested/cartesian_screened_hartree_correction_runtests.jl`.

Dependencies:

- [Screened Hartree residual-density formalism](screened_hartree_residual_density.md);
- [Atomic HF reference packets](atomic_hf_reference_packets.md);
- [Protected additive atomic reference correction](protected_additive_reference_correction.md)
  for molecular protected-localized construction.

Permission summary: consume represented, converged reference determinants and
same-basis `V_IDA`, `J0_G`, and `E0_G`; validate representation, symmetry,
and direct field algebra; and return

```text
Delta_J0 = J0_G - Diagonal(V_IDA * q0)
C = 0.5 * q0' * V_IDA * q0 - 0.5 * E0_G
```

as an in-memory `ScreenedHartreeCorrection`. These terms belong to direct
electron-electron/Hartree accounting, not physical kinetic/nuclear `H1`.

Exact/density-fit oracle fields retain strict energy identities. For an
ordinary fitted-potential field, report
`Tr(P0*J0_fit)-E0_fit` rather than rejecting the result solely at `1e-8 Ha`.
The active correction IDs approve the narrow source/test change needed to
separate that reported approximation from strict representation, finiteness,
symmetry, convergence, and derivative/algebra failures. Additive consumers may
forward total and self/cross consistency diagnostics under their own active
IDs.

Exclusions: public driver/default behavior, corrected artifacts, solver
integration, exchange, EGOI, row-gauge substitutions, source or interaction
transforms, `C' V C`, and Cr2 production claims.

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

Status: implemented expert input/provenance facility.

Purpose: expose a low-cognitive-overhead expert knob for PQS/WL parent mapping
shape. Cr and Cr2 evidence shows that scanning only `core_spacing` along the
standard `s = sqrt(Z * core_spacing)` path is too restrictive for expert
consumers. The repo should expose the controlled scalar and record provenance;
it should not decide or tune the optimal value.

Approved user-facing convention:

- optional positive `s_factor`;
- default `s_factor = 1.0`;
- omitted `s_factor` and explicit `1.0` preserve current behavior;
- `standard_s = sqrt(Z * core_spacing)`;
- `effective_s = s_factor * standard_s`;
- one-center WL/atom mapping is literal:
  `AsinhMapping(c = core_spacing, s = effective_s)`;
- `core_spacing` remains the near-core physical scale.

Provenance must record:

- `mapping_s_factor`;
- `mapping_s_standard`;
- `mapping_s_effective`;
- `mapping_c` / `mapping_d` / `core_spacing` as already appropriate.

For multicenter PQS mapping, doer may apply the analogous per-center mapping
strength factor into the combined inverse-sqrt construction only if the
semantics are unambiguous. The implementation report must state exactly how
`s_factor` maps to the combined fit and how provenance records per-center or
per-axis effective values. If this is ambiguous, implement only the
one-center path and report the exact multicenter design question.

Approved source surface:

- `src/mappings.jl`;
- `src/pqs_source_box_route_driver_helpers.jl`;
- `src/cartesian_base_hamiltonian.jl`;
- `bin/cartesian_ham_builder.jl` only if needed for normal expert input;
- `src/cartesian_protected_ladder_bundle.jl` only to preserve/read recipe
  provenance.

Guardrails:

- this is an expert knob, not a default or optimization policy;
- do not add element-table defaults or automatic tuning;
- do not revive public `d`, public `parent_mapping_d`, public
  `parent_mapping_Z`, or route-specific mapping controls;
- do not reinterpret this as "smaller `core_spacing` is bad";
- do not change `Vee`, EGOI, rho0/P0, solver workflow, protected-localized
  convention, or residual/injection selection policy.

### HP-PQS-MAP-SFACTOR-TEST-01 — mapping `s_factor` validation

Status: implemented validation evidence.

Approved validation:

- `git diff --check`;
- package load;
- default H/H2 or small base artifact/readback path unchanged with omitted
  `s_factor`;
- explicit one-center atom with `s_factor != 1` records provenance and changes
  the mapping;
- small multicenter smoke if the multicenter path supports the knob;
- no Cr2 production run.

Failure rule: if multicenter combined-invsqrt mapping cannot unambiguously
support the same `s_factor` semantics, implement the one-center path only and
report the exact blocker before touching CR2 production scripts.

### HP-PQS-COULOMB-ACCURACY-FN-01 - producer-wide Coulomb accuracy policy

Status: compact/high producer policy implemented; fixed standard tier and
narrow canonical-driver exposure approved for implementation.

Canonical design: `coulomb_accuracy_policy.md`.

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
- delete `_cartesian_base_ida_hamiltonian(...)` if its focused caller scan
  remains empty, or require an explicit carried expansion argument if it is
  live. It must not independently select compact accuracy.

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

### HP-ROUTE-RECIPE-FN-01 — family-selective route recipe cleanup

Approved source files:

```text
src/pqs_source_box_route_driver_helpers.jl
src/cartesian_base_hamiltonian.jl
```

Approved behavior:

- `cartesian_recipe(route_inputs)` may construct only the subrecipe selected by
  `route_inputs.route_family`;
- for `route_family = :pqs_source_box`, route inputs must not require inactive
  `white_lindsey_*` fields; the produced recipe may set the inactive
  `white_lindsey` subrecipe to `nothing` while retaining the existing field
  name for caller compatibility;
- for `route_family = :white_lindsey_low_order`, explicit White-Lindsey route
  support must be preserved and the selected `white_lindsey` subrecipe must
  continue to be built from the existing WL route fields; the inactive
  `source_box` subrecipe may be `nothing` if no live WL caller requires it;
- `_cartesian_base_route(kind)` in `src/cartesian_base_hamiltonian.jl` may
  remove unused `white_lindsey_*` fields because the live base producer route
  uses `route_family = :pqs_source_box`;
- existing precomposed recipes that already provide `source_box` and
  `white_lindsey` fields may remain accepted if that compatibility path is
  still live, but it must not force new PQS-only route inputs to carry inactive
  WL vocabulary.

This ID preserves real WL/PQS algorithm differences while removing inactive WL
route-family fields from the current PQS base producer contract. It does not
approve canonical-driver changes, numerical kernel changes, terminal lowering
policy changes, shellification behavior changes, materialization or artifact
schema changes, route-stage diagnostics, status/report expansion, deletion of
WL materialization, or source files outside the two approved files.

Line budget: at most `80` added `src` lines, with net simplification expected.

Failure rule: if `cartesian_recipe(...)` cannot be made family-selective
without broader route-driver, report, materialization, or stage-object changes,
make no source commit and report the blocker.

### HP-ROUTE-RECIPE-TEST-01 — route recipe cleanup validation

Approved validation:

- `git diff --check`;
- package load;
- H atom/base artifact readback;
- H2 base artifact readback;
- compact H2 supplemented facade or driver path;
- focused route recipe smoke for explicit `:white_lindsey_low_order` if still
  practical, or a report of the exact live test/tool callers that block further
  WL route-input cleanup.

Existing committed tests may be adjusted only where they directly construct
route inputs that now no longer need inactive family fields. Known direct
route-recipe tests such as
`test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl` may drop inactive
WL fields if required by the source cleanup. No new committed test file, Cr2
run, driver workflow change, route diagnostic, or physics-reference scalar is
approved by this ID.

### HP-ROUTE-INV-FN-01 — retained-unit route inventory type cleanup

Approved source file:

```text
src/pqs_source_box_route_driver_helpers.jl
```

Approved cleanup targets:

- `_pqs_source_box_route_driver_named_tuple_from_units(...)`;
- runtime-keyed retained-unit inventory fields derived from unit labels,
  including `source_boxes`, `source_dimensions`, `retained_counts`, and
  `ranges`;
- runtime-keyed `pair_family_counts = NamedTuple{families}(...)`;
- same-file internal consumers that currently expect those runtime-keyed
  `NamedTuple` shapes.

Approved replacements:

- vector-backed records or tables with stable field names;
- stable dictionaries keyed by unit or pair-family labels only where lookup by
  label is genuinely needed;
- helper accessors that hide the storage shape from same-file callers;
- compact summaries that expose counts/order without encoding route size in the
  concrete type.

The retained-unit vector remains the ordered inventory authority. Unit labels
and pair-family labels may remain data values, but they must not become type
parameters.

This ID does not approve edits to `RawProductBoxPlan.source_mode_indices`,
`source_mode_column_indices`, `source_mode_indices(...)`,
`TerminalLoweringPlan.available_contracts`, `TerminalLoweringPlan.contracts`,
or `RetainedUnitTransformContractPlan.contracts`. It also does not approve
public input `NamedTuple` changes, fixed `NTuple{3,Int}` coordinate/dimension
changes, artifact sidecar table changes, numerical kernels, route recipe
behavior changes, shellification, terminal lowering, terminal basis, Residual
Gaussian, raw-block changes, canonical driver changes, Hamiltonian object
changes, matrix-key changes, reader changes, public API/export changes,
report/status/payload expansion, compatibility adapters, new committed tests,
Cr2 runs, or Cr2-specific workflow.

Line budget: at most `120` added `src` lines, with net simplification expected.
Failure rule: if the cleanup requires source files outside the approved file,
broader route/stage rewiring, public API changes, artifact changes, or an
adapter that preserves the old runtime-keyed type shape, make no source commit
and report the blocker.

### HP-ROUTE-INV-TEST-01 — retained-unit route inventory cleanup validation

Approved validation:

- `git diff --check`;
- package load;
- H2 base artifact write/readback through the existing reader;
- H2 supplemented artifact write/readback through the existing reader or
  canonical driver path;
- focused search confirming no `NamedTuple{unit_keys}` or
  `NamedTuple{families}` route inventory remains in
  `src/pqs_source_box_route_driver_helpers.jl`;
- no Cr2 run.

Existing committed tests may be adjusted only if they directly assert the old
runtime-keyed inventory shape. No new committed test file, Cr2 fixture,
driver-input fixture, benchmark, or route-diagnostic test is approved by this
ID.

### HP-RAW-SRCMODE-FN-01 — raw product source-mode inventory cleanup

Approved source files:

```text
src/cartesian_raw_product_sources/records.jl
src/cartesian_raw_product_sources/source_mode_indices.jl
src/cartesian_raw_product_sources/summaries.jl
```

Approved narrow consumer files, only as required by the storage change:

```text
src/cartesian_retained_unit_transform_contracts/unit_contracts.jl
src/cartesian_final_basis_realization/pqs_terminal_basis_realization.jl
src/cartesian_base_hamiltonian.jl
```

Approved cleanup targets:

- `RawProductBoxPlan.source_mode_indices::Tuple{Vararg{NTuple{3,Int}}}`;
- `RawProductBoxPlan.source_mode_column_indices::Tuple{Vararg{Int}}`;
- `source_mode_indices(...)` / source-mode summary accessors only to the extent
  required to hide vector-backed storage from approved callers;
- same-file and listed narrow consumers that currently depend on the
  tuple-backed storage shape.

Approved replacement:

- vector-backed source-mode coordinate storage;
- vector-backed source-mode column storage, or no stored column vector when the
  column sequence is exactly `1:count` and accessors provide the same ordered
  column numbers;
- stable accessors preserving deterministic source-mode order, mode values,
  length, indexing/iteration where currently used, and retained-rule parity.

The fixed `NTuple{3,Int}` coordinate and dimension values remain valid. The
variable-length source-mode inventory must not be encoded in `RawProductBoxPlan`
field types or accessor return types. Accessor compatibility means same facts
and order, not the old concrete `Tuple{Vararg{...}}` return shape.

This ID does not approve terminal-lowering `contracts` /
`available_contracts` tuple cleanup, retained-unit transform-contract tuple
cleanup outside the listed narrow consumer wiring, broad pair-block/source-box
rewrites, numerical kernel changes, route semantic changes, public API/export
changes, canonical driver changes, Hamiltonian object changes, matrix-key
changes, reader changes, artifact schema changes, route-stage objects,
report/status/payload expansion, persistent caches, compatibility layers that
preserve the old tuple-backed shape, new committed tests, Cr2 runs, or
Cr2-specific workflow.

Line budget: at most `150` added `src` lines, with net simplification expected.
Failure rule: if vectorizing the raw product plan requires source files outside
the approved surfaces, broad pair-block/source-box rewrites, public API or
artifact changes, numerical changes, or compatibility layers preserving the old
tuple-backed shape, make no source commit and report the exact caller/blocker.

### HP-RAW-SRCMODE-TEST-01 — raw product source-mode inventory validation

Approved validation:

- `git diff --check`;
- package load;
- H2 base artifact write/readback through the existing reader;
- H2 supplemented artifact write/readback through the existing reader;
- H2 R3 endpoint;
- focused raw-product source order and retained-rule parity;
- manifest source-mode and final-basis source-relation inspection;
- focused search confirming `RawProductBoxPlan` no longer stores source-mode
  inventories as `Tuple{Vararg{...}}`;
- no Cr2 run.

Existing committed tests may be adjusted only if they directly assert the old
tuple-backed source-mode inventory shape. No new committed test file, Cr2
fixture, driver-input fixture, benchmark, or route-diagnostic test is approved
by this ID.

### HP-CONTRACT-VEC-FN-01 — contract-plan vector cleanup

Approved source files:

```text
src/cartesian_terminal_lowering/contracts.jl
src/cartesian_terminal_lowering/selection.jl
src/cartesian_terminal_lowering/summaries.jl
src/cartesian_retained_unit_transform_contracts/records.jl
src/cartesian_retained_unit_transform_contracts/unit_contracts.jl
src/cartesian_retained_unit_transform_contracts/summaries.jl
```

Approved narrow consumer files, only as required by the storage change:

```text
src/pqs_source_box_route_driver_helpers.jl
src/cartesian_base_hamiltonian.jl
src/cartesian_final_basis_realization/pqs_terminal_basis_realization.jl
```

Approved cleanup targets:

- `TerminalLoweringPlan.available_contracts::Tuple{Vararg{TerminalLoweringContract}}`;
- `TerminalLoweringPlan.contracts::Tuple{Vararg{TerminalLoweringContract}}`;
- `RetainedUnitTransformContractPlan.contracts::Tuple{Vararg{RetainedUnitTransformContract}}`;
- same-file and listed narrow consumers that currently depend on those
  tuple-backed plan field shapes.

Approved replacement:

- vector-backed terminal-lowering available-contract storage;
- vector-backed terminal-lowering selected-contract storage;
- vector-backed retained-unit transform-contract storage;
- stable accessors preserving ordered contract facts and current behavior:
  `available_contracts(plan)`, `selected_contracts(plan)`, `contracts(plan)`,
  and `transform_contracts(plan)`.

Accessor compatibility means same ordered facts, iteration order, selected
contract behavior, transform-contract behavior, and summaries. It does not
mean preserving variable-length `Tuple` concrete field types or accessor return
types.

This ID does not approve changing
`source_cpbs::Tuple{Vararg{CoordinateProductBox}}`, raw product source-mode
storage, retained-unit route inventories, public input `NamedTuple`s, fixed
coordinate/product-box value objects, numerical kernels, route semantic
changes, shellification behavior changes, public API/export changes, canonical
driver changes, Hamiltonian object changes, matrix-key changes, reader changes,
artifact/manifest schema changes, route-stage objects,
report/status/payload expansion, persistent caches, compatibility layers that
preserve the old tuple-backed plan field types, new committed tests, Cr2 runs,
or Cr2-specific workflow.

Line budget: at most `150` added `src` lines, with net simplification expected.
Failure rule: if vectorizing the plan inventories requires source files outside
the approved surfaces, broad route/stage rewrites, public API or artifact
changes, numerical changes, or compatibility layers preserving the old
tuple-backed plan field types, make no source commit and report the exact
caller/blocker.

### HP-CONTRACT-VEC-TEST-01 — contract-plan vector cleanup validation

Approved validation:

- `git diff --check`;
- package load;
- H2 base artifact write/readback through the existing reader;
- H2 supplemented artifact write/readback through the existing reader;
- H2 R3 endpoint;
- focused terminal-lowering contract order parity;
- focused retained-unit transform-contract order parity;
- focused search confirming targeted plan inventories no longer store
  contracts as `Tuple{Vararg{...}}`;
- no Cr2 run.

Existing committed tests may be adjusted only if they directly assert the old
tuple-backed contract-plan field shape. No new committed test file, Cr2
fixture, driver-input fixture, benchmark, or route-diagnostic test is approved
by this ID.

### HP-ROUTE-STAGE-TYPE-FN-01 — route/stage type-surface cleanup

Approved source files:

```text
src/pqs_source_box_route_driver_helpers.jl
src/cartesian_terminal_shellification_geometry.jl
```

Approved cleanup targets:

- `_pqs_source_box_route_driver_terminal_lowering_contract_inventory_from_plan`;
- `cartesian_units` route/stage return surfaces that carry oversized
  compatibility inventories;
- `_pqs_source_box_route_driver_transform_stage_low_order_summary`;
- `cartesian_transforms` route/stage return surfaces that carry oversized
  compatibility inventories;
- `_cartesian_terminal_shellification_region_unit_inventory`;
- related terminal-region lowering inventory summary surfaces in
  `src/cartesian_terminal_shellification_geometry.jl` only where the same
  runtime-sized type-surface pattern appears.

Approved replacement/deletion shapes:

- delete stale route/stage compatibility inventories with no active approved
  caller;
- replace remaining runtime-sized `NamedTuple` / `Tuple` carriers with
  vector-backed compact internal objects, stable dictionaries, accessors, or
  smaller summaries;
- shrink wide internal stage return signatures only where all live approved
  callers can be updated within the approved source files;
- preserve deterministic terminal shellification/lowering order and existing
  behavior.

Required preservation:

- H2 base artifact/readback behavior;
- H2 supplemented artifact/readback behavior;
- deterministic terminal shellification/lowering order;
- existing public driver contract;
- existing artifact schema and manifest behavior;
- existing numerical matrices.

This ID does not approve source files outside the approved boundary, driver
changes, artifact schema or manifest changes, public API/export changes,
numerical kernel changes, matrix value changes, raw-block changes, Residual
Gaussian/MWG/IDA semantic changes, route semantic changes, shellification
behavior changes, route diagnostic/status/report expansion, broad route-stage
redesign, new public contracts, PackageCompiler/PrecompileTools/sysimage or
precompile workload work, new committed tests, Cr2 runs, or Cr2-specific
workflow. No compatibility adapter may preserve the old runtime-sized type
surface merely under a new name.

Line budget: at most `200` added `src` lines, with net simplification expected.
Failure rule: if cleanup requires source files outside the approved boundary,
broad route-stage redesign, new public contracts, artifact changes, numerical
changes, or a precompile/sysimage mechanism, make no source commit and report
the exact blocker.

### HP-ROUTE-STAGE-TYPE-TEST-01 — route/stage type-surface cleanup validation

Approved validation:

- `git diff --check`;
- package load;
- H2 base artifact write/readback through the existing reader;
- H2 supplemented artifact write/readback through the existing reader;
- H2 R3 endpoint if the pass touches terminal realization behavior;
- focused terminal shellification/lowering order parity;
- focused scan for newly introduced `NamedTuple{...}`, variable-size
  `Tuple(...)`, `Tuple{Vararg{...}}`, and runtime-keyed inventories in the
  approved files;
- optional Be2 q5 compile/timing comparison after correctness passes;
- no Cr2 run.

Existing committed tests may be adjusted only if they directly assert the old
stale compatibility inventory shape. No new committed test file, Cr2 fixture,
driver-input fixture, benchmark, route-diagnostic test, or precompile workload
is approved by this ID.

### HP-ROUTE-STAGE-CARRIER-FN-01 — route/stage carrier cleanup

Approved source files:

```text
src/pqs_source_box_route_driver_helpers.jl
src/pqs_source_box_diatomic_complete_core_shell.jl
```

Optional source file, only if directly required to slim the approved path:

```text
src/cartesian_final_basis_realization/pqs_terminal_basis_realization.jl
```

Approved cleanup targets:

- `cartesian_shells` stage carrier and return signature;
- `cartesian_units` stage carrier and return signature;
- `cartesian_transforms` stage carrier and return signature;
- terminal topology support-region planning in
  `src/pqs_source_box_diatomic_complete_core_shell.jl`;
- terminal retained-rule planning in
  `src/pqs_source_box_diatomic_complete_core_shell.jl`;
- terminal realization plan carriers in
  `src/cartesian_final_basis_realization/pqs_terminal_basis_realization.jl` only
  where directly required to avoid reintroducing a large stage-carried plan
  shape through the approved route/stage path.

Approved replacement/deletion shapes:

- stop carrying giant shellification, route-skeleton, support-plan,
  retained-rule-plan, and terminal-plan `NamedTuple` / tuple shapes across the
  approved stage function signatures;
- replace necessary carriers with compact typed/vector-backed records, stable
  dictionaries, accessors, or smaller summaries;
- recompute small derived summaries from canonical objects inside the approved
  path where that is simpler than carrying wide stage payloads;
- delete stale compatibility carriers with no active approved caller;
- preserve deterministic terminal support, shellification, and lowering order.

Route skeleton construction semantics are not changed by this ID. Edits to
`src/pqs_source_box_route_driver_skeletons.jl` are not approved.

Required preservation:

- H2 base artifact/readback behavior;
- H2 supplemented artifact/readback behavior;
- H2 R3 endpoint if terminal realization is touched;
- deterministic terminal support/shellification/lowering order;
- existing public driver contract;
- existing artifact schema and manifest behavior;
- existing route semantics and numerical matrices.

This ID does not approve source files outside the approved boundary, edits to
`src/pqs_source_box_route_driver_skeletons.jl`, driver changes, artifact schema
or manifest changes, public API/export changes, numerical kernel changes,
matrix value changes, raw-block changes, Residual Gaussian/MWG/IDA semantic
changes, route semantic changes, shellification behavior changes, route
diagnostic/status/report expansion, broad route-stage redesign, new public
contracts, PackageCompiler/PrecompileTools/sysimage or precompile workload
work, new committed tests, Cr2 runs, or Cr2-specific workflow. No compatibility
adapter may preserve the old runtime-sized carrier merely under a new name.

Line budget: at most `250` added `src` lines, with net simplification expected.
Failure rule: if cleanup requires source files outside the approved boundary,
broad route-stage redesign, public API changes, artifact changes, numerical
changes, or precompile/sysimage machinery, make no source commit and report the
exact blocker.

### HP-ROUTE-STAGE-CARRIER-TEST-01 — route/stage carrier cleanup validation

Approved validation:

- `git diff --check`;
- package load;
- H2 base artifact write/readback through the existing reader;
- H2 supplemented artifact write/readback through the existing reader;
- H2 R3 endpoint if terminal realization is touched;
- focused terminal support/shellification/lowering order parity;
- focused scan for newly introduced runtime-sized `NamedTuple{...}`,
  `Tuple(...)`, `Tuple{Vararg{...}}`, and runtime-keyed inventories in the
  approved files;
- optional Be2 q5 post-cleanup compile/timing comparison after correctness
  passes;
- no Cr2 run.

Existing committed tests may be adjusted only if they directly assert the old
stale carrier shape. No new committed test file, Cr2 fixture,
driver-input fixture, benchmark, route-diagnostic test, or precompile workload
is approved by this ID.

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

## Approved For R1 One-Center Base Atoms

This section approves only the explicit origin-centered one-center
all-electron atom relaxation recorded in `r1_one_center_base_atoms.md`. It
extends the existing base facade scope without adding a new public function,
new export, new artifact schema, new route vocabulary, or supplemented atom
authority.

### HP-R1-ATOM-FN-01 — explicit one-center all-electron base atom facade

Approved source file:

```text
src/cartesian_base_hamiltonian.jl
```

Approved behavior:

- accept exactly one origin-centered atom in the existing
  `cartesian_base_hamiltonian(system; basis, hamfile)` call shape;
- require explicit vector-valued `atom_symbols`, `nuclear_charges`,
  `atom_locations`, and explicit integer `nup`, `ndn`;
- require finite positive integer-valued nuclear charge supplied by the caller;
- require neutral all-electron count
  `nup + ndn == round(Int, only(nuclear_charges))`;
- treat the atom symbol as provenance/user labeling only, not as a source of
  charge, spin, basis, or ECP defaults;
- keep required one-center basis fields `ns`, `core_spacing`, and `radius`
  after `HP-COMP-NS-*` normalization;
- treat public `d`, if temporarily accepted, as a deprecated compatibility
  alias that must equal resolved `core_spacing`.

This ID does not approve translated atoms, element lookup/default tables,
inferred charge or spin, ECP, solver workflow, supplemented atom construction
under the base facade, public API redesign, or new artifact fields. Supported
supplemented one-center atoms are governed separately by `HP-COMP-SUPPATOM-*`.

### HP-R1-ATOM-WIRE-01 — one-center atom shared workflow wiring

Approved source file:

```text
src/cartesian_base_hamiltonian.jl
```

Approved behavior:

- map public `only(system.nuclear_charges)` to the existing private
  White-Lindsey atomic mapping `Z`;
- map the resolved public `core_spacing` to the private White-Lindsey
  `parent_mapping_d`;
- keep `reference_spacing`, `tail_spacing`, and box/domain controls separate
  from `core_spacing`;
- feed atom geometry/shellification normalization into the same terminal-basis,
  one-body, IDA, `CartesianIDAHamiltonian`, artifact-writing, and provenance
  machinery used by the base producer;
- preserve existing `HP-R1-ART-01` `producer_provenance/` keys with
  `route = :one_center_pqs_base`.

Atoms and diatomics must share the same producer workflow after the narrow
geometry/shellification differences. This ID does not approve an atom-only
Hamiltonian builder, parallel atom materialization path, atom route-stage
object, atom report/status/payload object, or metadata/provenance carrier used
as algorithmic data.

Line budget for `HP-R1-ATOM-FN-01` plus `HP-R1-ATOM-WIRE-01`: at most `80`
added `src` lines. If implementation needs source edits outside
`src/cartesian_base_hamiltonian.jl`, changes to private materialization
owners, atom-only materialization, new artifact keys, translated atoms,
ECP behavior, solver workflow, element lookup/default tables, committed
fixtures/tests, route/report/status/payload expansion, or supplemented atom
work outside `HP-COMP-SUPPATOM-*`, stop and request a new docs-only amendment.

### HP-R1-ATOM-TEST-01 — one-center base atom validation

Approved validation:

- existing origin-centered H public facade endpoint remains unchanged,
  now expressed as `core_spacing = 0.3`, `reference_spacing = 1.0`, and
  internal `parent_mapping_d = core_spacing`;
- optional ignored/user-run Be or Cr one-center base atom artifact
  write/readback using explicit charge, spin sectors, origin geometry, and
  basis controls;
- finite/symmetric `K`, unit `U_A`, and IDA `V` for ignored/user-run non-H
  atom checks;
- clear `ArgumentError` for translated atom input, mismatched temporary `d`,
  noninteger or nonpositive charge, nonneutral electron count, or
  element-table/default requests where practical.

No new committed test file, committed non-H atom fixture, public non-H
reference scalar, solver run, supplemented atom endpoint under this base-atom
lane, ECP gate, translated-atom gate, or driver change is approved by this ID.
Supported supplemented atom validation is governed by
`HP-COMP-SUPPATOM-TEST-01`.

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

## Approved For Residual Gaussian Domain Migration

Current Residual Gaussian domain algorithm authority is
`residual_gaussian_domain_module.md`. This registry records the approved module
files and function surfaces only.

The Residual Gaussian module does not own compact supplemented artifact writing
or `supplement_provenance/`. Artifact/facade hooks remain outside RG unless a
later amendment names a real duplication or consumer reason.

### HP-RG-FILE-01 — Residual Gaussian module files

Approved internal module and files:

```text
src/cartesian_residual_gaussians/CartesianResidualGaussians.jl
src/cartesian_residual_gaussians/residual_basis.jl
src/cartesian_residual_gaussians/augmented_operators.jl
src/cartesian_residual_gaussians/mwg_interaction.jl
```

No public export is approved.

### HP-RG-OBJ-01 — residual Gaussian basis object

Approved domain object: a numerical residual Gaussian basis object carrying base
dimension, candidate count, residual dimension, candidate owner indices,
residual source owner indices, owner retained counts, residual occupations,
cutoff/tolerance policy, selection/orientation/sign rules, and final
`T_G::Matrix{Float64}` / `T_A::Matrix{Float64}` transforms.

It must not be a status/result payload and must not carry route metadata,
report fields, status flags, artifact data, MWG descriptors, or public API
state.

### HP-RG-FN-01 — residual Gaussian basis construction

Approved production name:

```julia
build_residual_gaussian_basis(...)
```

Candidate owner indices are required. The current algorithm is the owner-local
residual occupation construction defined in
`residual_gaussian_domain_module.md`.

### HP-RG-FN-02 — exact augmented operator transformation

Approved production name:

```julia
transform_augmented_operator(...)
```

This owns exact `[G,A] -> [G,R]` transformation for kinetic, uncharged
by-center nuclear attraction, and first/second Cartesian moment matrices.

### HP-RG-FN-03 — moment-matched Gaussian descriptors

Approved production name:

```julia
moment_matched_gaussians(...)
```

Descriptors are computed from final merged residual functions and exact moment
matrices. They are residual-containing interaction descriptors only, not exact
residual-GTO Coulomb integrals.

### HP-RG-FN-04 — residual IDA interaction assembly

Approved production name:

```julia
assemble_residual_ida_interaction(...)
```

This owns residual-containing MWG/IDA blocks `V_GM` and `V_MM`, combined with
unchanged base `V_GG`. `V_GM` uses weight-aware final-basis density
normalization for PQS shell blocks.

### HP-RG-WIRE-01 — migration from terminal residual file

`src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl` may keep
small compatibility wrappers for live callers and artifact/facade hooks. Moved
physics helpers should be deleted from that file once callers use RG-domain
helpers directly.

### HP-RG-TEST-01 — migration validation

Approved validation is the existing standalone H2 residual-GTO/MWG endpoint
with augmented dimension `505`, self-Coulomb `0.4574161883692301`, exact
one-body/moment checks, independent weight-aware `V_GM` check, and optional
ignored Be2 usability/performance measurement when a source pass changes the
interaction path or facade wiring.

No new committed test file, Be2 committed gate, or Cr2 full
Hamiltonian/artifact/facade validation is approved.

### HP-RG-ORTHO-FN-01 — residual final-orthogonality robustness

Approved source files:

```text
src/cartesian_residual_gaussians/residual_basis.jl
src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl
```

`residual_basis.jl` is the primary owner. The terminal residual file is
approved only for narrow existing internal keyword plumbing if an approved
tolerance/check option must be passed through a compatibility entry point.

Approved target: make final residual orthogonalization and final
`R' S R - I` identity validation robust for small floating-point overshoots
with healthy owner-local selection and final merge spectra. The motivating
strict N2 q5 p10 case at `core_spacing = 0.042857` has `max |G' S R| =
1.776e-14`, `max |R' S R - I| = 1.673e-10`, retained counts `9,9`, and final
merge eigenvalues `7.232e-2 .. 1.928`.

Allowed changes:

- stable symmetric final merge normalization/check;
- explicitly symmetrized final residual-overlap validation;
- combined absolute/relative final identity check
  `err_RR <= 1.0e-10 + 1.0e-10 * max(1, scale_RR)`;
- no public API unless the option is already routed through internal keywords.

This ID does not approve blind broad tolerance relaxation, residual-selection
semantic changes, global residual selection, occupation-cutoff changes,
negative-eigenvalue tolerance changes, final merge eigenvalue flooring, width
filtering as a conditioning repair, MWG/IDA/nuclear/raw-block changes, artifact
schema changes, driver changes, status/report fields, public API/export
changes, new committed tests, Cr2 workflow, or source files outside the two
approved files.

Failure rule: if the strict N2 case requires changing residual selection,
supplement construction, or final-basis construction, make no source commit and
report the blocker.

### HP-RG-ORTHO-TEST-01 — residual final-orthogonality validation

Approved validation:

- existing H2 residual-GTO/MWG endpoint;
- H2 base/supplemented readback if touched through the facade or compatibility
  file;
- ignored strict N2 q5 p10 residual audit or artifact smoke at
  `core_spacing = 0.042857`;
- one passing N2 comparison at `core_spacing = 0.05` or `0.075`;
- report `max |G' S R|`, `max |R' S R - I|`, retained owner counts, and final
  merge eigenvalue range/condition.

No committed fixture/test, Cr2 full Hamiltonian, Cr2 artifact, Cr2 facade
support, driver workflow, artifact schema change, solver/RHF, ECP, or EGOI
work is approved.

### HP-RG-IDTOL-FN-01 — residual final-identity tolerance default

Approved source files:

```text
src/cartesian_residual_gaussians/residual_basis.jl
src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl
```

`residual_basis.jl` is the primary owner. The terminal residual file is
approved only for narrow compatibility keyword default plumbing if needed.

Approved target: set the default final residual `R' S R` identity validation
tolerance to `1.0e-8` when owner-local selection, final merge metric checks,
and `G' S R` orthogonality remain healthy. This is a final
validation/cleanup tolerance only. It is not a residual direction-selection
criterion.

Evidence: Be atom cc-pV5Z `lmax = 1`, `ns = 5`,
`core_spacing = 0.075`, `radius = 20`, `Z = 4`, `nup = 2`, `ndn = 2`
failed only with `max |R' S R - I| = 2.183e-10` against an allowed error of
about `2.000e-10`. The same run had retained residual count `21`, minimum
retained occupation `6.151e-6`, final merge condition `1.0`, and
`max |G' S R| = 1.776e-14`.

Required policy:

- keep the then-current `residual_occupation_cutoff = 1.0e-8` for the Be
  tolerance pass only; this older production default is superseded by
  `HP-RG-CUTOFF-FN-01`;
- keep width/zeta filtering explicit and user-controlled;
- keep owner-local metric checks, final merge metric checks, and `G' S R`
  orthogonality checks active;
- do not drop retained directions to pass final identity validation.

This ID does not approve driver changes, artifact schema/provenance/reader/
manifest changes, residual-selection algorithm changes, width/zeta filtering
default changes, owner grouping changes, merge metric failure-rule changes,
MWG/IDA convention changes, Gaussian raw-block changes, terminal-basis changes,
WL/PQS route changes, shellification changes, Hamiltonian assembly changes,
committed tests/fixtures, Cr2 workflow, or source files outside the two
approved files.

Failure rule: if Be cc-pV5Z cannot pass by changing only the final identity
tolerance default, make no source commit and report the exact blocker.

### HP-RG-IDTOL-TEST-01 — residual final-identity tolerance validation

Approved validation:

- Be atom cc-pV5Z `lmax = 1` residual audit/artifact path passes with the same
  `21` retained residual directions;
- Be atom cc-pVDZ `lmax = 1` still passes;
- H2 residual-GTO/MWG endpoint remains unchanged;
- report `max |R' S R - I|`, allowed tolerance, retained count, minimum
  retained occupation, final merge condition, and `max |G' S R|`;
- no Cr2 run.

No committed fixture/test, driver workflow, artifact schema change,
solver/RHF, ECP, EGOI, Cr2 full Hamiltonian, Cr2 artifact, or Cr2 facade
support is approved.

### HP-RG-CUTOFF-FN-01 — residual occupation cutoff and identity tolerance defaults

Approved source files:

```text
src/cartesian_residual_gaussians/residual_basis.jl
src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl
```

`residual_basis.jl` is the primary owner. The terminal residual file is
approved only for narrow compatibility keyword default plumbing if needed.

Approved target: set the default Residual Gaussian owner-local residual
selection cutoff and final identity validation tolerance to:

```text
residual_occupation_cutoff = 5.0e-8
identity_atol = 5.0e-8
```

This was the production default after the Cr atom marginal-direction pass. The
residual occupation cutoff is superseded by `HP-RG-CUTOFF-FN-02`; the final
identity validation tolerance remains `5.0e-8`.

Evidence: Cr atom PQS `basis_ns = 9`, `map_ns = 11`, `lmax = 1` retained a
marginal residual direction at occupation `3.637e-8`. If the production policy
is to discard such marginal residual directions, the default cutoff must say
so explicitly rather than preserving the direction through the older
`1.0e-8` cutoff.

Approved behavior for this now-superseded default:

- discard owner-local residual directions below `5.0e-8` by default;
- keep the final `R' S R` identity validation default aligned at
  `identity_atol = 5.0e-8`;
- preserve explicit caller overrides where already supported;
- keep width/zeta filtering explicit and user-controlled;
- keep owner-local metric checks, negative-eigenvalue tolerances, final merge
  metric checks, and `G' S R` orthogonality checks active;
- keep owner grouping, final merge failure rules, MWG/IDA conventions,
  artifact schema/provenance/reader/manifest, driver workflow, and public API
  unchanged.

This ID supersedes the older `1.0e-8` production defaults recorded under
`HP-RG-IDTOL-FN-01`. It does not change the residual-selection algorithm; it
changes the default retained-occupation policy.

Forbidden:

- residual selection algorithm changes;
- global raw-candidate residual selection;
- global raw-column pivoted-Cholesky residual selection;
- owner grouping changes;
- negative-eigenvalue tolerance changes;
- final merge metric failure-rule changes or merge eigenvalue flooring;
- width/zeta filtering default changes or width filtering as conditioning
  repair;
- MWG/IDA, nuclear, raw-block, exact-operator, terminal-basis, WL/PQS route,
  shellification, Hamiltonian assembly, artifact schema, reader, manifest,
  driver, public API/export, solver/RHF, ECP, EGOI, or Cr2 workflow changes;
- committed fixtures/tests except the exact existing H2 endpoint assertion
  update named under `HP-RG-CUTOFF-TEST-01`.

Failure rule: if the Cr atom case cannot pass or cleanly drop the marginal
direction by changing only the two approved defaults, make no source commit in
the later implementation pass and report the exact blocker.

### HP-RG-CUTOFF-TEST-01 — residual cutoff/tolerance validation

Approved validation:

- Cr atom PQS `basis_ns = 9`, `map_ns = 11`, `lmax = 1` residual construction
  passes or cleanly drops the marginal `s4` direction at occupation
  `3.637e-8` as intended;
- Be atom cc-pV5Z still passes;
- H2 residual-GTO/MWG endpoint remains unchanged;
- exactly update the existing committed H2 endpoint test
  `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl` so both cutoff
  assertions expect `5.0e-8` instead of `1.0e-8`: the in-memory
  `residual.occupation_cutoff` assertion and the artifact/provenance
  `values[:occupation_cutoff]` assertion;
- report retained counts, minimum retained occupation, `max |G' S R|`,
  `max |R' S R - I|`, allowed tolerance, and final merge condition;
- no Cr2 run.

No other committed fixture/test, driver workflow, artifact schema change,
solver/RHF, ECP, EGOI, Cr2 full Hamiltonian, Cr2 artifact, or Cr2 facade
support is approved.

### HP-RG-CUTOFF-FN-02 — production residual cutoff tightening

Status: implemented current production default.

Approved source files:

```text
src/cartesian_residual_gaussians/residual_basis.jl
src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl
```

`residual_basis.jl` is the primary owner. The terminal residual file is
approved only for narrow compatibility keyword default plumbing if needed.

Approved target:

```text
residual_occupation_cutoff = 1.0e-6
identity_atol = 5.0e-8
```

Evidence: Cr2 residual spectra show the worst low-H1 modes are built from
marginal owner-local residual directions with occupations around
`1.27e-7` to `8.98e-7`. The current `5.0e-8` cutoff retains those directions;
`1.0e-6` drops `6` directions per owner in the cited audit. Broad residual
widths still matter, but broad width alone is not the first production rule
because one-center atoms can have broad candidates without the same bad
`H1_RR` sector.

Approved behavior:

- set the default owner-local residual occupation cutoff to `1.0e-6`;
- keep the final `R' S R` identity validation default at
  `identity_atol = 5.0e-8`;
- preserve explicit caller overrides where already supported;
- keep width/zeta filtering explicit and user-controlled;
- keep owner-local metric checks, negative-eigenvalue tolerances, final merge
  metric checks, and `G' S R` orthogonality checks active;
- keep owner grouping, final merge failure rules, MWG/IDA conventions,
  artifact schema/provenance/reader/manifest, driver workflow, and public API
  unchanged.

This ID supersedes only the `HP-RG-CUTOFF-FN-01` residual occupation default.
It does not change the residual-selection algorithm, identity tolerance,
negative-eigenvalue tolerances, or merge policy.

Forbidden:

- residual selection algorithm changes;
- kinetic or `H1_RR` spectral guards;
- global raw-candidate residual selection;
- global raw-column pivoted-Cholesky residual selection;
- owner grouping changes;
- negative-eigenvalue tolerance changes;
- final merge metric failure-rule changes or merge eigenvalue flooring;
- width/zeta filtering default changes or width filtering as conditioning
  repair;
- MWG/IDA, nuclear, raw-block, exact-operator, terminal-basis, WL/PQS route,
  shellification, Hamiltonian assembly, artifact schema, reader, manifest,
  driver, public API/export, solver/RHF, ECP, EGOI, Cr2 workflow, Cr2 artifact,
  or full HF changes;
- committed fixtures/tests except the exact existing H2 endpoint assertion
  update named under `HP-RG-CUTOFF-TEST-02`.

Failure rule: if the Cr2 residual-only audit does not cleanly remove the
identified low-occupation directions or if low-`H1_RR` ghost modes remain
after the cutoff change, make no further source changes in this lane. Report
the residual spectra and request separate kinetic/`H1_RR` spectral-guard
authority if needed.

### HP-RG-CUTOFF-TEST-02 — production residual cutoff validation

Status: implemented validation evidence.

Approved validation:

- Cr2 residual-only audit, not full HF or a Cr2 artifact/workflow run;
- Cr2 owner retained counts should drop from `68 + 68` to `62 + 62`;
- recompute and report residual spectra:
  - `min eig(K_RR)`;
  - `min eig(H1_RR)`;
  - low-mode candidate composition;
- if low-H1 ghost modes remain, stop and request separate kinetic/`H1_RR`
  spectral-guard authority;
- Be high-zeta residual construction still passes;
- H2 residual-GTO/MWG endpoint remains unchanged except for the approved
  default cutoff/provenance value;
- exactly update the existing committed H2 endpoint test
  `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl` so both cutoff
  assertions expect `1.0e-6` instead of `5.0e-8`: the in-memory
  `residual.occupation_cutoff` assertion and the artifact/provenance
  `values[:occupation_cutoff]` assertion.

No other committed fixture/test, driver workflow, artifact schema change,
solver/RHF, ECP, EGOI, full HF, Cr2 full Hamiltonian, Cr2 artifact, or Cr2
facade support is approved.

### HP-RG-NUMCOMP-FN-01 — numerical-complete residual basis and additive consumer

Status: implemented internal opt-in facility at `b2da7070c`.

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

Canonical record:
[Default-off direct-G residual injection](residual_gaussian_injection_hybrid.md).

Evidence: ignored Cr/Cr2 probes and manager running-log history established
that direct-`G` replacement was well conditioned but did not remove the
tested low two-owner residual sector. The audit authorizes no source,
artifact, workflow, or renewed measurement work.

### HP-RG-INJECT-FN-01 — default-off direct-G injection compatibility

Status: implemented preservation-only internal compatibility facility.

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

Purpose: established mandatory occupied protection and capture-spectrum
selection on Be/Ne before the source-backed helper was approved. Evidence is
retained in manager running-log Passes 323-324 and summarized by the canonical
contract below.

### HP-RG-OCC-FIRST-INJECT-FN-01 / HP-RG-OCC-FIRST-INJECT-TEST-01 — occupied-first injection geometry

Status: implemented internal facility.

Owner: `CartesianResidualGaussians` occupied-first geometry in
`residual_basis.jl`.

Canonical contract:

- [Occupied-first injection geometry](occupied_first_injection.md)

Approved source surface:

- `src/cartesian_residual_gaussians/residual_basis.jl`;
- optional read-only consumption of already-owned coefficients from
  `src/cartesian_reference_density/atomic_hf_reference_packets.jl` and
  `src/cartesian_external_gto_import.jl`, without changing either contract.

Approved test surface:

- `test/misc/runtests.jl` for the tiny pre/post and malformed-capture contract;
- `test/nested/cartesian_occupied_first_injection_runtests.jl` for the real
  bounded Be/Ne packet-driven PQS gate and terminal due diligence.

Dependencies:

- `HP-PQS-ATOMREF-PACKET-*` or `HP-REP-XGTO-IMPORT-*` supplies identified
  occupied coefficients;
- `HP-RG-PROTECT-ADDREF-*` owns the separate protected-localized consumer over
  `M = [G, R_compact]`.

Permission summary: validate physical mixed-overlap geometry, make supplied
`Y_occ` mandatory, distinguish pre-inclusion capture from post-inclusion
recovery, and select/report optional supplement directions by capture cutoff.
Weak rejected directions do not become MWG residual channels.

Current boundary: `occupied_first_injection_geometry(...)` is implemented and
tested but is not wired into the protected-localized builder and is not a
direct substitute for staged protected-original geometry.

Non-goals: screened-Hartree changes, protected-builder composition under these
IDs, EGOI, shell-local injection, fake-RDM changes, artifacts, public workflow,
automatic defaults, solver, exchange, row-gauge rho0/P0, or Cr/Cr2 claims.

### HP-RG-PROTECT-ADDREF-FN-01 / HP-RG-PROTECT-ADDREF-TEST-01 - protected additive atomic reference correction

Status: implemented narrow internal, opt-in facility. This is the first
real protected-localized consumer of occupied-first reference geometry. It is
not public workflow, artifact, solver, or production Cr2 authority.

Purpose: build a protected-localized homonuclear molecular member whose basis
contains the full span of all placed converged atomic packet occupied spaces,
then assemble the existing screened direct-Hartree correction in native `L`
order.

This path is internal and explicitly opt-in. Protected members built without
placed reference packets must preserve current geometry, `H1_L`, `Vee_L`,
native ordering, and artifact behavior.

The intended private composition seam is conceptually:

```text
_plb_build_additive_reference_member(recipe, stages, placements)
    -> (member, correction)
```

Each normalized placement carries `packet`, `owner_index`, `center`, and
explicit `supplement_indices`. It must share the current member-build core;
`_plb_build_member(recipe, stages)` remains the unchanged no-reference path.

Required construction:

1. Build ordered compact-first `R_compact` once and define
   `M = [G,R_compact]`.
2. Validate each packet's stored overlap fingerprint exactly, then embed its
   original occupied block `Y_a` by exact owner-local atom/basis, count,
   indices, placement, labels, angular powers, and column order. Accept the
   mapped overlap block only when its matrix infinity-norm difference from the
   packet overlap is at most the unchanged `1e-10`; mapped raw-byte hash
   equality is diagnostic, not required.
3. Form a full-rank `S_AA`-orthonormal union `Y_union` for basis protection.
4. Make `Y_union` mandatory first, then add current compact-original protected
   directions after orthogonalizing them against that union.
5. Apply the existing staged representability and fake-RDM policy only to the
   remaining optional supplement complement.
6. Preserve every original `Y_a` and occupation vector separately. The
   orthonormalized union is not used to define `P0`.
7. Represent each original packet block in native `L` order by final-basis
   cross overlap and form `P0_L = sum_a C_aL*n_a*C_aL'` and
   `q0_L = diag(P0_L)`. Validate every packet block separately; do not globally
   orthogonalize packet blocks when forming `P0_L`.
8. Build and sum placed ordinary fitted-potential `GG/GA/AA` Hartree blocks,
   reject packets carrying retired moment-polish provenance, transform through
   the existing protected fixed-sector helper, then localize
   `J0_L = W' * J0_F * W`.
9. Build the compact-convention no-half reference Coulomb energy
   `E0 = sum_a E_aa + 2*sum_{a<b}E_ab`, including explicit cross terms.
10. Call the existing screened-Hartree core with native-order `Vee_L`, `J0_L`,
    `P0_L/q0_L`, and `E0` to return the existing
    `ScreenedHartreeCorrection`. Do not duplicate the formula in the ladder
    owner.

The staged protected geometry must consume the exact already-built residual
object rather than reconstruct compact selection. If needed, this lane permits
one internal vector-backed field on `CartesianResidualGaussianBasis`:

```text
compact_source_candidate_indices::Union{Nothing,Vector{Int}}
```

It is the sorted unique accepted-source set for ordered compact-first MGS, not
a one-to-one label for final residual columns after merge cleanup. It is
`nothing` for selection rules without native accepted-source semantics. Do not
parse labels. Delete or delegate the current duplicate compact-selection block
in staged protected geometry. No artifact field is approved.

Mandatory occupied directions are not optional cutoff candidates. If the
union is not full-rank/recoverable or is not stably representable by `M` under
the active staged representability threshold, stop as insufficient
compact-main-basis support. Do not discard the direction, convert it into an
RG, or weaken the threshold.

Reference algebra is additive:

```text
P0 = sum_a P_a
q0 = diag(P0)
J0 = sum_a J_a
E0 = sum_a E_aa + 2*sum_{a<b} E_ab
```

The original per-packet occupied blocks define `P0`. Packet density fits define
self/cross energies. Packet fitted potentials are only the fast representation
of the same placed reference field.

Report `Tr(P0*J0_fit)-E0_fit` as a fitted-potential consistency approximation.
For each packet report `Tr(P_a*J_a)-E_aa`; for each pair report
`Tr(P_a*J_b)+Tr(P_b*J_a)-2E_ab`. The decomposition must sum to the total within
numerical assembly tolerance, but its magnitude is not a `1e-8 Ha` rejection
gate. Exact/density-fit oracle identities remain strict.

Approved source surface:

- `src/cartesian_residual_gaussians/residual_basis.jl`;
- `src/cartesian_residual_gaussians/augmented_operators.jl`;
- `src/cartesian_gaussian_raw_blocks/mixed_hartree_blocks.jl`;
- `src/cartesian_reference_density/atomic_hf_reference_packets.jl`;
- `src/cartesian_reference_density/screened_hartree_correction.jl`;
- `src/cartesian_protected_ladder_bundle.jl` for narrow internal composition
  only, with no public recipe or artifact change.

The neutral raw-block owner may add one internal explicit spherical Gaussian
potential entry point that reuses existing factor and mixed/self block kernels
for fitted-potential `GG/GA/AA`. Packet or ladder files must not duplicate
analytic Gaussian loops. Existing exact density-fit mixed-Hartree blocks remain
the small-system oracle.

No new source file, public export, artifact/schema field, or persistent
workflow object is approved. Related diagnostics must be nested in compact
records rather than copied as a flat stage field cloud. Target source budget is
at most `350` added lines across the approved files, with duplicate compact
selection deletion reported separately.

For this embedding-equivalence amendment only, source edits are limited to
`src/cartesian_reference_density/atomic_hf_reference_packets.jl` and the
existing private additive-reference caller if nested diagnostic forwarding is
directly required. The other implemented addref surfaces are not reopened.

For the moment-polish retirement only, the existing private additive caller
may reject retired packets and return the total/self/cross consistency
diagnostics. Packet fit deletion and correction acceptance behavior remain in
their active module owners. No geometry, `H1_L`, `Vee_L`, placement, or
ordinary no-reference behavior change is approved.

Approved committed test surface:

- `test/misc/runtests.jl` for a tiny mandatory-union/additive-density contract;
- `test/nested/cartesian_screened_hartree_correction_runtests.jl` for additive
  block validation and anchor algebra.

Do not add a new committed test file or binary fixture.

The existing nested test may also cover the owner-local embedding distinction:
exact packet self-integrity and structural mapping remain hard failures, while
a numerically equivalent translated/reconstructed overlap block may have a
different mapped fingerprint. Differences above `1e-10`, corrupt packet
fingerprints, reordered labels/powers, wrong owners, and wrong centers must
fail. Cr2 preflight follows only after focused source tests pass.

First end-to-end acceptance is an ignored, source-backed, physically padded
Be2 construction with two converged Be core `2e` cc-pV5Z, `lmax = 1` packets.
Use driver-style padding of at least `10` bohr and inspect terminal due
diligence. One ignored `tmp/work/*.jl` probe and durable text/TSV output under
`/Users/srw/dmrgtmp` are allowed; generated packets, Hamiltonians, and matrix
fixtures are not committed. Required evidence:

- one shared compact residual object for geometry and operators;
- omitted-reference protected-member parity with the current path;
- mandatory union Gram/rank and roundoff per-packet recovery through `L`;
- per-packet trace `2` and total `P0/q0` charge `4`;
- explicit `E_AA`, `E_BB`, both cross orderings, and
  `E0 = E_AA + E_BB + 2E_AB`;
- finite/symmetric placed raw blocks, `J0_F`, `J0_L`, and `Delta_J0`;
- strict derivative/algebra checks;
- total fitted-potential consistency error and matching self/cross
  decomposition, reported without a `1e-8 Ha` magnitude gate;
- current optional staged-selection diagnostics after mandatory inclusion;
- exact confirmation that unscreened `H1_L` and `Vee_L` were not mutated;
- parent bounds, axis counts, padding/radius, final dimension, retained counts,
  shell/slab topology, warning flags, and phase timings.

No endpoint energy or SCF assertion is required. The earlier polish-assisted
padded Be2 energy result is historical false-start evidence; regenerate
ordinary Be/Ne/Cr packets and rerun the bounded construction before further
consumption. The additive construction remains implemented, but this does not
authorize a production claim or repo Cr2 test.

Explicit exclusions:

- protected one-center atom compactness or removal of the current two-owner
  compactness assumption;
- counterpoise artifacts or retention of separated kinetic/unit-nuclear
  matrices;
- compact/high transfer helpers;
- corrected protected-localized artifact variants;
- public driver/API/export/default changes;
- solver, HF, MP2-NO, or production Cr2 workflow;
- `Vee` transformation, `C' V C`, or four-index interactions;
- EGOI, exchange, rho0 row-gauge, residual-selection, injection-threshold, or
  mapping-default changes;
- endpoint or publication claims.

Before any future compact/high transfer authority, require a same-commit audit
of `X`, `S_AA`, compact residual geometry, protected `G_L/A_L`, native
ordering, and exact final-basis cross overlap. A historical final-row mismatch
does not establish a Coulomb-accuracy basis change.

Decision rule: if exact packet self-integrity, exact structural owner-local
mapping, numerical mapped-overlap equivalence, mandatory recovery, neutral
fitted-potential blocks, compact cross energy, and existing correction reuse
cannot all fit this in-memory surface, stop and report the exact missing native
fact. Do not solve the blocker with public/artifact/solver wiring or a second
screened-Hartree formula.

### HP-RG-PROTECT-INJECT-FN-01 / HP-RG-PROTECT-INJECT-TEST-01 — staged protected-original geometry

Status: implemented internal, default-off facility.

Owner and canonical contract:

- `CartesianResidualGaussians` geometry owner;
- [Protected-localized basis convention](protected_localized_basis.md).

Approved source surface:

- `src/cartesian_residual_gaussians/residual_basis.jl`.

Approved validation/evidence surfaces:

- `tmp/work/cr2_source_backed_staged_protected_geometry_probe.jl`;
- `docs/src/developer/reports/cr2_staged_subspace_filter_870498b54/README.md`;
- no dedicated committed unit-test file.

Dependencies: the already-built ordered compact-first residual object,
`X_GA`, `S_AA`, and supplement owner/label/center metadata.

Permission summary: build protected and broad original subspaces, apply
separate Gaussian Gram, representability, and fake-RDM subspace gates, and
return transform-ready `Z`, `B`, `Q_perp`, `F`, and geometry diagnostics.
Rejected broad directions never become MWG residual channels.

Non-goals: public wiring/defaults, artifacts, exact operator or interaction
construction under these IDs, solver/HF, selection-policy changes, or Cr2
production claims.

### HP-RG-PROTECT-ONEBODY-AUDIT-01 — protected fixed-sector one-body audit

Status: completed historical measurement; not active source authority.

Evidence: manager running-log Passes 254-255 and
`docs/src/developer/reports/cr2_protected_onebody_audit_eaf05a38c/README.md`.
The audit established the dataflow later implemented under the source/test IDs
below.

### HP-RG-PROTECT-ONEBODY-FN-01 / HP-RG-PROTECT-ONEBODY-TEST-01 — exact protected one-body transform

Status: implemented internal, default-off facility.

Owner and canonical contract:

- `CartesianResidualGaussians` operator owner;
- [Protected-localized basis convention](protected_localized_basis.md).

Approved source surfaces:

- `src/cartesian_residual_gaussians/augmented_operators.jl`;
- `src/cartesian_residual_gaussians/residual_basis.jl` only for
  transform-ready protected geometry.

Approved validation/evidence surfaces:

- `tmp/work/cr2_protected_onebody_dense_source_replay.jl`;
- `docs/src/developer/reports/cr2_protected_onebody_audit_eaf05a38c/README.md`;
- no dedicated committed unit-test file.

Dependencies: `HP-RG-PROTECT-INJECT-FN-01` geometry and exact one-body
`GG/GA/AA` blocks.

Permission summary: construct dense exact fixed-sector kinetic, per-center
unit nuclear, and assembled `H1_F` matrices with orthogonality, symmetry, and
low-spectrum diagnostics.

Non-goals: matrix-action frameworks, public wiring/defaults, artifacts,
interaction rotation, solver/HF, residual-policy changes, screened-reference
work, or Cr2 production claims.

### HP-RG-PROTECT-VEE-AUDIT-01 — protected interaction decision audit

Status: completed historical measurement; not source authority.

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

Test/evidence paths: `test/ida/cartesian_ida_hamiltonian_runtests.jl` and
`docs/src/developer/pqs_manager_running_log.md` Pass 299; implementation
evidence is commit `fd105b751`.

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

Test/evidence paths: `test/ida/cartesian_ida_hamiltonian_runtests.jl` and
`docs/src/developer/pqs_manager_running_log.md` Pass 299. No dedicated
committed protected-artifact test file exists; commit `fd105b751` records the
accepted bounded write/readback and rejection smokes.

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

Test/evidence paths: `test/ida/cartesian_ida_hamiltonian_runtests.jl` and
`docs/src/developer/pqs_manager_running_log.md` Pass 301; implementation
evidence is commit `3fe2af697`.

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

Test/evidence paths: `test/ida/cartesian_ida_hamiltonian_runtests.jl` and
`docs/src/developer/pqs_manager_running_log.md` Pass 301. No dedicated
committed row-locality test file exists; commit `3fe2af697` records the
accepted center, inverse-permutation, sector, spread, and legacy-no-locality
smokes.

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

Owner and evidence:
[rho0 and reference-density correction history](rho0_reference_density_matrix.md).
No source or test authority remains.

### HP-RHO0-REFDENS-AUDIT-01 - fixed-P0 audit

Status: completed and superseded historical measurement evidence.

Owner and evidence:
[rho0 and reference-density correction history](rho0_reference_density_matrix.md).
Candidate source IDs `HP-RHO0-REFDENS-FN-01` and
`HP-RHO0-REFDENS-ERI-01` remain unapproved.

### HP-RHO0-REFDENS-MIXH-AUDIT-01 - exact mixed-Hartree seam audit

Status: completed historical measurement evidence. Its durable result is
the implemented neutral MIXH/FEXACT family below.

### HP-RHO0-MIXH-GG-FN-01 - exact mixed-Hartree GG block

Status: implemented.

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

Committed consumer test:
`test/nested/cartesian_screened_hartree_correction_runtests.jl`.

Compact gate: package load; bounded H/Be/Be2 finite/symmetric output;
angular and off-diagonal same-center pair coverage; dense-oracle spots; no
Cr/Cr2. Historical acceptance is manager Pass 284 and commit `efaee93f6`.

### HP-RHO0-MIXH-GAAA-FN-01 - exact mixed-Hartree GA/AA blocks

Status: implemented.

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

Test ownership: no dedicated committed test file. The accepted compact
ignored-probe gate covered prior `GG` parity, finite `GA`, symmetric `AA`,
angular cases, and dense-oracle spots. Historical acceptance is manager
Pass 286 and commit `daac231d0`.

### HP-RHO0-MIXH-FEXACT-FN-01 - protected exact-Hartree transform

Status: implemented.

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

Test ownership: no dedicated committed test file. The accepted compact
ignored-probe gate covered H/Be/Be2 raw-block replay, fixed/localized finite
symmetry, dimension and protected-geometry checks, and dense-oracle spots.
Historical acceptance is manager Pass 288 and commit `40a6f7e99`.

### HP-RHO0-FAPP-AUDIT-01 - approximate fixed-P0 Fock audit

Status: completed historical measurement evidence.

### HP-RHO0-FAPP-FN-01 / HP-RHO0-FAPP-TEST-01 - approximate IDA energy/Fock seam

Status: implemented, caller-free, and dormant retirement candidates.

Owner and source:
`src/cartesian_ida_hamiltonian.jl`.

Canonical lifecycle:
[rho0 and reference-density correction history](rho0_reference_density_matrix.md).
There is no committed live caller or dedicated committed test. These IDs
authorize no new caller, source work, correction policy, or public surface.

### HP-RHO0-ANCHOR-FN-01 / HP-RHO0-ANCHOR-TEST-01 - old full-interaction anchor

Status: superseded historical evidence with no authority.

The old `Delta_F0_alpha/beta` interpretation is not a Hartree correction
contract and must not be revived.

### HP-RHO0-CORR-AUDIT-01 - corrected-Hamiltonian audit

Status: completed and superseded historical measurement evidence.

### HP-RHO0-JANCHOR-FN-01 / HP-RHO0-JANCHOR-TEST-01 - direct-Hartree anchor

Status: source-backed but superseded in use; dormant retirement candidates.

Source:

- `src/cartesian_ida_hamiltonian.jl`;
- `src/cartesian_residual_gaussians/augmented_operators.jl`.

Current replacement:
`src/cartesian_reference_density/screened_hartree_correction.jl`, governed
by [screened Hartree correction assembly](screened_hartree_correction_assembly.md)
and [screened Hartree residual density](screened_hartree_residual_density.md).

These IDs authorize no new caller, source work, artifact, public workflow,
solver, exchange, or Cr/Cr2 work.

### HP-RHO0-XPAIR-AUDIT-01 - exchange/direct pairing question

Status: approved but deferred measurement-only question.

Owner and evidence:
[rho0 and reference-density correction history](rho0_reference_density_matrix.md).
It is not a current blocker or source lane and authorizes no tracked source,
test, artifact, public workflow, solver, exchange implementation, or Cr/Cr2
work.

Candidate future IDs, not approved:

- `HP-RHO0-REFDENS-FN-01`;
- `HP-RHO0-REFDENS-ERI-01`.

## Approved For Cartesian Gaussian Raw-Block Nuclear Owner

This section approves only the neutral uncharged by-center nuclear raw-block
slice recorded in `cartesian_gaussian_raw_blocks_nuclear.md`. It does not
approve a broad Gaussian raw-block framework.

### HP-CGRB-FILE-01 — neutral Cartesian Gaussian raw-block module files

Approved internal module and files:

```text
src/cartesian_gaussian_raw_blocks/CartesianGaussianRawBlocks.jl
src/cartesian_gaussian_raw_blocks/nuclear_blocks.jl
```

`src/GaussletBases.jl` may add only the internal include needed to load the
module, with include-order changes limited to immediate dependency/caller
needs. No public export is approved.

### HP-CGRB-FN-01 — exact uncharged Gaussian nuclear raw blocks

Approved internal kernel family:

```text
cartesian Gaussian nuclear raw blocks by center
```

The exact Julia names may follow local module style, but they must be neutral:
no `r3`, `residual`, `qw`, route-stage, report, or status vocabulary.

Approved numerical outputs:

- parent-supplement `G-A` matrices by nuclear center;
- supplement-supplement `A-A` matrices by nuclear center;
- uncharged unit attraction convention, `U_A = -1/r_A`.

Approved construction details:

- analytic one-dimensional nuclear factor construction;
- reuse across unique nuclear coordinates;
- upper-triangular `A-A` assembly and mirroring;
- function-local scratch/workspace reuse;
- term-first contraction over the Gaussian expansion.

The kernel must not apply physical nuclear charges, perform terminal
projection, transform into residual bases, create overlap/kinetic/moment
blocks, assemble Hamiltonians, or create persistent caches/bundles.

### HP-CGRB-FN-02 — nuclear one-dimensional axis-family reuse

Approved source owner:

```text
src/cartesian_gaussian_raw_blocks/nuclear_blocks.jl
```

Approved optimization target: reorganize exact uncharged nuclear raw-block
construction around unique supplement one-dimensional axis families rather
than flattened 3D orbital labels.

Approved concepts:

- unique supplement axis-family inventory independent of flattened 3D orbital
  labels;
- integer map `orbital_axis_family[orbital, axis] -> family_id`;
- unique `G-A` table keys
  `(axis, supplement_axis_family, nuclear_axis_coordinate)`;
- unique `A-A` table keys
  `(axis, canonical(left_family, right_family), nuclear_axis_coordinate)`;
- transpose/orientation flags when canonical `A-A` family order is reversed;
- term-first filling of each required one-dimensional table at most once per
  Coulomb Gaussian term;
- reuse of those tables across all 3D orbitals and orbital pairs that reference
  the same axis families;
- coupled primitive-pair contraction
  `sum_pq c_p c_q Ix[p,q] Iy[p,q] Iz[p,q]`;
- function-local workspaces and integer lookup plans only.

Independent contraction of x/y/z axis tables into separate scalar contractions
is forbidden. The kernel must not introduce persistent caches, metadata,
status/report fields, route objects, payload structs, public API/export,
artifact changes, Residual Gaussian algorithm changes, Qiu-White route
semantic changes, overlap/kinetic/moment migration, Cr2 facade support, or Cr2
artifact workflow.

### HP-CGAI-FN-01 — optional Cartesian Gaussian axis helper

Approved source owner:

```text
src/cartesian_gaussian_axis_integrals.jl
```

Status: superseded as a performance endpoint. It remains an optional helper
surface only if needed by `HP-CGRB-FN-02`.

Optional internal helper concept:

```julia
_cartesian_gaussian_axis_integral_table!(
    destination,
    left_exponents,
    left_centers,
    left_powers,
    left_prefactors,
    right_exponents,
    right_centers,
    right_powers,
    right_prefactors,
    term;
    factor_exponent = 0.0,
    factor_center = 0.0,
)
```

The helper fills an already allocated destination matrix and must return the
same values as `_cartesian_gaussian_axis_integral_table(...)` without
allocating the result matrix. The existing scalar
`_cartesian_gaussian_axis_integral(...)` behavior is unchanged. The allocating
helper may delegate to the in-place helper if that preserves behavior cleanly.
The same owner may also add a specialized nonallocating nuclear-factor scalar
integral for the `:factor` term if needed by the `HP-CGRB-FN-02` family-reuse
kernel.

Allowed consumer surface:

```text
src/cartesian_gaussian_raw_blocks/nuclear_blocks.jl
```

That file may consume the optional helper only inside the exact uncharged
by-center Gaussian nuclear raw-block construction under `HP-CGRB-FN-02`. Do
not treat result-matrix allocation removal as the accepted Cr2 optimization
target. No public API, export, persistent cache, raw-block payload,
metadata/status/report field, artifact change, route object, Residual Gaussian
algorithm change, Qiu-White route semantic change, overlap/kinetic/moment
migration, Cr2 facade, or Cr2 artifact workflow is approved by this ID.

### HP-CGRB-WIRE-01 — Residual Gaussian and Qiu-White rewiring

Approved caller rewiring surfaces:

```text
src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl
src/ordinary_qw_raw_blocks.jl
src/ordinary_qw_operator_assembly.jl
```

The implementation sequence is binding:

1. Extract current nuclear `G-A`/`A-A` behavior preserving conventions.
2. Rewire Residual Gaussian and Qiu-White callers to the neutral kernel.
3. Delete duplicate route-local nuclear loops once parity is established.
4. Optimize allocation inside the neutral owner only after extraction parity.

No Qiu-White route objects, Residual Gaussian selection logic, augmented
operator transforms, terminal projection, parent construction, artifact
workflow, report/status/payload fields, public API, or Cr2 facade/artifact
workflow may be added under this ID.

### HP-CGRB-TEST-01 — nuclear extraction validation

Approved validation:

- existing H2 Residual Gaussian endpoint unchanged;
- Be2 Residual Gaussian endpoint unchanged, as ignored validation if needed;
- Cr2 q4 exact nuclear `G-A`/`A-A` blocks match the current implementation at
  roundoff, as ignored measurement only;
- one small Qiu-White nuclear parity fixture.

Approved committed standalone parity file, if needed:

```text
test/nested/cartesian_gaussian_raw_blocks_nuclear_runtests.jl
```

This file must not be added to `test/runtests.jl` without a later amendment.
It must not validate Cr2 workflow, artifacts, public API, route statuses,
report mirrors, payload fields, or Residual Gaussian internals.

## Approved For Cartesian Gaussian Raw-Block Non-Nuclear Owner

This section approves only the non-nuclear raw-block slice recorded in
`cartesian_gaussian_raw_blocks_non_nuclear.md`. It extends the existing neutral
Cartesian Gaussian raw-block owner under new `HP-CGRB-NN-*` IDs. It must not be
implemented under `HP-CGRB-FN-02`.

### HP-CGRB-NN-FILE-01 — non-nuclear raw-block file

Approved owner file:

```text
src/cartesian_gaussian_raw_blocks/non_nuclear_blocks.jl
```

Allowed module plumbing:

```text
src/cartesian_gaussian_raw_blocks/CartesianGaussianRawBlocks.jl
```

Only the include needed to load `non_nuclear_blocks.jl` is approved there.
Root include changes in `src/GaussletBases.jl` are not expected and are not
approved unless a later amendment identifies a real include-order blocker.

### HP-CGRB-NN-FN-01 — exact non-nuclear Gaussian raw blocks

Approved internal kernel family:

```text
cartesian Gaussian non-nuclear raw blocks
```

The exact Julia names may follow local module style, but they must be neutral:
no `r3`, `residual`, `qw`, route-stage, report, or status vocabulary.

Approved numerical outputs:

- parent-supplement `G-A` overlap, kinetic, coordinate moments, and second
  moments;
- supplement-supplement `A-A` overlap, kinetic, coordinate moments, and second
  moments.

Approved construction details:

- analytic one-dimensional table construction;
- unique supplement axis-family reuse;
- canonical `A-A` family-pair table keys and orientation handling;
- upper-triangular `A-A` assembly and mirroring;
- function-local scratch/workspace reuse;
- coupled product-axis contraction preserving existing Qiu-White values;
- reuse of once-built overlap `G-A` for residual setup mixed overlap
  `X = G' S A` and exact augmented-operator assembly when both are built in
  the same local construction call.

The kernel may return a compact fixed-field internal result containing only the
approved raw matrices. It must not be a status object, route stage, report
payload, metadata carrier, persistent cache, broad provider bundle, or artifact
data.

### HP-CGRB-NN-WIRE-01 — Residual Gaussian and Qiu-White non-nuclear rewiring

Approved caller rewiring surfaces:

```text
src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl
src/ordinary_qw_raw_blocks.jl
src/ordinary_qw_operator_assembly.jl
```

The implementation sequence is binding:

1. Extract current Qiu-White non-nuclear `G-A`/`A-A` behavior preserving
   conventions.
2. Rewire Residual Gaussian exact-operator construction and residual mixed
   overlap setup to consume the neutral output.
3. Rewire Qiu-White consumers to the neutral output.
4. Delete duplicate route-local non-nuclear loops once parity is established.
5. Optimize allocation inside the neutral owner only after extraction parity.

No nuclear raw-block changes, final-basis `G-G` product-matrix optimization,
terminal projection, Residual Gaussian algorithm changes, augmented-operator
transform changes, Qiu-White semantic changes, Qiu-White route objects, parent
construction, persistent cache, report/status/payload fields, public API,
artifact workflow, or Cr2 facade/artifact workflow may be added under this ID.

Post-extraction clarification: the accepted diatomic Qiu-White non-nuclear
path consumes the neutral owner. Remaining QW-local non-nuclear cross/self
helpers used by `_qwrg_atomic_cartesian_blocks_3d`, atomic QW reference
operators, factor-term output, hybrid sidecars, dense-parent probes, and
CPB/provider surfaces are retained live reference/sidecar/provider surfaces.
They are not dead duplicates approved for deletion under this ID. Future
rewiring must first name the affected sidecar/provider surfaces and, if needed,
approve neutral ownership for factor blocks such as `factor_ga`/`factor_aa`.

### HP-CGRB-NN-TEST-01 — non-nuclear extraction validation

Approved validation:

- existing H2 Residual Gaussian endpoint unchanged;
- Be2 Residual Gaussian endpoint unchanged, as ignored validation if needed;
- Cr2 q4 non-nuclear `G-A`/`A-A` overlap, kinetic, coordinate moment, and
  second-moment blocks match the current implementation at roundoff, as
  ignored measurement only;
- residual setup mixed overlap `X` matches the current construction at
  roundoff;
- one small Qiu-White non-nuclear parity fixture.

Approved committed standalone parity file, if needed:

```text
test/nested/cartesian_gaussian_raw_blocks_non_nuclear_runtests.jl
```

This file must not be added to `test/runtests.jl` without a later amendment.
It must not validate Cr2 workflow, artifacts, public API, route statuses,
report mirrors, payload fields, or Residual Gaussian internals.

## Approved For R3 Terminal G-G Product Matrices

This section approves only the terminal final-basis `G-G` product-matrix
optimization recorded in `r3_terminal_gg_product_matrices.md`. It is owned by
`CartesianFinalBasisRealization`, not by `CartesianGaussianRawBlocks`.

### HP-R3GG-FN-01 — R3/RG terminal G-G product-matrix optimization

Approved owner:

```text
Owner module: CartesianFinalBasisRealization
```

Approved source files:

```text
src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl
src/cartesian_final_basis_realization/pqs_terminal_one_body.jl
```

The first implementation should prefer editing only
`pqs_terminal_residual_gto.jl`. Edits to `pqs_terminal_one_body.jl` are
approved only for a small internal terminal-product workspace or multi-product
helper needed to reuse function-local buffers across consecutive product
assemblies.

Approved product matrices:

- kinetic `K_GG`;
- coordinate moments `x_GG`, `y_GG`, `z_GG`;
- second moments `x2_GG`, `y2_GG`, `z2_GG`.

Approved implementation shapes:

- accumulate the three kinetic-axis product contributions into one destination;
- reuse an already constructed base Hamiltonian kinetic `G-G` block in the
  same-construction path when available and validated equal;
- build coordinate and second-moment `G-G` products one axis at a time and
  transform immediately;
- share function-local scratch/workspace across consecutive terminal product
  assemblies;
- delete or simplify `_r3a_product_matrix(...)` when replaced and no live
  caller remains.

This ID does not approve `G-A`/`A-A` raw-block changes, nuclear raw-block
changes, unit-nuclear `U_A` Gaussian-sum changes, terminal basis realization
changes, residual Gaussian algorithm changes, Qiu-White semantic changes,
IDA/MWG changes, parent construction, persistent caches, metadata,
report/status/payload fields, public API/export, artifact changes, Cr2 facade
support, or Cr2 artifact workflow.

Line budget: at most `100` added `src` lines total. If implementation needs a
broad product-operator framework, persistent workspace/cache object, files
outside the approved source files, or a public/internal payload, stop and
request a new docs-only amendment.

### HP-R3GG-TEST-01 — terminal G-G product validation

Approved validation:

- existing H2 Residual Gaussian endpoint unchanged;
- Be2 Residual Gaussian usability/performance measurement unchanged except for
  allowed timing/allocation improvement;
- Cr2 q4 `K_GG`, coordinate moment `G-G`, and second-moment `G-G` products
  match the current construction at roundoff as ignored validation;
- augmented exact operators remain finite and symmetric;
- base `G-G` block equality checks in the existing H2 endpoint still pass;
- Cr2 q4 exact-operator allocation is remeasured after parity.

No new committed test file is approved by this ID.

## Approved For R3 Unit-Nuclear U_GG Gaussian Sum

This section approves only the terminal final-basis unit-nuclear `U_GG`
Gaussian-sum optimization recorded in
`r3_unit_nuclear_ugg_gaussian_sum.md`. It is owned by
`CartesianFinalBasisRealization`, not by `CartesianGaussianRawBlocks`.

### HP-R3UN-FN-01 — R3/RG unit-nuclear U_GG Gaussian-sum optimization

Approved owner:

```text
Owner module: CartesianFinalBasisRealization
```

Approved source files:

```text
src/cartesian_final_basis_realization/pqs_terminal_one_body.jl
src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl
```

The first implementation should prefer `pqs_terminal_one_body.jl`. Edits to
`pqs_terminal_residual_gto.jl` are approved only for narrow R3 exact-operator
caller wiring needed to use function-local scratch or the optimized helper.

Approved target functions:

```text
_accumulate_terminal_gaussian_sum!
_terminal_gaussian_sum_action
```

Approved implementation shapes:

- reuse function-local scratch/workspace across Gaussian-sum terms and center
  calls;
- accumulate terminal Gaussian-sum contributions in-place into the caller's
  destination;
- reduce avoidable allocation in factor lookup and terminal Gaussian-sum action
  construction;
- add small internal scratch arguments or file-local helpers only if they remain
  inside `CartesianFinalBasisRealization` and create no persistent state;
- simplify or delete allocation-heavy helper code inside the targeted
  Gaussian-sum path after parity.

This ID does not approve neutral raw-block changes, terminal kinetic/moment
`G-G` product changes, residual Gaussian selection/orientation/transform
changes, MWG/IDA changes, Qiu-White semantic changes, route/stage setup
cleanup, raw-block setup cleanup, parent construction, terminal basis
realization changes, persistent caches/workspaces, broad Gaussian-sum
frameworks, metadata/report/status/payload fields, artifacts, public
API/export, Cr2 facade support, or Cr2 artifact workflow.

Line budget: at most `100` added `src` lines total. If implementation needs a
broad Gaussian-sum framework, persistent cache/workspace object, files outside
the approved source files, or source edits outside the terminal unit-nuclear
`U_GG` path, stop and request a new docs-only amendment.

### HP-R3UN-TEST-01 — unit-nuclear U_GG validation

Approved validation:

- existing H2 Residual Gaussian endpoint unchanged;
- Be2 Residual Gaussian facade/readback unchanged except for allowed
  timing/allocation improvement;
- Cr2 q4 exact-operator audit reports before/after unit-nuclear `U_GG`
  allocation and total wrapper allocation;
- Cr2 q4 unit-nuclear `U_GG` block replay parity and final exact augmented
  operator parity at roundoff;
- exact operators remain finite and symmetric.

No new committed test file is approved by this ID.

## Approved For R3 Same-Construction Base K/U Reuse

This section approves only narrow reuse of already-built same-construction
base final-basis kinetic and unit-nuclear blocks in supplemented residual-GTO
/ MWG exact augmented operators. It is an orchestration reuse lane, not a
terminal product, Gaussian-sum, raw-block, residual-basis, or interaction
algorithm lane.

Evidence recorded before approval: a replay that reused same-construction base
`K_GG` and unit `U_GG[A]` blocks had exact operator delta `0.0` and reduced the
exact augmented-operator replay to `0.8620s / 1237.136 MiB`.

### HP-R3BASE-FN-01 — same-construction base K/U reuse

Approved owner:

```text
Owner module: CartesianFinalBasisRealization plus narrow caller wiring
```

Approved source files:

```text
src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl
src/cartesian_base_hamiltonian.jl
```

Approved behavior:

- `pqs_terminal_residual_gto_augmented_products(...)`, or its approved caller
  wrapper, may accept a trusted same-construction base kinetic matrix and use
  it as the `G-G` kinetic block for `transform_augmented_operator`;
- `pqs_terminal_residual_gto_augmented_unit_nuclear(...)`, or its approved
  caller wrapper, may accept trusted same-construction unit nuclear
  `U_GG[A]` matrices and use them as the `G-G` unit blocks for
  `transform_augmented_operator`;
- `cartesian_residual_gto_mwg_hamiltonian(...)` and staged helpers in
  `src/cartesian_base_hamiltonian.jl` may pass `base_ham.kinetic` and
  `base_ham.nuclear_attraction_unit_by_center` into the augmented operator
  construction;
- current behavior must be preserved when trusted base blocks are not supplied.

Trust condition:

- the base Hamiltonian, terminal basis realization, parent bundles, residual
  basis, and supplement must come from the same
  `cartesian_base_working_basis(...)` construction path;
- implementation must validate matrix dimensions and center count before
  reuse;
- no provenance payload, metadata proof, report field, status object, or
  persistent cache is required or approved for this trust check.

This ID does not approve public API/export changes, canonical-driver changes,
raw-block changes, residual selection/orientation/transform changes, MWG/IDA
convention changes, terminal product or Gaussian-sum kernel rewrites,
persistent cache/workspace objects, metadata/status/report/artifact schema
fields, route/stage setup cleanup, committed tests, Cr2 workflow, or source
files outside the two approved files.

Line budget: target under `100` added `src` lines. If trusted
same-construction provenance cannot be guaranteed by local call shape plus
dimension/center validation, or if implementation needs public payloads,
metadata, or stage objects, make no source commit and report the blocker.

### HP-R3BASE-TEST-01 — base K/U reuse validation

Approved validation:

- `git diff --check`;
- package load;
- existing H2 R3 endpoint unchanged;
- Be2 supplemented facade/readback unchanged except allowed
  timing/allocation improvement;
- Cr2 exact-operator attribution audit or focused ignored replay showing base
  `K_GG` / unit `U_GG[A]` reuse parity and allocation effect;
- final exact operators finite and symmetric.

No new committed test file, Cr2 artifact, Cr2 workflow, public API/export,
driver workflow change, metadata/status/report field, or artifact schema
change is approved by this ID.

### HP-R3BASE-DRV-WIRE-01 — canonical driver K/U reuse call-site wiring

Approved source file:

```text
bin/cartesian_ham_builder.jl
```

Approved behavior:

- in supplemented mode only, pass `base_ham.kinetic` into
  `cartesian_residual_gto_augmented_products(...)` as `base_kinetic`;
- in supplemented mode only, pass
  `base_ham.nuclear_attraction_unit_by_center` into
  `cartesian_residual_gto_augmented_unit_nuclear(...)` as
  `base_unit_nuclear`;
- keep the current public inputs, hooks, timing labels, visible stage sequence,
  artifact schema, and driver contract unchanged.

This ID is only call-site wiring so the canonical driver uses the already
approved same-construction base K/U reuse path. It does not approve source or
kernel changes, diagnostics, new hooks, new timing labels, public API/export
changes, artifact changes, tests/fixtures, Cr2 workflow, or edits outside
`bin/cartesian_ham_builder.jl`.

Failure rule: if the driver call-site update needs any visible driver contract
change, make no source commit and report the blocker.

### HP-R3BASE-DRV-TEST-01 — driver K/U reuse validation

Approved validation:

- `git diff --check`;
- package load;
- H2 supplemented driver artifact/readback;
- Be2 supplemented driver artifact/readback if practical;
- no Cr2 run.

No new committed test file, fixture, diagnostic, hook, timing label, public
input, artifact schema, or Cr2 workflow is approved by this ID.

## Approved For Canonical Cartesian Driver Usability

This section approves only the compact artifact-producing canonical driver
workflow recorded in `cartesian_driver_usability_workflow.md`. It is workflow
authority over approved producer surfaces, not algorithm, kernel, solver,
artifact-schema, or diagnostic authority.

### HP-DRV-FILE-01 — canonical driver file

Approved file:

```text
bin/cartesian_ham_builder.jl
```

No other `bin`, `tools`, `src`, `test`, or committed driver-input fixture file
is approved by this ID.

### HP-DRV-FN-01 — compact functional driver workflow

Approved invocation shape:

```text
julia --project=. bin/cartesian_ham_builder.jl [input.jl] [key=value ...]
```

Approved behavior:

- visible editable defaults near the top of the driver;
- optional trusted local Julia input file for project-specific defaults;
- later command-line `key=value` overrides;
- visible public `system`, `basis`, and optional `supplement` contract
  construction before calling an approved facade;
- compact normalized run summary;
- coarse user-facing phase timing;
- visible physics-level construction stages through the staged producer surface;
- base or supported supplemented Hamiltonian construction through approved
  producer surfaces;
- artifact write;
- optional readback check.

Approved configuration concepts are `basisname`, `system`, base `basis`,
optional `supplement`, `nesting`, `hamfile`, `padding`, `check_file`,
`print_contract`, `print_timing`, and `expected_dimension`.

Compact summary printing and artifact readback checks remain allowed workflow
behavior, but they are not open-ended hooks and must not introduce route,
diagnostic, artifact-schema, or solver controls.

`basisname = nothing` selects base mode. `basisname !== nothing` selects a
supported supplemented mode and is the visible supplement basis label. The
original driver workflow lane covered supplemented diatomics only;
`HP-COMP-SUPPATOM-*` separately approves relaxing the old `Natom == 1`
rejection.

`padding` is a public physical box-padding control. For one-center atoms it
maps to the base facade `radius`. For z-axis diatomics it maps to the existing
facade extents as padding around the two nuclei; under the current origin-based
z-axis contract this means `xmax_parallel = max(abs(z_i)) + padding` and
`xmax_transverse = padding`.

Approved hooks are only `check_file`, `print_contract`, `print_timing`, and
`expected_dimension`. They may support human expert review and
Codex-controlled artifact checks. They must not expose route internals,
stop-after stages, raw-block switches, allocation probes, artifact schema
dumps, solver controls, Cr2-specific workflow, or private helper calls.
`check_file` may contain compact public contract facts, artifact path, final
dimension, expected-dimension result, readback deltas, and coarse timing only.

This ID does not approve private route-stage controls, stop-after internals,
ladder probes, stage markers, fixture hacks, diagnostic knobs, underscored
package helper calls from the driver, raw-block provider switches,
report/status/payload dumps, metadata field clouds, allocation probes,
benchmark harness behavior, solver/RHF/ECP/EGOI/HamV6 workflow, public
API/export changes, artifact schema changes or dumps, committed test files,
committed driver-input fixtures, unsupported atom/supplement combinations, or
Cr2-specific workflow support. Generic explicit homonuclear z-axis Cr2 stress through
`HP-R3U-ZDI-WIRE-01` is separate ignored/user-run validation authority, not
driver-owned Cr2 support.

Line budget: at most `150` added `bin` lines. If implementation needs a parser
framework, source files outside the approved driver and staged producer
surfaces, committed input fixtures, route-stage diagnostics,
status/report/payload expansion, artifact schema changes, or Cr2-specific
workflow support, stop and request a new docs-only amendment.

### HP-DRV-NEST-FN-01 — construction-family driver input

Approved source file:

```text
bin/cartesian_ham_builder.jl
```

Approved behavior:

- add a visible public driver input `nesting`;
- accepted values are `:pqs` and `:wl`;
- default is `nesting = :pqs`;
- `nesting = :pqs` means the PQS source-box construction family;
- `nesting = :wl` means the White-Lindsey low-order construction family;
- include `nesting` in public contract construction, optional
  `print_contract`, and optional `check_file` output as a public contract fact.

This is a first-class construction-family choice, not a diagnostic route
switch. The driver must not expose internal route-family names, route
skeletons, retained-rule plans, raw-block switches, stop-after controls,
diagnostic knobs, old route-stage labels, private helper calls, allocation
probes, or route reports.

This ID does not approve public API/export changes, artifact schema changes,
stage-label changes, solver/ECP workflow, Cr2-specific behavior, broad driver
diagnostics, committed fixtures/tests, or source files outside the canonical
driver.

### HP-DRV-NEST-WIRE-01 — construction-family route mapping

Approved source files:

```text
bin/cartesian_ham_builder.jl
src/cartesian_base_hamiltonian.jl
```

Approved behavior:

- map public `nesting = :pqs` to the existing internal `:pqs_source_box`
  route family;
- map public `nesting = :wl` to the existing internal
  `:white_lindsey_low_order` route family;
- keep route skeletons, retained rules, raw-block switches, stop-after
  controls, diagnostics, and internal route-stage vocabulary hidden;
- preserve the existing public stage labels, Hamiltonian object, matrix keys,
  artifact schema, driver hooks, and solver-free workflow;
- reject unsupported combinations with clear `ArgumentError`s.

Supplemented `nesting = :wl` is governed by `HP-COMP-SUPPWL-*` for the
supported homonuclear z-axis diatomic composition cell. Unsupported geometry
or supplement combinations must still reject clearly rather than adding broad
White-Lindsey route behavior.

This ID does not approve new route algorithms, route-skeleton construction
changes, White-Lindsey materialization deletion, terminal lowering policy
changes, shellification behavior changes, numerical kernel changes, raw-block
changes, Residual Gaussian/MWG/IDA changes, artifact/provenance schema changes,
public API/export changes, committed tests, Cr2-specific workflow, or source
files outside the two approved files.

Line budget: at most `80` added source/bin lines, with net simplification
preferred where old hidden assumptions can be removed.

Failure rule: if `nesting = :wl` cannot produce a small base artifact/readback
through the existing White-Lindsey low-order route without broader route or
materialization work, make no source commit and report the exact blocker.

### HP-DRV-NEST-TEST-01 — construction-family validation

Approved validation:

- `git diff --check`;
- package load;
- current default `nesting = :pqs` base driver or facade artifact/readback;
- current default `nesting = :pqs` supplemented H2 driver/facade path if
  supplemented-mode input plumbing is touched;
- one small `nesting = :wl` base artifact/readback using a currently supported
  base geometry;
- explicit negative check or ignored smoke showing unsupported supplemented
  `nesting = :wl` combinations fail clearly outside the supported
  `HP-COMP-SUPPWL-*` cell;
- no Cr2 run.

No new committed test file, committed input fixture, artifact schema
validation, solver run, Cr2-specific driver run, or broad White-Lindsey
workflow validation is approved.

### HP-DRV-STAGE-FN-01 — visible physics-stage producer surface

Approved source file:

```text
src/cartesian_base_hamiltonian.jl
src/pqs_source_box_low_order_materialization.jl
src/cartesian_final_basis_realization/pqs_terminal_one_body.jl
src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl
```

Approved purpose: expose a small non-exported, non-underscored,
driver-facing staged producer surface so the canonical driver can execute and
time visible physics-level construction stages without calling private
underscored helpers.

`src/cartesian_base_hamiltonian.jl` remains the primary driver-facing owner.
The lower-level files listed above are approved only for behavior-preserving
operator-class stage factoring in their existing domains; they are not
approved for new algorithms, raw-block changes, or diagnostics.

Approved visible stages are:

- construct public `system`, `basis`, and optional `supplement`;
- build the base working basis / terminal realization;
- build base product/moment operators;
- build base unit-nuclear attraction operators;
- build base electron-electron / IDA interaction;
- assemble the base Hamiltonian;
- load or build the Gaussian supplement basis when `basisname !== nothing`;
- build residual Gaussian augmentation;
- build augmented product/moment operators;
- build augmented unit-nuclear attraction operators;
- build augmented electron-electron / residual-MWG interaction;
- assemble the supplemented Hamiltonian;
- write and check the artifact.

Approved physical operator classes are:

- product/moment operators: kinetic `K`, Cartesian coordinate moments
  `x`/`y`/`z`, and second moments `x^2`/`y^2`/`z^2` where present;
- unit-nuclear attraction: uncharged by-center `U_A` / `Vnuc` matrices before
  applying physical nuclear charges;
- electron-electron interaction: base localized IDA `Vee` and supplemented
  residual-MWG/IDA `Vee`.

The first and last stages remain driver/writer responsibilities. This ID
approves source factoring needed for the base working-basis/terminal
realization, base product/moment, base unit-nuclear, base IDA, Gaussian
supplement, residual augmentation, augmented product/moment, augmented
unit-nuclear, residual-MWG/IDA, and Hamiltonian assembly stages.

The staged surface may factor the existing `cartesian_base_hamiltonian(...)`
and `cartesian_residual_gto_mwg_hamiltonian(...)` bodies so that those facades
can remain wrappers over the same implementation. It may return existing
domain objects and small fixed-key ephemeral stage products required by the
next approved stage.

The staged surface must be a set of separate named construction-stage
functions. It must not be a single opaque replacement wrapper that hides the
same construction sequence under a new name. The canonical driver must be able
to bind visible local variables for the base realization, base products,
base unit-nuclear operators, base `Vee`, base Hamiltonian, supplement basis,
residual augmentation, augmented products, augmented unit-nuclear operators,
augmented `Vee`, and final Hamiltonian assembly.

This ID does not approve public exports, public API redesign, route-stage
objects, reports, status/result payloads, metadata field clouds, runtime-keyed
field groups, persistent caches, raw-block switches, allocation probes,
per-kernel timing frameworks, solver/ECP workflow, artifact schema changes, or
source files outside the four approved paths listed above.

Line budget: at most `200` added `src` lines across the approved source files.
If the staged surface requires a new module, source files outside the four
approved paths, a broad payload object, committed tests, new artifact keys,
raw-block changes, or kernel rewrites, stop and request a new docs-only
amendment.

### HP-DRV-STAGE-WIRE-01 — canonical driver staged wiring

Approved source file:

```text
bin/cartesian_ham_builder.jl
```

The canonical driver may call the `HP-DRV-STAGE-FN-01` staged producer surface
as separate top-level stage calls and assign local variables using the approved
physics-stage names. It may print or record coarse user-facing timings for
those stages. Replacing the current facade call with one all-in-one staged
wrapper call is not approved for the canonical driver.

Driver timing should expose the three physical operator classes separately:
product/moment, unit-nuclear, and electron-electron. These timings are
user-facing stage timings only. They must not become allocation probes,
raw-block timing controls, per-kernel instrumentation, or diagnostic stop
points.

This ID does not approve calls from the driver to underscored package helpers,
old route stages such as `cartesian_parent`, `cartesian_shells`,
`cartesian_units`, `cartesian_pair_terms`, or `cartesian_assembly`, raw-block
provider switches, stop-after controls, route diagnostics, allocation probes,
artifact schema dumps, solver controls, Cr2-specific workflow, or new
committed fixtures/tests.

### HP-DRV-STAGE-TEST-01 — staged driver validation

Approved validation:

- package load;
- H atom or H2 base driver run with visible base-stage timing/summary;
- H2 supplemented driver run with visible supplement/residual/operator/
  Hamiltonian stage timing/summary;
- artifact write/readback still passes for those runs;
- `expected_dimension`, `print_contract`, and `check_file` behavior still
  uses only public contract and coarse stage facts.

No committed test file, committed input fixture, Cr2-specific driver run,
solver run, or diagnostic harness is approved by this ID.

### HP-DRV-INV-FN-01 — canonical driver terminal-region inventory

Status: implemented.

Approved source files:

```text
bin/cartesian_ham_builder.jl
src/cartesian_base_hamiltonian.jl
```

Optional only if a compact accessor is directly required:

```text
src/pqs_source_box_route_driver_helpers.jl
src/cartesian_final_basis_realization/CartesianFinalBasisRealization.jl
src/cartesian_final_basis_realization/terminal_face_product_blocks.jl
```

Problem:

Large accidental identity sectors, such as the Cr2 z-end slab blowup, are
hard to notice from the current canonical driver output. Normal driver users
should see a compact basis-region inventory without running ignored probes and
without receiving a route-debug dump.

Approved behavior:

- print a compact terminal-region / shellification inventory as part of
  canonical driver output;
- include it for base construction;
- for supplemented construction, report at least the base terminal inventory
  plus final supplemented dimension;
- keep the output bounded and human-facing;
- preserve the existing driver stage sequence and public inputs;
- preserve artifact schema, matrix keys, reader behavior, and readback checks.

Minimum useful columns:

- region key/index or compact label;
- region kind;
- lowering kind or final realization kind;
- shell index or explicit unavailable status;
- support row count;
- retained/final column count;
- compression ratio;
- identity versus compact/product realization;
- index ranges for each axis, `x = i:j`, `y = k:l`, `z = m:n`;
- physical coordinate ranges for each axis, `x`, `y`, and `z`;
- slab normal axis, side, thickness, and stack index/count when applicable.

The summary should also print total base final dimension, supplemented final
dimension when applicable, and a clear count or visible rows showing any
direct identity slab sectors if they exist.

The geometry columns are part of the canonical inventory contract. They are
needed to catch shellification errors where every region is compact but the
z-axis slab stack is emitted only after the final shared shell instead of being
interleaved with the angular-balanced shell steps. Physical `x`/`y` ranges are
required, not only `z`, because the angular-balance rule compares the
transverse physical scale against the bond-axis margin.

This ID does not approve route skeleton exposure, source-mode inventories,
pair inventories, raw-block details, all-row listings, full metadata dumps,
recursive route-stage dumps, new driver inputs, flags, stop-after controls,
route switches, diagnostic switches, solver settings, broad status/report
payloads, artifact schema changes, reader changes, public API/export changes,
numerical construction changes, shellification changes, terminal lowering,
retained-unit changes, transform-contract changes, terminal-realization
changes, Residual Gaussian, MWG, IDA, Hamiltonian assembly, raw-block changes,
Cr2-specific workflow, committed Cr2 fixtures, or committed tests.

Line budget: target at most `80` added `src`/`bin` lines. This should be
formatting plus compact accessor work, not a reporting subsystem.

Failure rule: if the driver cannot print this from existing compact stage or
final-basis summaries without adding a broad payload, artifact fields, or a
route-report framework, make no source commit and report the missing summary
seam.

### HP-DRV-INV-TEST-01 — terminal-region inventory validation

Status: implemented validation coverage.

Approved validation:

- `git diff --check`;
- package load;
- bounded H2 or Be2 driver run showing the summary for `nesting = :pqs`;
- bounded H2 or Be2 driver run showing the summary for `nesting = :wl`;
- supplemented smoke if the printed summary touches supplemented-stage
  objects;
- artifact/readback deltas unchanged;
- output remains bounded and excludes source modes, pair inventories,
  raw-block details, all-row listings, and full metadata;
- output includes shell index or explicit unavailable status, index ranges for
  all axes, physical coordinate ranges for all axes, and slab stack facts when
  applicable;
- no Cr2 run required; optional user-side Cr2 run only.

No committed test file, committed driver-input fixture, Cr2-specific driver
run, artifact schema validation, solver run, or diagnostic harness is approved
by this ID.

### HP-DRV-SHELLDD-FN-01 — terminal shellification due-diligence report

Status: implemented. The canonical base working construction carries the
in-memory report and the canonical driver prints it.

Design note:

```text
docs/src/developer/designs/cartesian_hamiltonian_producer/terminal_shellification_due_diligence.md
```

Problem: compact terminal inventory rows can show that a region is compact
without showing whether the derived basis setup is what the user thinks it is
or whether the shell basis is physically adequate. The H2+
`ns = 5`/`ns = 7` audit showed a `complete_shell_1` with physical side lengths
`3.464 x 3.464 x 6.646`, source shape `(5, 5, 5)`, expected aspect-balanced
shape `(5, 5, 10)`, retained `98`, and expected retained scale `178`.
That was inadequate basis construction, and it should be visible in ordinary
producer/driver due diligence before interpreting energies or residual/
injection behavior.

Approved source files:

```text
src/cartesian_base_hamiltonian.jl
bin/cartesian_ham_builder.jl
```

Optional only if a compact accessor is directly required:

```text
src/pqs_source_box_route_driver_helpers.jl
```

Approved behavior:

- add one helper/report surface for terminal shellification due diligence;
- extend or wrap `_cartesian_terminal_inventory_rows(...)`;
- join existing terminal inventory rows with terminal retained-rule
  plan/support records;
- gather normalized system/geometry context and parent-axis summaries from
  existing staged producer objects;
- gather gausslet/IDA weight statistics only from existing weights already
  present in the construction path;
- produce an in-memory/report object first;
- expose the report from canonical driver/producer workflows;
- keep warning flags advisory by default;
- keep output bounded and row-oriented;
- preserve existing compact terminal-region inventory behavior unless it is
  intentionally extended by this report;
- preserve all numerical construction behavior.

Required report sections:

- normalized system and geometry context;
- parent axes, physical box, 1D center locations, and gausslet/IDA weight
  statistics;
- final-basis dimension and compression accounting;
- shell-by-shell terminal region table.

Required system/geometry fields:

- geometry kind, nesting, atom symbols, nuclear charges, `nup`, and `ndn`;
- validated atom locations;
- bond axis and bond length for z-axis diatomics;
- box-center convention and parent physical bounds/lengths;
- padding/radius-derived extents;
- snapped nuclear indices and physical snap errors;
- `core_spacing`, `reference_spacing`, and `tail_spacing` summary;
- parent axis counts.

Required parent-axis and weight fields:

- per-axis count, physical min/max/length, and bounded center preview or full
  per-axis center table when practical;
- min/median/max spacing and nearest spacing/index at each nucleus;
- core/tail region index spans when available;
- per-axis or aggregate gausslet/IDA weight count, sum, min/max, absolute sum,
  negative count, near-zero count/threshold, and large-weight warning where
  available.

These weight summaries are diagnostics only. They are not residual integral
weights, not MWG weights, and not automatic proof of quadrature quality.

Required dimension/accounting fields:

- parent grid size;
- direct/core, complete-shell, slab, compact-product, and identity columns;
- base final dimension;
- supplement, residual, and augmented dimensions when present;
- compression by class;
- large identity-sector count.

Required shell table fields:

- terminal order/key;
- role, region kind, and shell index;
- owner/contact/shared classification;
- index outer and inner boxes and shapes;
- physical bounds and physical side lengths for `x`, `y`, and `z`;
- physical aspect ratios;
- actual `source_mode_shape`;
- expected aspect-balanced `source_mode_shape` when applicable;
- source-mode count;
- retained count and final column range;
- lowering kind, retained rule, and realization rule/status;
- slab normal axis, side, thickness, stack index, and stack count when
  applicable;
- warning flags and warning summary.

Initial advisory warning flags should include rectangular physical shells
represented by cubic source modes, expected source shape larger than actual,
retained count below aspect-balanced scale, large identity sectors, missing
shell index, missing physical bounds, missing source-mode shape, slab rows
without native metadata, unavailable expected-shape diagnostics,
axis-center-table truncation, gausslet-weight anomalies, and padding/radius
values not reflected in the derived physical box when that is suspicious.

Contract:

- repo/driver workflows must expose this terminal due-diligence report for
  Cartesian/PQS terminal bases;
- consumers are expected to inspect it before interpreting energies, residual
  behavior, or injection behavior;
- warning flags are advisory diagnostics unless a caller/test/later policy
  explicitly enforces them.

Forbidden:

- artifact schema/provenance/reader changes;
- public input, public semantics, or driver contract changes;
- shellification policy changes;
- source-mode selection changes;
- aspect-balanced source-mode implementation;
- terminal lowering, retained-unit, terminal-realization, RG/MWG/IDA,
  Hamiltonian, raw-block, solver, or Cr2 workflow changes;
- route skeleton exposure, source-mode inventory dumps, pair inventories,
  raw-block details, all-row support listings, full metadata dumps, or
  recursive route-stage reports;
- dense coefficient, transform, pair, or raw support dumps;
- automatic failure on warning flags unless a later policy approves it.

Failure rule: if the due-diligence report cannot be built by extending/wrapping
`_cartesian_terminal_inventory_rows(...)` and compact accessors without adding
a broad report/payload framework, artifact fields, or shellification policy
changes, make no source commit and report the missing seam.

Line budget: target at most `180` added `src`/`bin` lines. This is one report
surface, not a new reporting subsystem.

Open follow-up: aspect-balanced complete-shell source modes are likely a
separate source-policy fix. Do not mix that with this reporting lane.

### HP-DRV-SHELLDD-TEST-01 — terminal shellification due-diligence validation

Status: implemented validation coverage.

Approved validation:

- `git diff --check`;
- package load if source is touched;
- bounded H2 or H2+ driver/producer smoke showing due-diligence report
  sections and shell rows;
- focused row inspection showing a rectangular physical shell warning when an
  existing bounded fixture has one;
- focused inspection of normalized system/geometry, axis/box summaries, and
  gausslet/IDA weight statistics;
- ordinary compact terminal inventory output remains bounded;
- artifact/readback matrix deltas unchanged if artifact writing is exercised;
- no Cr2 run required.

No committed fixtures or tests are approved by default.

### HP-PQS-ASPECTSHELL-FN-01 / HP-PQS-ASPECTSHELL-TEST-01 — PQS complete-shell aspect-aware source modes

Status: implemented internal construction policy and completed bounded
validation.

Owner: PQS terminal low-order route enrichment and multilayer shell
realization.

Canonical contract:
[PQS complete-shell aspect-aware source modes](pqs_complete_shell_aspect_source_modes.md).

Implemented source surfaces:

- `src/pqs_source_box_route_driver_helpers.jl`;
- `src/pqs_multilayer_shell_region_plan.jl`;
- `src/pqs_multilayer_shell_source_plan.jl`;
- `src/cartesian_base_hamiltonian.jl` for due-diligence shape reporting.

Permission summary: for bond-aligned z-axis diatomic PQS shared complete
shells, select angular-band source dimensions after shellification and before
lowering/retained/support records are frozen; carry one authoritative
`(q,q,L)` shape through multilayer realization and due diligence.

Validation/evidence: bounded H2/H2+ replay, noncubic retained-count checks,
finite/symmetric artifact/readback, and removal of the stale cubic-source
warning. Evidence is recorded in manager running-log Pass 247.

Non-goals: new source work under a completed lane, public inputs, WL policy,
shell ownership, thin-slab/direct-core changes, artifacts, RG/MWG/IDA,
injection, solver workflow, or Cr2 production claims.
### HP-DRV-TEST-01 — driver workflow validation

Approved validation:

- package load;
- public contract construction and optional `print_contract`/`check_file`
  output for at least one base run when driver construction code changes;
- H2 base driver run writes a `CartesianIDAHamiltonian` artifact and optional
  readback passes;
- H2 supplemented driver run writes a supplemented `CartesianIDAHamiltonian`
  artifact with approved compact `supplement_provenance/` and optional readback
  passes;
- optional ignored Be2 usability run if the implementation touches
  supplemented mode.

Validation input files, if needed, must be ignored `tmp/work` files. No
committed test file, committed driver-input fixture, Cr2-specific driver run,
or solver run is approved by this ID.

## Approved For Canonical Driver Atom Workflow

This section approves only the base atom workflow recorded in
`cartesian_driver_atom_workflow.md`. It is driver authority over the existing
base facade, not new atom physics, Residual Gaussian, artifact-schema, or
solver authority.

### HP-DRV-ATOM-FN-01 — explicit base atom driver workflow

Approved source file:

```text
bin/cartesian_ham_builder.jl
```

Approved behavior:

- accept explicit one-center atom input in `mode = :base`;
- require `atom_symbols`, `nuclear_charges`, `atom_locations`, `nup`, and
  `ndn`;
- require exactly one center at `(0.0, 0.0, 0.0)`;
- require finite positive explicit nuclear charge and neutral all-electron
  count `nup + ndn == round(Int, only(nuclear_charges))`;
- pass explicit one-center base `basis` fields, including `core_spacing`, to
  the existing base facade;
- allow visible, easily edited driver/project defaults such as
  `core_spacing = 0.3` and template `padding`, while treating them as explicit
  resolved input values that may be overridden for quick tests;
- use clear `ArgumentError`s for unsupported atom workflow inputs where
  practical.

Current driver validation remains origin-centered H. This driver ID does not
approve changing `src/cartesian_base_hamiltonian.jl`; producer-side
one-center atom support is governed separately by `HP-R1-ATOM-*`.

### HP-DRV-ATOM-WIRE-01 — driver atom-to-base-facade wiring

Approved source file:

```text
bin/cartesian_ham_builder.jl
```

Approved behavior:

- the canonical driver may call
  `cartesian_base_hamiltonian(system; basis, hamfile)` for one-center base
  atom construction;
- artifact write/readback uses the existing base facade and existing
  `producer_provenance/` schema;
- no package-internal route-stage helper, terminal basis object, raw-block
  provider, report/status/payload object, or new artifact field is approved.

This base-atom driver wiring ID did not itself approve supplemented atom
Hamiltonians. Supported origin-centered one-center supplemented atom wiring is
now governed by `HP-COMP-SUPPATOM-*`. If a requested atom is outside the
existing base or supported supplemented facade scope, the implementation must
stop at a clear unsupported-input error rather than adding broader atom
construction.

Line budget for `HP-DRV-ATOM-FN-01` plus `HP-DRV-ATOM-WIRE-01`: at most `80`
added `bin` lines, with no committed test, tool, or input-fixture file.

### HP-DRV-ATOM-TEST-01 — base atom driver validation

Approved validation:

- package load;
- origin-centered H base driver artifact write/readback with explicit system
  and one-center basis fields;
- optional ignored negative checks for non-origin atom input, nonneutral
  electron count, mismatched temporary `d` if accepted, or unsupported atom
  input.

No supplemented atom endpoint was approved by this base-atom driver test ID;
supported supplemented atom validation is owned by
`HP-COMP-SUPPATOM-TEST-01`. No translated-atom gate, committed atom fixture,
new committed test file, solver run, artifact-schema validation, or broader
base-atom validation is approved by this ID.

### HP-DRV-ATOM-CLEAN-01 — remove hidden atom `d` driver residue

Approved source file:

```text
bin/cartesian_ham_builder.jl
```

Approved behavior:

- remove the hidden one-center atom basis field `d = vars[:core_spacing]`;
- keep the visible driver atom basis in terms of `ns`, `core_spacing`,
  `radius`, and existing optional public fields only;
- keep public inputs, defaults, overrides, hooks, timing labels, visible stage
  sequence, artifact schema, and driver contract unchanged.

This ID exists because the producer no longer requires public `d` for
one-center atoms. It does not approve source/kernel changes, diagnostics, new
hooks, new timing labels, public input changes, committed tests/fixtures, Cr2
workflow, old `:white_lindsey_low_order` retirement, test/tool route-input
cleanup, or edits outside `bin/cartesian_ham_builder.jl`.

Failure rule: if removing the hidden `d` field requires any visible driver
contract change or producer/source change, make no source commit and report the
blocker.

## Completed Complete-Core-Shell RHF Retirement

Status: completed and closed by `28e9b2c84`. The entries below are historical
deletion authority and no longer authorize source or test work.

This section records the deletion lane in
`complete_core_shell_rhf_retirement.md`. The old complete-core-shell RHF stack
was stale route-era workflow machinery, not a live Cartesian Hamiltonian
producer path.

### HP-RETIRE-CCS-RHF-FN-01 — remove stale RHF payload stack

Status: completed and closed by `28e9b2c84`.

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

## Approved Measurement-Only Authority

These entries authorize ignored measurement/probe work only. They do not
authorize production source edits, committed tests, source files, persistent
objects, metadata/report/status/payload fields, artifacts, public API, or Cr2
workflow support.

### HP-COMP-ANGBOX-AUDIT-01 — angular-balanced shellification geometry audit

Approved scope:

- use ignored `tmp/work` probes only to measure z-axis diatomic
  shellification against the angular-balanced molecular box rule;
- report parent axis physical endpoints and counts, snapped nuclear indices,
  core boxes, molecular inner box, each proposed shared-shell expansion,
  transverse physical scale, low/high longitudinal margins from outer nuclei,
  angular-balance ratios, planned non-boundary and boundary z-extension slab
  stacks, and residual outer mismatch if any;
- classify whether the CR2-style thickness-5 axial slabs are planned angular
  z-extension stacks or unexplained fallback outer mismatch;
- recommend a later source lane only if exact files, functions, forbidden
  surfaces, validation, and failure rules are clear.

This ID does not approve production source edits, shellification repair,
thin-slab lowering changes, driver changes, artifact/schema/reader changes,
route skeleton changes, RG/MWG/IDA/Hamiltonian/raw-block changes, public API,
committed tests/fixtures, Cr2-specific workflow, or Cr2 Hamiltonian runs.

### HP-R3REM-AUDIT-01 — remaining exact-operator allocation audit

Approved scope:

- measure the Cr2 q4 R3/RG exact augmented-operator allocation remaining after
  `954c86cd` and the terminal `G-G` product-workspace optimization;
- separate total wrapper allocation from neutral raw-block construction,
  terminal `G-G` product buffers, unit-nuclear `U_GG` Gaussian-sum
  construction, exact augmented nuclear transforms, route/stage setup, and
  audit/replay overhead;
- use ignored `tmp/work` probes only, with H2/Be2 sanity if needed.

Required outcome:

- classify the dominant remaining allocation bucket;
- recommend a future source lane only if the owner, files, functions,
  forbidden surfaces, validation gates, line budget, deletion/simplification
  expectation, and failure rule are specific enough for a separate docs-only
  amendment.

This ID does not extend `HP-R3GG-FN-01` and does not approve unit-nuclear
`U_GG` Gaussian-sum optimization, route/raw-block setup cleanup, final-basis
`G-G` changes, `G-A`/`A-A` raw-block changes, residual Gaussian algorithm
changes, IDA/MWG changes, Qiu-White semantic changes, parent or terminal-basis
changes, persistent caches/workspaces, artifacts, public API/export, Cr2 facade
support, or Cr2 artifact workflow.

## Rejected Or Deferred

### HP-RES-01 — terminal basis build result — rejected

Do not introduce a persistent terminal-basis result wrapper. The realizer
returns `CartesianTerminalBasisRealization` on success.

### HP-CHANGE-01 — return shell overlap from existing shell plan — rejected/deferred

This can be a helper detail under HP-FN-00, but it is not standalone authority.

### HP-OBJ-03 — generic build-result wrapper — rejected

Do not introduce `CartesianHamiltonianBuildResult`, another payload, or a broad
status wrapper around `CartesianIDAHamiltonian`.

### HP-TEST-01 — new committed terminal smoke — rejected

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
