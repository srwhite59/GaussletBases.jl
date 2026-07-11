# Common Terminal Shell Decomposition

Status: implemented common shellification, angular-z-extension, neutral
face-product, and compact thin-slab subsystem. The registry records the
`HP-COMP-SHELLGEOM-*`, `HP-COMP-ANGBOX-*`, `HP-COMP-FACEPROD-*`,
`HP-COMP-THINSLAB-*`, and `HP-COMP-THINSLAB-META-*` source/test surfaces.

## Problem

PQS and White-Lindsey now share the public composition model:

```text
geometry -> nesting -> optional supplement
```

They should also share the first geometric operation. For one-center atoms,
that operation is:

```text
parent lattice + nuclear center + direct core side
-> direct nucleus-centered core
-> shell 1, shell 2, ...
-> owned shell support rows
```

That operation is not PQS-specific or White-Lindsey-specific. If each
construction family computes its own shells, or if route-family variables
change shell/core ownership before lowering, the two families will diverge for
accidental reasons.

## Two Geometry Levels

There are two different geometry levels. They must not be flattened into one
route-stage idea.

### Common terminal shell decomposition

This is shared by PQS and White-Lindsey:

```text
parent lattice
+ nuclear centers
+ public/common direct core side
-> direct core regions
-> terminal shell regions
-> owned support rows
```

It owns coverage, ordering, direct-core centering, and shell-owned supports.
It must be route-family-free.

### Family-specific retained-construction geometry

This starts only after common shell regions exist.

PQS consumes the common shell and builds retained functions using a full local
source box:

```text
common shell support
+ full source CPB
-> boundary COMX/product-mode selection
-> restrict rows to shell-owned support
-> shell-local Lowdin
```

White-Lindsey consumes the common shell and decomposes its boundary into local
units:

```text
common shell support
-> faces / edges / corners / strata
-> product-of-1D contractions
-> compact WL retained functions
```

The full source CPB in PQS and the face/edge/corner units in White-Lindsey are
second-level retained-construction geometry. They are not separate first-step
shellification algorithms.

## Decision

Common terminal shell decomposition must be a shared route-family-free
operation. Route-family code may consume the common shell records, but it must
not recompute shell ownership or reinterpret owned support rows.

For both one-center atoms and z-axis diatomics, PQS and White-Lindsey must call
the same common shell decomposition function with the same first-step
arguments when the public system, parent axes, nuclear centers, direct core
side, and public `ns` match. The construction family must not alter those
arguments before common shell decomposition.

Public `ns` is the common user-facing size. Direct nucleus-centered core side
comes from public `ns` under `HP-COMP-NSCORE-*`. PQS may derive `q = ns` for
PQS retained/source-mode policy. White-Lindsey may derive its inner side
`ns - 2` for WL boundary contraction policy. Neither PQS `q` nor WL inner side
is the common shell/core ownership authority.

## Approved IDs

- `HP-COMP-SHELLGEOM-FN-01` - common terminal shell decomposition audit and
  narrow cleanup.
- `HP-COMP-SHELLGEOM-TEST-01` - validation gates.
- `HP-COMP-SHELLGEOM-DIAT-FN-01` - z-axis diatomic same-function/same-argument
  common shell entry cleanup.
- `HP-COMP-SHELLGEOM-DIAT-TEST-01` - diatomic parity validation gates.
- `HP-COMP-OUTERMM-FN-01` / `HP-COMP-OUTERMM-TEST-01` - superseded
  outer-mismatch-only subset; do not implement as a separate source lane.
- `HP-COMP-FACEPROD-FN-01` - neutral compact face-product terminal helper.
- `HP-COMP-FACEPROD-TEST-01` - face-product helper reuse/parity validation.
- `HP-COMP-THINSLAB-FN-01` - z-axis diatomic thin-slab stack compact lowering
  repair for both PQS and White-Lindsey.
- `HP-COMP-THINSLAB-TEST-01` - thin-slab stack compact lowering validation
  gates.

## Approved Source Surface

Approved later source files:

```text
src/cartesian_shellification/terminal_geometry.jl
src/pqs_source_box_route_driver_helpers.jl
```

`src/cartesian_shellification/terminal_geometry.jl` owns the common
route-family-free shell/core region decomposition.

`src/pqs_source_box_route_driver_helpers.jl` is approved only for narrow
caller plumbing and summary/provenance wording needed to pass common shell
inputs before selecting PQS or White-Lindsey lowering.

## Approved Behavior

The later source pass may:

- audit whether the first-step shell/core region construction is already
  identical for `nesting = :pqs` and `nesting = :wl`;
- remove route-family branching from common terminal shell decomposition if it
  exists;
- rename or locally clarify shellifier parameters and summary labels so common
  shell geometry does not appear to be governed by PQS `q` or WL inner side;
- keep direct-core side tied to public `ns` through `HP-COMP-NSCORE-*`;
- keep deterministic shell region order, owned support rows, and coverage
  checks unchanged;
- leave PQS and White-Lindsey lowering/realization families separate after
  common shell records are produced.

For one-center atoms, same public system, parent extent, and `ns` must produce
the same direct core and shell-owned support regions before family-specific
lowering.

For z-axis diatomics, this lane may audit whether central-gap and shared-shell
planning use the same common shell decomposition. It must not change central
gap/contact policy unless the fix is the same route-family-free shell input
cleanup and does not touch lowering, retained units, or WL/PQS realization.

## Diatomic Same-Function/Same-Argument Requirement

`HP-COMP-SHELLGEOM-DIAT-FN-01` promotes the z-axis diatomic part from audit to
narrow cleanup authority.

For a fixed public z-axis diatomic system, parent axes, public `ns`, direct
core side, nuclear centers, and bond axis, `nesting = :pqs` and `nesting = :wl`
must enter `raw_terminal_geometry(...)` or its common replacement with the same
first-step arguments. In particular:

- direct core side is the public-`ns` value from `HP-COMP-NSCORE-*`;
- public `ns` is the common shell-size input;
- PQS `q` and the WL inner side are family-specific retained-construction
  inputs after common shell records exist;
- central-gap, contact-core, shared-shell, and outer-mismatch region ownership
  are common shell decomposition facts, not retained-construction facts.

This authority allows only the caller plumbing and parameter naming needed to
make the diatomic common shellifier entry route-family-free. If
`raw_terminal_geometry(...)` currently uses a parameter named `q` for common
central-gap or shared-shell decisions, the source pass should rename or
reinterpret that parameter as common `ns` at the shellifier boundary rather
than passing PQS `q` or WL inner side.

This does not approve changing the central-gap/contact algorithm itself. It
approves only making both construction families use the same algorithm with
the same common inputs before lowering.

## Angular-Balanced Molecular Boxes

The shared z-axis diatomic shellifier should not create large axial leftovers
merely because rectangular index-layer growth keeps the bond-axis margin too
short relative to the transverse scale, and then stops when a transverse axis
reaches the parent boundary. That is a shellification geometry problem, not a
lowering problem.

For each shared molecular shellification step around a z-axis diatomic,
shellification should compute an angular-balanced target box from the outer
nuclei in physical parent-axis coordinates. In each bond-axis/transverse
plane, compare:

```text
longitudinal margin = physical distance from the outer nucleus to the box end
transverse scale    = physical distance from the bond axis to the box side
```

The target box should keep the longitudinal margin comparable to the selected
transverse scale. If the `x` and `y` transverse scales differ, use the smaller
scale as the conservative angular-resolution guard unless a later amendment
approves a different convention. This is the operational outer-nucleus
`45` degree rule. It must be evaluated in physical coordinates, not raw index
counts.

The current audit shows the ordinary index-layer shared-shell bodies are
underextended in `z` relative to the angular-balance rule. When the
angular-balanced target for a shared-shell step requires bond-axis-only
extension beyond the ordinary body, shellification should emit that difference
as native axial thin-slab stacks. The ordinary body plus planned z-extension
slabs, not the ordinary body alone, realizes the angular-balanced target
coverage. These pieces must not be left as mysterious outer-mismatch identity
regions at the end. The thin-slab concept is common and applies to:

- central midpoint slabs between atom-local regions;
- planned non-boundary angular z-extension slabs produced inside a larger
  parent growth step;
- planned boundary angular z-extension slabs at the low/high parent boundary;
- unexpected outer-mismatch fallback slabs, which should become rare and
  diagnostic after angular-balanced shellification.

Planned angular z-extension slabs should carry native metadata instead of
requiring later code to parse labels such as `z_low_outer_mismatch_slab`:

```text
slab_kind = :angular_z_extension_slab
slab_normal_axis
slab_side
slab_thickness
slab_stack_index
slab_stack_count
bond_axis
reference_nucleus_index
angular_balance_rule = :outer_nucleus_45_degree
longitudinal_margin_physical
transverse_scale_physical
angular_extension_physical
```

For a planned angular z-extension, total extension thickness larger than `ns`
should be split into multiple ordered compact slab units, each with thickness
`<= ns`, unless a later whole-block compression policy is approved. An
unplanned fallback slab with thickness `> ns` remains a setup/shellification
failure under the existing thin-slab guardrail.

`HP-COMP-ANGBOX-FN-01` approves this shellification repair in
`src/cartesian_shellification/terminal_geometry.jl`, with only narrow summary
or caller plumbing if directly required. `HP-COMP-ANGBOX-TEST-01` approves
ignored geometry probes for H2/Be2/Cr2-style z-axis diatomics. Lowering planned
z-extension slabs remains deferred to `HP-COMP-THINSLAB-*`.

## Neutral Face-Product Helper

Compact thin slabs and White-Lindsey facets share a numerical shape: two
active axes use retained one-dimensional contractions, while the normal axis is
fixed to one or more parent indices. A thickness-1 slab is one face-like
product block. A thickness-`t` slab is an ordered stack of face-like blocks.

The reusable coefficient assembly is not PQS-owned and not White-Lindsey-owned.
`HP-COMP-FACEPROD-FN-01` approves a private/module-internal neutral helper in:

```text
src/cartesian_final_basis_realization/terminal_face_product_blocks.jl
```

with the include in:

```text
src/cartesian_final_basis_realization/CartesianFinalBasisRealization.jl
```

and narrow consumers in:

```text
src/cartesian_final_basis_realization/white_lindsey_terminal_basis_realization.jl
src/cartesian_final_basis_realization/pqs_terminal_basis_realization.jl
```

The helper should reuse the existing numerical primitives:

```text
_nested_doside_1d(...)
_nested_face_product(...)
```

It should support normal axes `:x`, `:y`, and `:z`; one fixed normal-axis
index for a face-like block; an ordered stack of fixed normal-axis indices for
a thickness-`t` slab; and a caller-supplied retained count, normally public
`ns`.

White-Lindsey facet terminal realization should be refactored to use this
neutral helper. That is the reuse proof: if the helper cannot serve current WL
facets and future thin slabs without changing numerical semantics, the source
patch should stop and report whether the blocker is the helper signature,
support-record shape, retained-unit metadata, or terminal-realization
ownership.

Do not put the shared helper in
`white_lindsey_terminal_basis_realization.jl`; do not create a PQS-specific
thin-slab projection path; do not pretend thin slabs are WL boundary strata
for naming convenience; and do not duplicate `_nested_face_product(...)`
assembly.

## Thin-Slab Lowering

Common shell decomposition owns thin-slab support regions, but it does not
make those rows retained basis functions.

The audited bad paths are:

```text
:direct_midpoint_slab
-> :direct_slab_identity_cpb
-> direct retained unit
-> full identity terminal block

:outer_mismatch_slab
-> :direct_boundary_slab_identity_cpb
-> direct retained unit
-> full identity terminal block
```

For z-axis diatomics this is not approved under either `PQSLowering` or
`WhiteLindseyLowering`. A thin slab is a boundary-face-like object, not a
direct core and not a real shell. Midpoint slabs, planned non-boundary
z-extension slabs, planned boundary z-extension slabs, and fallback
outer-mismatch slabs must be lowered through the same compact thin-slab
function for both construction families, with the same terminal region, public
`ns`, slab normal axis, slab thickness, and native source/support facts, when
those facts are sufficient.

This sameness rule is deliberately limited to thin slab stacks. Real shell
regions still diverge after common shellification: PQS uses full source-box
shell projection, while White-Lindsey uses face/edge/corner product-of-1D
contractions. Direct nucleus-centered core and atom-contact core sectors remain
identity sectors.

The compact unit slice should retain the oriented scale:

```text
ns x ns x 1
```

after standard one-dimensional COMX/product compression. They must not become
full identity slabs. The `1` is along the slab normal. The support rows remain
owned and disjoint, but the retained functions are compact slab functions, not
the support rows themselves.

An unplanned outer-mismatch fallback region of thickness `t <= ns` should be
decomposed or realized as an oriented stack of compact one-slice slab
functions, with scale about `t * ns * ns`. Planned angular z-extensions with
total thickness greater than `ns` should be chunked into ordered slab units of
thickness `<= ns`; an unplanned fallback slab with thickness greater than `ns`
must still stop and report the condition. A future policy could approve a
whole-block `ns x ns x ns` compression, or treat that fallback case as a setup
error, but this lane does not choose silently between those options.

Approved source files for this repair are:

```text
src/cartesian_terminal_lowering/selection.jl
src/cartesian_terminal_lowering/region_contracts.jl
src/cartesian_retained_units/lower_contract_units.jl
src/cartesian_retained_unit_transform_contracts/unit_contracts.jl
src/cartesian_final_basis_realization/pqs_terminal_basis_realization.jl
src/cartesian_final_basis_realization/white_lindsey_terminal_basis_realization.jl
```

Optional files are limited to native slab-axis/thickness metadata or existing
summary/record plumbing if directly required:

```text
src/cartesian_shellification/terminal_geometry.jl
src/pqs_source_box_route_driver_helpers.jl
src/pqs_source_box_diatomic_complete_core_shell.jl
```

This repair deliberately treats WL the same as PQS for thin slabs. It approves
only the narrow retained-unit, transform-contract, and terminal-realization
work needed to consume the shared thin-slab retained object. If existing
compact slab machinery cannot compact the slab from the available facts,
source work must stop and report the missing native fact. Do not infer slab
normal or thickness from role strings such as `z_low_outer_mismatch_slab`; add
native metadata under the optional shellification surface if needed.

## Forbidden

This amendment does not approve:

- driver changes;
- public input changes;
- route skeleton redesign;
- terminal lowering redesign beyond the narrow common thin-slab repair in
  `HP-COMP-THINSLAB-FN-01`;
- angular-balanced shellification source changes outside
  `src/cartesian_shellification/terminal_geometry.jl` except narrow approved
  summary/caller plumbing;
- retained-unit record changes beyond the shared thin-slab retained object;
- retained-unit transform changes beyond the shared thin-slab transform
  contract;
- PQS source-box retained-mode realization changes;
- WL face/edge/corner coefficient or retained-basis changes;
- direct-core parity changes beyond `HP-COMP-NSCORE-*`;
- central-gap/contact policy redesign;
- artifact schema, manifest, reader, or provenance expansion;
- Hamiltonian, one-body, IDA, MWG, Residual Gaussian, raw-block, solver, ECP,
  or Cr2 workflow changes;
- old WL materialization revival;
- committed tests or fixtures.

## Validation

`HP-COMP-SHELLGEOM-TEST-01` approves only:

- `git diff --check`;
- package load;
- focused audit showing the one-center atom common shell decomposition is
  route-family-free for the same public `ns`, parent extent, and center;
- same-`ns` PQS/WL one-center atom direct core and shell-owned support counts
  match before family-specific lowering;
- same-`ns` PQS/WL one-center atom base artifact/readback still works for a
  bounded fixture;
- H2 or Be2 smoke to confirm the diatomic path still constructs;
- existing H2 Residual Gaussian endpoint smoke only if touched code crosses
  supplemented path;
- no Cr2 run.

`HP-COMP-SHELLGEOM-DIAT-TEST-01` additionally approves:

- focused audit showing z-axis diatomic PQS/WL calls enter the common
  shellifier with the same parent axes, nuclear centers, direct core side,
  public `ns`, and bond axis;
- same-`ns` PQS/WL z-axis diatomic direct core, central-gap/contact,
  shared-shell, and outer-mismatch region counts match before
  family-specific lowering for a bounded H2 or Be2 fixture;
- base artifact/readback smoke for the same bounded fixture under
  `nesting = :pqs` and `nesting = :wl` if the WL retained-basis path remains
  available.

`HP-COMP-THINSLAB-TEST-01` additionally approves:

- focused audit showing `:direct_midpoint_slab` and `:outer_mismatch_slab`
  do not lower to `:direct_slab_identity_cpb` or
  `:direct_boundary_slab_identity_cpb` under either lowering family;
- focused audit showing PQS and White-Lindsey call the same compact thin-slab
  lowering function with the same region/public-`ns` inputs for matched
  slab regions or stack slices;
- bounded H2 or Be2 artifact/readback under `nesting = :pqs` and
  `nesting = :wl` where midpoint slabs and outer-mismatch slabs are present or
  explicitly probed;
- retained thin-slab count is compact relative to support count, with normal
  oriented unit-slice target `ns x ns x 1` and thickness-`t <= ns`
  outer-mismatch scale about `t * ns * ns`;
- existing H2 Residual Gaussian endpoint smoke if touched code crosses
  supplemented construction;
- optional CR2 user-run inventory only, not a committed gate.

No committed test file, committed fixture, driver contract test,
solver/RHF/ECP/EGOI validation, route-diagnostic validation, artifact schema
validation, or Cr2 fixture is approved.

`HP-COMP-ANGBOX-TEST-01` approves only ignored geometry probes that report:

- parent axis physical endpoints and counts;
- snapped nuclear indices and direct/core boxes;
- molecular inner box;
- each proposed shared-shell expansion;
- transverse physical scale;
- low/high longitudinal margins from the outer nuclei;
- angular-balance ratios;
- planned non-boundary and boundary z-extension slab stacks;
- residual outer mismatch, if any.

The audit must classify whether CR2-style axial slabs are planned angular
z-extension stacks or unexplained fallback outer mismatch. It should show
planned z-extension support, zero residual z mismatch after classification,
and PQS/WL geometry parity. No artifact/readback is required while lowering is
intentionally deferred. This validation does not approve committed tests,
Cr2 Hamiltonian runs, or driver changes.

## Failure Rule

If common shell decomposition cannot be made route-family-free without
changing terminal lowering, retained-unit records, PQS retained-mode
realization, WL boundary coefficient construction, route skeleton semantics,
artifact schema, or driver inputs, make no source commit and report the exact
blocker.

If making z-axis diatomic PQS/WL use the same shellifier with the same
first-step arguments requires changing the central-gap/contact algorithm rather
than just its route-family-independent inputs, stop and request a separate
docs-only amendment.

If the shared compact thin-slab lowering cannot be built from existing native
facts without broad route skeleton changes, artifact/schema changes, driver
inputs, or general real-shell policy changes, make no source commit and report
the missing fact. If slab thickness exceeds `ns`, stop and request a separate
whole-block compression or setup-error policy.
