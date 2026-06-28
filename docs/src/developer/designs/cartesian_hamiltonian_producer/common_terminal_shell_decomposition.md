# Common Terminal Shell Decomposition

Status: approved narrow audit/source authority under
`HP-COMP-SHELLGEOM-FN-01` and `HP-COMP-SHELLGEOM-TEST-01`, with the
z-axis diatomic same-function/same-argument correction approved under
`HP-COMP-SHELLGEOM-DIAT-FN-01` and `HP-COMP-SHELLGEOM-DIAT-TEST-01`, and
narrow common thin-slab stack lowering repair approved under
`HP-COMP-THINSLAB-FN-01` and `HP-COMP-THINSLAB-TEST-01`.

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
direct core and not a real shell. `:direct_midpoint_slab` and
`:outer_mismatch_slab` must be lowered through the same compact thin-slab
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

An outer-mismatch region of thickness `t <= ns` should be decomposed or
realized as an oriented stack of compact one-slice slab functions, with scale
about `t * ns * ns`. If slab thickness exceeds `ns`, implementation must stop
and report the condition. A future policy could approve a whole-block
`ns x ns x ns` compression, or treat that case as a setup error, but this lane
does not choose silently between those options.

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
