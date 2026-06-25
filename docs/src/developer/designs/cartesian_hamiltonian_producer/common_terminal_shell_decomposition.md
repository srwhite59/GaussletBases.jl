# Common Terminal Shell Decomposition

Status: approved narrow audit/source authority under
`HP-COMP-SHELLGEOM-FN-01` and `HP-COMP-SHELLGEOM-TEST-01`.

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

Public `ns` is the common user-facing size. Direct nucleus-centered core side
comes from public `ns` under `HP-COMP-NSCORE-*`. PQS may derive `q = ns` for
PQS retained/source-mode policy. White-Lindsey may derive its inner side
`ns - 2` for WL boundary contraction policy. Neither PQS `q` nor WL inner side
is the common shell/core ownership authority.

## Approved IDs

- `HP-COMP-SHELLGEOM-FN-01` - common terminal shell decomposition audit and
  narrow cleanup.
- `HP-COMP-SHELLGEOM-TEST-01` - validation gates.

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

## Forbidden

This amendment does not approve:

- driver changes;
- public input changes;
- route skeleton redesign;
- terminal lowering redesign;
- retained-unit record changes;
- retained-unit transform changes;
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

No committed test file, committed fixture, driver contract test,
solver/RHF/ECP/EGOI validation, route-diagnostic validation, artifact schema
validation, or Cr2 fixture is approved.

## Failure Rule

If common shell decomposition cannot be made route-family-free without
changing terminal lowering, retained-unit records, PQS retained-mode
realization, WL boundary coefficient construction, route skeleton semantics,
artifact schema, or driver inputs, make no source commit and report the exact
blocker.

Separate follow-up: if the z-axis diatomic central-gap/contact policy turns
out to mix common shell decomposition with family-specific retained geometry,
that needs a later docs-only amendment. Do not hide that inside the atom shell
cleanup.
