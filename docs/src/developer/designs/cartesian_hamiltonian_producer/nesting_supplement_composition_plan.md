# Nesting/Supplement Composition

Status: implemented under `HP-ROUTE-RECIPE-FN-01`, `HP-COMP-BASEDIAT-*`,
`HP-COMP-SUPPWL-*`, `HP-COMP-SUPPATOM-*`, `HP-COMP-NS-*`, and
`HP-COMP-WLNS-*`. `HP-ROUTE-RECIPE-TEST-01` is completed evidence with no
continuing permission.

This page is the canonical contract for the supported geometry, nesting, and
supplement composition and for family-selective route recipes. The registry
owns ID lifecycle and source permission. Internal route/stage inventory and
carrier semantics are owned by
[Route/stage metadata](route_stage_metadata_contract.md).
Direct-core parity and one-center physical extent are owned separately by
[Public ns direct-core side parity](public_ns_core_side_parity.md) and
[R1 one-center base atoms](r1_one_center_base_atoms.md).

## Route Recipe

`cartesian_recipe(route_inputs)` selects one real route family:

- `:pqs_source_box` builds the `source_box` subrecipe and does not require
  inactive `white_lindsey_*` input fields;
- `:white_lindsey_low_order` builds the `white_lindsey` subrecipe and does not
  require inactive PQS source-box fields.

The inactive output subrecipe is `nothing`. Already-precomposed recipes with
both stable field names remain accepted only through the existing compatibility
path. This family selection removes inactive vocabulary; it does not merge the
PQS and White-Lindsey algorithms or make either family removable.

## Composition Contract

The producer composes three choices:

```text
geometry:   origin-centered atom | homonuclear z-axis diatomic
nesting:    :pqs | :wl
supplement: off | on
```

All eight bounded cells use the shared producer:

| Geometry | Supplement | `:pqs` | `:wl` |
| --- | --- | --- | --- |
| atom | off | implemented | implemented |
| atom | on | shared residual-GTO/MWG augmentation | shared residual-GTO/MWG augmentation |
| z-axis diatomic | off | implemented | native WL terminal realization |
| z-axis diatomic | on | shared residual-GTO/MWG augmentation | shared residual-GTO/MWG augmentation |

The WL diatomic path remains subject to its compact retained-basis and boundary
count contracts. Composition does not imply equal PQS/WL dimensions or general
molecular support.

## Public Size Rules

Public `ns` is the requested cube/source/nesting size. Route-local `q` is
derived only after nesting selection:

```text
nesting = :pqs  -> q = ns
nesting = :wl   -> q = ns - 2
```

Legacy public `q`, where still accepted, is normalization compatibility only.
If both values are supplied, they must agree with the selected nesting. New
inputs and provenance use `ns`.

White-Lindsey z-axis diatomics reject normalized `ns < 4` before route
construction because the complete-shell inner box requires the stricter floor;
the generic WL floor remains `ns >= 3`. WL retained support may saturate across
working `ns` values when physical parent extent dominates. Changing `ns` can
still change route-local `q`, decomposition, and row order without changing
final support. Equal `ns` does not promise equal PQS/WL dimensions.

Direct nucleus-centered core side is independent of route-local `q`:

```text
direct_core_side = isodd(ns) ? ns : ns + 1
```

This oddization does not apply to boundary retained counts.

## Geometry And Supplement

Base diatomics require:

- two equal atom symbols used as provenance labels;
- two equal, finite, positive, integer-valued nuclear charges;
- two finite, distinct centers on the Cartesian z axis;
- explicit nonnegative `nup` and `ndn` satisfying
  `nup + ndn == sum(nuclear_charges)`.

No element table supplies charge, electron count, spin, basis, or geometry.

After either nesting family produces a `CartesianTerminalBasisRealization`,
supplemented atoms and diatomics share this downstream path:

```text
terminal basis
-> supplement loading
-> owner-local residual Gaussian selection
-> exact augmented one-body operators
-> residual MWG/IDA interaction
-> Hamiltonian assembly and existing artifact
```

One-center inputs use `legacy_atomic_gaussian_supplement(...)`; two-center
inputs use `legacy_bond_aligned_diatomic_gaussian_supplement(...)`. The shared
residual-GTO/MWG producer does not branch into a second WL implementation.

For one-center atoms, `basis.radius` is physical parent-extent authority.
Parent counts follow the atomic mapping and spacing policy and are bounded
below by the direct-core side. `ns` is resolution/nesting input, not box extent.

## Ownership And Failure

`src/cartesian_base_hamiltonian.jl` owns input normalization, family-selective
base route inputs, public `ns`/derived `q`, WL rejection, one-center sizing, and
shared supplemented composition. `src/pqs_source_box_route_driver_helpers.jl`
owns route recipe selection. `bin/cartesian_ham_builder.jl` owns visible driver
input and dispatch, not a second composition algorithm. Shared residual
augmentation remains in
`src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl`.

Reject heteronuclear, translated, non-z-axis, or general molecular geometry;
nonneutral or noninteger-charge all-electron systems; inconsistent `ns`/`q`;
WL diatomic `ns < 4`; and supplement basis-count or owner mismatch before
expensive construction where practical. Do not add geometry-specific
Hamiltonian builders when the shared terminal and residual-GTO/MWG boundaries
express the case.

## Non-Goals

This contract does not authorize translated atoms, heteronuclear or generally
oriented molecules, ECP, EGOI, solver/RHF, exchange, Cr2-specific workflow,
route or shellification changes, terminal-lowering changes, numerical policy
changes, raw-block or RG/MWG/IDA interaction-policy changes, new driver
controls, route diagnostics, status/report payloads, public exports, artifact
schemas, readers, fixtures, compatibility layers, WL materialization deletion,
or retirement of route/stage machinery. WL compact retained-basis,
shellification, and thin-slab behavior remain separate subsystem contracts.
