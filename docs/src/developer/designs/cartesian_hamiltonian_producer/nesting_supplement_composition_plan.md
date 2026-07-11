# Nesting/Supplement Composition

Status: implemented for the bounded composition and input lanes owned here:
`HP-COMP-BASEDIAT-*`, `HP-COMP-SUPPWL-*`,
`HP-COMP-SUPPATOM-*`, `HP-COMP-NS-*`, and
`HP-COMP-WLNS-*`.

This page is the canonical contract for how supported geometry, nesting, and
supplement choices compose. Direct-core side parity is owned separately by
[Public ns direct-core side parity](public_ns_core_side_parity.md), and
one-center physical extent is owned by
[R1 one-center base atoms](r1_one_center_base_atoms.md).

## Composition Contract

The producer accepts three successive choices:

```text
geometry:   origin-centered atom | homonuclear z-axis diatomic
nesting:    :pqs | :wl
supplement: off | on
```

All eight bounded cells are implemented through the shared producer:

| Geometry | Supplement | `:pqs` | `:wl` |
| --- | --- | --- | --- |
| atom | off | implemented | implemented |
| atom | on | implemented through shared residual-GTO/MWG augmentation | implemented through shared residual-GTO/MWG augmentation |
| z-axis diatomic | off | implemented | implemented through native WL terminal records |
| z-axis diatomic | on | implemented through shared residual-GTO/MWG augmentation | implemented through shared residual-GTO/MWG augmentation |

The WL diatomic path remains subject to the separately owned compact retained
basis and boundary-count contracts. Implemented composition does not imply
identical PQS/WL dimensions or general molecular support.

## Public Size Rules

Public `ns` is the requested cube/source/nesting size. Route-local `q` is
derived only after selecting the nesting family:

```text
nesting = :pqs  -> q = ns
nesting = :wl   -> q = ns - 2
```

Legacy public `q`, where still accepted, is normalization compatibility only.
When both `ns` and `q` are supplied, they must agree with the selected
nesting. New inputs, examples, and provenance use `ns`.

For WL z-axis diatomics, normalized `ns < 4` is rejected before route
construction. The generic WL normalization floor remains `ns >= 3`; the
stricter diatomic floor is required by the complete-shell inner box.

WL retained support may saturate across working `ns` values when physical
parent extent dominates. A changed `ns` can still change route-local `q`,
block decomposition, and row ordering without changing the final retained
support. This is not an ignored-input bug, and equal `ns` does not promise
equal PQS/WL dimensions.

Direct nucleus-centered core side is derived separately from public `ns`:

```text
direct_core_side = isodd(ns) ? ns : ns + 1
```

That oddization does not apply to boundary retained counts. The full rule is in
[Public ns direct-core side parity](public_ns_core_side_parity.md).

## Geometry And Supplement Rules

Base diatomics require explicit:

- two equal atom symbols used as provenance labels;
- two equal, finite, positive, integer-valued nuclear charges;
- two finite, distinct centers on the Cartesian z axis;
- explicit nonnegative `nup` and `ndn` with
  `nup + ndn == sum(nuclear_charges)`.

No element table supplies charge, electron count, spin, basis, or geometry.

After geometry and nesting produce a
`CartesianTerminalBasisRealization`, supplemented atoms and diatomics use the
same downstream path:

```text
terminal basis
-> supplement loading
-> owner-local residual Gaussian selection
-> exact augmented one-body operators
-> residual MWG/IDA interaction
-> Hamiltonian assembly and existing artifact
```

One-center inputs use `legacy_atomic_gaussian_supplement(...)`. Two-center
inputs use `legacy_bond_aligned_diatomic_gaussian_supplement(...)`.
White-Lindsey and PQS differ upstream of the terminal-basis boundary; the
residual-GTO/MWG producer does not branch into a second WL implementation.

For one-center atoms, `basis.radius` is physical parent-extent authority.
Parent counts are derived from the atomic mapping and spacing policy, then
bounded below by the direct-core side. `ns` is resolution/nesting input, not
a direct box-size substitute. See
[R1 one-center base atoms](r1_one_center_base_atoms.md).

## Source Ownership

Implemented source owners are:

- `src/cartesian_base_hamiltonian.jl` for system/basis normalization, base
  diatomic validation, public `ns` and derived `q`, WL diatomic rejection,
  one-center parent sizing, and shared supplemented composition;
- `bin/cartesian_ham_builder.jl` for the canonical visible `ns` input and
  atom/diatomic supplemented dispatch;
- `src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl` for
  the shared residual-GTO augmentation consumed by both nesting families.

The driver does not own a separate composition algorithm.

## Validation And Evidence

No target lane added a dedicated committed fixture. Accepted bounded
construction/readback and rejection evidence is recorded in
`docs/src/developer/pqs_manager_running_log.md`:

| Contract | Implementation | Accepted evidence |
| --- | --- | --- |
| base homonuclear diatomic | `095a89d41` | Pass 139 |
| supplemented WL diatomic | `4cfb47ace` | Pass 141 |
| supplemented atom | `2e9818c90` | Pass 142 |
| public `ns` normalization | `5f0185a16` | Pass 145 |
| atom physical parent sizing | `18d683575` | Pass 146 |
| WL diatomic `ns` guard | `50327dc1e` | Pass 148 |

Current downstream regression owners include
`test/driver_public/cartesian_base_hamiltonian_runtests.jl` for the public
base facade and
`test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl` for the
supplemented path. They do not replace the historical per-lane evidence above.

## Failure Behavior

Reject unsupported or inconsistent input before expensive construction where
practical. In particular, reject:

- heteronuclear, translated, non-z-axis, or general molecular geometry;
- nonneutral or noninteger-charge all-electron systems;
- inconsistent `ns`/`q`;
- WL diatomic `ns < 4`;
- supplement basis count or owner mismatch.

Do not add a geometry-specific Hamiltonian builder when the shared terminal
basis and residual-GTO/MWG boundaries can express the case.

## Non-Goals

These contracts do not authorize:

- translated atoms, heteronuclear or generally oriented molecules;
- ECP, EGOI, solver/RHF, exchange, or Cr2-specific workflow;
- route, shellification, terminal-lowering, raw-block, RG/MWG/IDA, or
  interaction-policy changes;
- new driver controls, route diagnostics, status/report payloads, public
  exports, artifact schemas, readers, fixtures, or compatibility layers.

The WL compact retained-basis, shellification, and thin-slab families remain
separate subsystem contracts.
