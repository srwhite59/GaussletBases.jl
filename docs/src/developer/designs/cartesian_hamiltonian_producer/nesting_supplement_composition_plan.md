# Nesting/Supplement Composition Plan

This note freezes the target shape for the canonical Cartesian Hamiltonian
producer after `nesting = :wl` became a real construction-family input.

It is a planning and authority-boundary amendment. The initial explicit
composition cells are now promoted under the IDs listed below; deferred
geometry, solver, ECP, export, and Cr2-specific work remain candidate-only.

## Target Contract

The long-term supported producer shape is three successive user choices that
compose:

```text
geometry:   atom | z-axis diatomic
nesting:    :pqs | :wl
supplement: off | on
```

The driver may expose these as visible public construction choices, but it must
not encode the current partial implementation as the permanent contract. The
implementation target is a common staged construction:

```text
public system / basis / optional supplement contract
-> geometry normalization
-> nesting-specific parent and terminal-basis realization
-> common base product, unit-nuclear, IDA, and Hamiltonian assembly
-> optional residual-GTO/MWG augmentation
-> existing Hamiltonian artifact
```

`nesting` may affect parent construction, shellification, retained rules, and
terminal-basis realization. After it produces a
`CartesianTerminalBasisRealization`, downstream base operators, residual
Gaussian augmentation, IDA/MWG interaction, and artifact writing should not
branch on whether the terminal basis came from PQS or White-Lindsey.

`supplement` may affect supplement loading, residual Gaussian basis selection,
exact augmented operator transforms, MWG descriptors, residual-containing IDA
blocks, and supplemented artifact provenance. It must not resurrect route
reports, route-stage controls, pair/assembly public stages, or driver-specific
helper choreography.

## Current Matrix

Current status is intentionally explicit. Unsupported cells must fail clearly;
they must not be hidden by driver defaults or mislabeled artifact provenance.

| Geometry | Supplement | `nesting = :pqs` | `nesting = :wl` |
| --- | --- | --- | --- |
| atom | off | implemented for explicit origin-centered all-electron base atoms, with H as the committed endpoint | implemented for one-center base atoms through the WL terminal-basis seam |
| atom | on | approved implementation lane through the common Residual Gaussian path | approved implementation lane through the common Residual Gaussian path |
| z-axis diatomic | off | implemented for explicit homonuclear z-axis all-electron inputs | mechanically implemented through native WL terminal records; compact retained-basis correction approved under `HP-WLDIAT-COMPACT-*` |
| z-axis diatomic | on | supported for explicit homonuclear z-axis diatomics through the residual-GTO/MWG path | approved implementation lane through the same RG/MWG boundary after WL base terminal realization |

## Common Boundary Rules

- Atoms and diatomics must share the same producer workflow after
  geometry/shellification normalization. Atom-only Hamiltonian builders,
  atom-only materialization paths, and atom-specific artifact shapes are
  forbidden.
- PQS and White-Lindsey must converge to the same terminal-basis boundary:
  a `CartesianTerminalBasisRealization` with disjoint owned terminal supports.
- Residual-GTO/MWG supplementation consumes the terminal basis, public system
  facts, parent bundles, supplement specification, and same-construction base
  operators. It must not consume route skeletons, route reports, or old WL
  H1/H1+J materialization objects.
- Artifact provenance records the public choices (`geometry`, `nesting`, and
  supplement state) truthfully. It must not infer route identity from helper
  file names or preserve PQS labels for WL artifacts.
- Public `ns` is the shared source/cube/nesting size input across composition
  cells. Route-local `q` is derived from `ns` and `nesting`, so the driver does
  not expose different cube-size meanings for PQS and White-Lindsey.
- White-Lindsey z-axis diatomics require enough public `ns` to form the
  complete-shell inner box used by the WL terminal shellification path.
  `HP-COMP-WLNS-*` approves early rejection of normalized `nesting = :wl`,
  `Natom = 2`, `ns < 4` in the base producer.
- WL diatomic retained support may saturate across `ns` ranges when the
  physical parent extent dominates. A changed public `ns` may still change
  route-local `q`, block decomposition, and row order without changing the
  final retained support set, so equal public `ns` is not by itself a fair
  PQS/WL retained-basis comparison.
- The current WL diatomic boundary-stratum identity realization is a
  mechanical endpoint, not the compact WL retained-basis contract. The
  `HP-WLDIAT-COMPACT-*` lane makes public `ns` the fair starting input for
  PQS/WL construction comparison by requiring compact retained columns from
  one-dimensional contractions, while still allowing legitimate WL-specific
  geometry/contact differences.
- The canonical driver remains compact and copyable. It can print public
  contracts and coarse physics-stage timings, but it must not grow stop-after
  controls, raw-provider switches, route diagnostics, allocation probes, or
  status/report payloads.

## Dependency Order

### 1. White-Lindsey Diatomic Base

Status: approved for implementation under `HP-COMP-WLDIAT-FN-01` and
`HP-COMP-WLDIAT-TEST-01`.

Goal: make `nesting = :wl`, `supplement = off` work for z-axis diatomic base
artifacts by producing native WL terminal records and the common
`CartesianTerminalBasisRealization`.

This should extend the WL terminal-basis seam from one-center atoms to
two-center z-axis diatomics. It should not adapt the old WL H1/H1+J
materialization path, change the canonical driver contract, or add route
diagnostics.

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

`src/cartesian_base_hamiltonian.jl` is approved only for narrow staged/facade
wiring needed by the WL z-axis diatomic base path and for truthful route
provenance value `:z_axis_diatomic_wl_base`. This is not an artifact schema
change.

Forbidden in this lane:

- driver public input changes or driver special cases;
- old White-Lindsey H1/H1+J materialization revival or adaptation;
- artifact schema changes, matrix-key changes, reader behavior changes, or
  manifest shape changes;
- Residual Gaussian, MWG/IDA, supplement, ECP, solver, or Cr2 workflow work;
- route diagnostics, stop-after controls, report/status/payload fields,
  raw-block switches, or route-stage labels;
- committed tests or committed driver fixtures.

Line budget: at most `250` added `src` lines, with deletion or simplification
of obsolete blocker-only WL diatomic guards expected where practical.

### 1a. White-Lindsey Diatomic Compact Retained Basis

Status: approved for implementation under `HP-WLDIAT-COMPACT-FN-01` and
`HP-WLDIAT-COMPACT-TEST-01`.

The current WL diatomic artifact path can be mechanically realized, but the
retained-basis shape is still an elongated shared-shell boundary-stratum
identity realization. That is not the intended compact WL retained basis.

Approved source files:

```text
src/cartesian_shellification/terminal_geometry.jl
src/cartesian_terminal_lowering/region_contracts.jl
src/cartesian_retained_units/lower_contract_units.jl
src/cartesian_retained_unit_transform_contracts/unit_contracts.jl
src/cartesian_final_basis_realization/white_lindsey_terminal_basis_realization.jl
src/pqs_source_box_route_driver_helpers.jl
```

The correction must preserve the WL faces/edges/corners and small boundary
unit model, but the terminal units must carry or realize compact retained
columns from products of one-dimensional contractions. It must not force a
persistent shell object after splitting, retain full-support identity rows as
the production compact basis, fake compactness by dropping rows, change the
driver, revive old WL H1/H1+J materialization, change artifacts, change PQS
behavior, touch raw blocks/RG/MWG/IDA, add diagnostics/status payloads, add
committed fixtures/tests, or run Cr2.

### 1b. Base Homonuclear Z-Axis Diatomic Validation

Status: approved for implementation under `HP-COMP-BASEDIAT-FN-01` and
`HP-COMP-BASEDIAT-TEST-01`.

Goal: relax the base producer's two-center validation from H2-only to explicit
homonuclear z-axis all-electron diatomics, using the same public
geometry/electron-count contract already used by the supplemented diatomic
facade.

Approved source file:

```text
src/cartesian_base_hamiltonian.jl
```

This lane keeps the basis contract unchanged and preserves both `nesting = :pqs`
and `nesting = :wl`, but it does not approve route/shellification/terminal
lowering changes. The WL path still depends on the separately approved
`HP-COMP-WLDIAT-*` terminal-record lane.

Forbidden in this lane:

- driver changes;
- supplement, Residual Gaussian, MWG/IDA, or artifact schema changes;
- route skeleton, shellification, terminal lowering, raw-block, reader,
  public API/export, solver/ECP, diagnostic/status/report, or Cr2-specific
  workflow changes;
- element lookup/default tables or inferred electron counts;
- heteronuclear, translated, non-z-axis, or general-geometry support.

Line budget: target under `60` added `src` lines.

### 2. Supplemented Atoms

Status: approved for implementation under `HP-COMP-SUPPATOM-FN-01` and
`HP-COMP-SUPPATOM-TEST-01`.

Goal: make `geometry = atom`, `supplement = on`, and either `nesting = :pqs`
or `nesting = :wl` work through the same Residual Gaussian path used by
supplemented diatomics. One-center residual selection is the one-owner case of
the same owner-local residual Gaussian algorithm.

Approved source files:

```text
src/cartesian_base_hamiltonian.jl
bin/cartesian_ham_builder.jl
src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl
```

The terminal residual file is optional and may be touched only if the existing
RG/MWG compatibility entry point exposes a direct one-owner genericity blocker.
The expected source work is to use `legacy_atomic_gaussian_supplement(...)` for
one-center supplement loading, keep the existing diatomic supplement loader for
two-center inputs, and relax only the canonical driver's supplemented
`Natom == 2` guard.

This lane must preserve the base atom validation, terminal basis construction,
residual Gaussian augmentation, exact augmented operators, residual MWG/IDA
interaction, base K/U reuse, assembly, writer, readback, manifest/provenance,
driver public inputs, hooks, spacing/layout, stage labels, and artifact
contract. It must not introduce a separate atom supplement algorithm,
atom-specific Hamiltonian builder, atom-only materialization path, atom-only
artifact schema, new driver inputs, route switches, diagnostics, or stop-after
controls.

### 3. Supplemented White-Lindsey

Status: approved for implementation under `HP-COMP-SUPPWL-FN-01` and
`HP-COMP-SUPPWL-TEST-01`.

Goal: allow `geometry = z-axis diatomic`, `supplement = on`, and
`nesting = :wl` to use the same supplemented homonuclear z-axis diatomic
facade/staged path as `nesting = :pqs`.

Approved source files:

```text
src/cartesian_base_hamiltonian.jl
src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl
```

The terminal residual file is optional and may be touched only if the existing
RG/MWG compatibility entry point exposes a direct genericity blocker. The
expected source work is removal of the early supplemented-WL blockers in
`src/cartesian_base_hamiltonian.jl` after proving the existing Residual
Gaussian/MWG path consumes the WL `CartesianTerminalBasisRealization`.

This lane must preserve the existing supplement contract, residual selection,
exact augmented operators, residual MWG/IDA interaction, base K/U reuse,
artifact keys, manifest/provenance, driver inputs, and stage labels. It should
not branch the driver on WL supplemented cases and should not add diagnostics
or route-stage switches.

## Composition IDs

Approved composition IDs authorize only the exact lanes named above.

- `HP-COMP-WLDIAT-FN-01` / `HP-COMP-WLDIAT-TEST-01`: approved WL z-axis
  diatomic base terminal-basis and artifact path.
- `HP-COMP-BASEDIAT-FN-01` / `HP-COMP-BASEDIAT-TEST-01`: approved base
  homonuclear z-axis diatomic validation relaxation.
- `HP-COMP-SUPPATOM-FN-01` / `HP-COMP-SUPPATOM-TEST-01`: supplemented
  one-center atom path through common Residual Gaussian augmentation.
- `HP-COMP-SUPPWL-FN-01` / `HP-COMP-SUPPWL-TEST-01`: supplemented
  White-Lindsey z-axis diatomic path through the common RG boundary.
- `HP-COMP-WLNS-FN-01` / `HP-COMP-WLNS-TEST-01`: WL z-axis diatomic `ns`
  early rejection and retained-support saturation wording.
- `HP-WLDIAT-COMPACT-FN-01` / `HP-WLDIAT-COMPACT-TEST-01`: WL z-axis
  diatomic compact retained-basis correction.

The initial explicit `atom | z-axis diatomic`, `:pqs | :wl`,
`supplement = off | on` composition lanes now all have approved implementation
authority, with WL diatomic compact retained-basis correction tracked
separately under `HP-WLDIAT-COMPACT-*`. Future geometry, physics, export,
solver, or Cr2-specific expansion still needs a separate docs-only amendment.

## Deferred

- translated atoms;
- non-z-axis diatomics and rotations;
- heteronuclear supplemented workflow beyond explicit future approval;
- ECP, EGOI, solver/RHF, or CR2-specific branches;
- public export/API redesign;
- new Hamiltonian wrapper or artifact matrix format;
- route diagnostics, private stop-after stages, report/status payloads, and
  broad driver switches.
