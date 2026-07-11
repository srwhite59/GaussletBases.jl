# R3 Residual-GTO/MWG Compatibility History

Status: historical compatibility index. R3-A, R3-B, and R3-C are implemented
milestones, but this page is not current numerical authority. Current residual
Gaussian algorithms are canonical in
[Residual Gaussian domain module](residual_gaussian_domain_module.md), the
supported internal facade is canonical in
[R3 usability supplemented workflow](r3_usability_supplemented_workflow.md),
and artifact behavior is canonical in
[Cartesian Hamiltonian artifact manifest](cartesian_hamiltonian_artifact_manifest.md).

## Purpose

R3 established that a Gaussian supplement can augment a Cartesian terminal
basis without creating a second Hamiltonian architecture:

```text
same base construction
+ residual supplement directions
-> exact augmented one-body operators
+ residual-containing MWG interaction
-> CartesianIDAHamiltonian
-> optional existing artifact with compact supplement provenance
```

The `R3-A`, `R3-B`, and `R3-C` labels remain useful for Git and manager-log
archaeology. They are not production module, object, or algorithm names.

## Final Dispositions

| Milestone | Durable result | Current owner |
| --- | --- | --- |
| R3-A | Residual basis plus exact augmented kinetic, by-center unit-nuclear, coordinate, and second-moment matrices | `CartesianResidualGaussians`; [domain contract](residual_gaussian_domain_module.md) |
| R3-B | Moment-matched residual descriptors, residual-containing MWG/IDA blocks, and an in-memory `CartesianIDAHamiltonian{Float64}` | `CartesianResidualGaussians`, with terminal compatibility assembly |
| R3-C | Compact `supplement_provenance/` attached to the existing Hamiltonian artifact | terminal/facade workflow; [artifact manifest](cartesian_hamiltonian_artifact_manifest.md) |

The compatibility IDs are:

- `HP-R3-OBJ-01`;
- `HP-R3-FN-01`;
- `HP-R3-FN-02`;
- `HP-R3-FN-03`;
- `HP-R3-ART-01`;
- `HP-R3-TEST-01`.

Their lifecycle, exact source/test surfaces, dependencies, and permissions are
recorded in [the registry](registry.md). This page records why those surfaces
exist and where their current law moved.

## Current Ownership

### Residual Numerical Authority

The current numerical implementation is in:

```text
src/cartesian_residual_gaussians/residual_basis.jl
src/cartesian_residual_gaussians/augmented_operators.jl
src/cartesian_residual_gaussians/mwg_interaction.jl
```

Those files own owner-local residual selection, the final inter-owner merge,
exact augmented-operator transformation, moment-matched residual descriptors,
and residual-containing IDA interaction assembly. Their formulas, tolerances,
failure behavior, and validation gates live only in the domain contract.

### Compatibility And Artifact Boundary

The historical terminal owner remains:

```text
src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl
```

It now provides construction checks, compatibility entry points, raw-block
composition, `CartesianIDAHamiltonian` assembly, and the compact supplemented
artifact writer. In particular,
`CartesianTerminalResidualGTOAugmentation` is an alias of the domain-owned
`CartesianResidualGaussianBasis`; it is not a second object contract.

The live compatibility assembly is
`pqs_terminal_residual_gto_augmented_hamiltonian(...)`. It delegates residual
interaction physics to `CartesianResidualGaussians` and returns the existing
Hamiltonian type directly.

### Supported Workflow

`src/cartesian_base_hamiltonian.jl` owns the non-exported
`cartesian_residual_gto_mwg_hamiltonian(...)` facade. It constructs the base
Hamiltonian, terminal basis, parent bundles, Gaussian supplement, residual
basis, exact augmented operators, and MWG interaction from one construction.
See the usability contract for its supported inputs and failure behavior.

## Durable Lessons

### Owner-Local Selection

The first global raw-candidate selection was superseded. Residual directions
are selected separately by physical owner and merged once afterward. Global
candidate Lowdin or global raw-column pivoted-Cholesky selection must not be
revived through an R3 compatibility helper. The numerical-complete `1e-10`
lane is a separate explicit policy, not a reinterpretation of ordinary R3.

### Exact One-Body Versus MWG

Kinetic, every uncharged by-center nuclear matrix, coordinates, and second
moments are exact transformations of the raw `[G,A]` blocks into `[G,R]`.
Only residual-containing electron-electron blocks use the separable
moment-matched Gaussian approximation. Treating the exact one-body transform
as MWG, or treating MWG as exact residual-GTO Coulomb, is a category error.

### Same Construction

The base Hamiltonian, terminal basis, parent bundles, supplement, and carried
Coulomb expansion must belong to one construction. Dimension-compatible
objects assembled independently are not interchangeable. The current code
checks residual and reused-matrix dimensions, center/charge counts, and PGDG
exponent parity before composing the augmented interaction. Trusted kinetic and
unit-nuclear block reuse, including its exact fallback, is canonical in
[R3 same-construction base reuse](r3_same_construction_base_reuse.md).

### Compact Provenance

R3-C deliberately reused the existing Hamiltonian artifact and added compact
supplement provenance. It did not authorize serialization of `T_G`, `T_A`,
dense moments, full residual spectra, matched-Gaussian arrays, or a new result
wrapper. Current artifact keys and compatibility behavior belong to the
artifact-manifest contract, not this history page.

### Rejected Architecture

The following remain rejected lessons, not alternative workflows:

- post-hoc augmentation of an opaque Hamiltonian;
- global raw-candidate residual selection;
- route-provider or supplement-preflight payloads;
- status/result/report wrappers around the Hamiltonian;
- exposing terminal objects, residual objects, pair factors, or provenance
  payloads through the supported facade;
- direct interaction rotation in place of the inherited IDA/MWG convention.

## Evidence Index

Key implementation and correction commits:

- `08d9a92eb`, `50a9039f6`, `26e31c1bc`: R3-A basis, exact operators, and
  standalone endpoint;
- `cbf22c20c`, `0846d8e37`, `c36581e3e`: R3-B interaction, normalization
  correction, and tracked endpoint;
- `bd1c2d37d`: R3-C supplemented provenance writer;
- `923edf1f7`: supported internal usability facade;
- `3e9af7117`: owner-local residual-selection correction;
- `f4c26804d`, `7c0651be5`, `7e45d5da1`, `8d9b501f1`: domain ownership and
  removal of obsolete R3-B wrappers.

Manager-log Passes 049-071 preserve the detailed implementation chronology,
including superseded fixture values and performance measurements. Those
numbers are evidence only. The current tracked endpoint is the committed test,
not the old `489`-dimension narrative retained in Git history.

## Current Validation Surface

The focused compatibility and facade gate is:

```text
test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl
```

It exercises residual geometry, exact augmented operators, independent
weight-aware interaction checks, the in-memory Hamiltonian, facade input
rejection, artifact write/readback, and compact provenance. It is standalone
and is not normal `Pkg.test` pressure.

At the Pass 377 baseline, its facade fixture has base dimension `487`, residual
dimension `18`, and augmented dimension `505`. Historical `489` values must
not be reintroduced as current targets.

## Boundaries

This compatibility record does not own or reopen:

- residual cutoffs, numerical-complete policy, injection, or localization;
- protected-localized, screened-Hartree, EGOI, or solver behavior;
- artifact-manifest extensions;
- public exports or general molecular geometry;
- driver defaults or Cr2-specific workflow.

Use the linked current contracts for those questions. Git history and the
manager running log preserve the removed donor inventories, pseudocode, line
budgets, probe chronology, and superseded implementation plans.
