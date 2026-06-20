# Cartesian PQS and IDA Overview

This page is the entry point for the Cartesian PQS and integral-diagonal
approximation (IDA) algorithm pages. It fixes common notation and the dependency
graph used by the public algorithm suite.

## Dependency Graph

The intended construction order is:

```text
1D gausslet / PGDG data
-> Cartesian product operators
-> PQS localized shell basis
-> localized IDA Hamiltonian
-> optional residual-Gaussian extension
-> optional deliberate orbital transformations
```

The basis-construction steps are separate from later orbital rotations. A
later frozen-core or natural-orbital transformation may rotate orbitals on
purpose, but it is not a repair step for basis construction.

## Notation

- `S`: overlap matrix.
- `K`: kinetic-energy matrix.
- `U_A`: uncharged attraction matrix for center `A`,
  `<chi_i|-1/r_A|chi_j>`.
- `Z_A`: nuclear charge at center `A`.
- `Vee`: two-index IDA electron-electron interaction matrix in the localized
  IDA basis.
- `C`: orbital coefficient matrix in a localized IDA basis.
- `L_s`: shell-local Lowdin transform for projected PQS shell `s`.
- `L_R`: residual-Gaussian Lowdin transform after projection out of the fixed
  gausslet/PQS space.

`T` is reserved for kinetic energy only. Orbital coefficient matrices are
called `C`, not `T`.

## Spaces and Dimensions

The base all-electron Cartesian IDA Hamiltonian has one localized working
basis. In that basis:

- `K` is `n x n`.
- each `U_A` is `n x n`.
- `Vee` is `n x n`.
- nuclear positions are stored as `ncenter x 3`, one row per center.
- orbital coefficient matrices `C` are `n x m` when a later calculation chooses
  `m` orbitals.

There is no public all-electron density-transform map in this base contract.
Under a deliberate orbital rotation, IDA contractions use elementwise products
of orbital coefficients in the same localized IDA basis.

## Algorithm Pages

- [Cartesian low-dimensional operator assembly](cartesian_low_dimensional_operator_assembly.md)
  defines how one-dimensional data build Cartesian overlap, kinetic, nuclear,
  moment, and IDA pair quantities.
- [PQS shell construction](pqs_shell_construction.md) defines the localized
  shell-projection algorithm, the bond-aligned diatomic atom-contact core
  hull rule, and the allowed shell-local Lowdin transforms.
- [Residual-Gaussian extension](residual_gaussian_extension.md) defines the
  one-sided residualization and residual-block Lowdin for Gaussian supplements.
- [IDA Hamiltonian and counterpoise](ida_hamiltonian_and_counterpoise.md)
  defines the public base Hamiltonian convention: `K`, separated unit-nuclear
  `{U_A}`, `Vee`, charges, positions, and nuclear repulsion.

## Allowed Orthogonalizations

Only two orthogonalization steps are part of the PQS plus residual-Gaussian
algorithm:

1. shell-local Lowdin transforms inside each projected PQS shell block;
2. residual-Gaussian Lowdin transforms inside the Gaussian residual block after
   projection out of the fixed gausslet/PQS space.

## Forbidden Operations

- No symmetric Lowdin over all assembled PQS core and shell functions.
- No symmetric Lowdin over gausslets/PQS plus residual Gaussians.
- No rotation of previously accepted shell blocks while constructing a later
  shell.
- No use of a global orthogonalization to hide a shell-projection error.
- No public all-electron split between a "final orbital basis" and a separate
  "density gauge" for IDA. The base all-electron IDA contract has one localized
  working basis.

## Numerical Invariants

The completed localized basis should satisfy:

```math
B^T S B \approx I.
```

For residual Gaussians appended to a fixed PQS/gausslet basis `G`, the expected
checks are:

```math
G^T S G \approx I,
\qquad
G^T S R \approx 0,
\qquad
R^T S R \approx I.
```

## Code Map

Current implementation surfaces are:

- `src/cartesian_final_basis_realization/pqs_source_shell_final_basis.jl`
  for shell-local projection and Lowdin realization.
- `src/cartesian_final_basis_realization/pqs_complete_core_shell_final_basis.jl`
  for common complete core/shell basis and IDA interaction helpers.
- `src/pqs_multilayer_complete_core_shell_h1.jl` for the common PQS H1 and
  density-interaction payload.
- `src/cartesian_ida_hamiltonian.jl` for the public one-basis IDA Hamiltonian
  object and minimal artifact reader/writer.
- `src/cartesian_gaussian_axis_integrals.jl` for shared Cartesian Gaussian axis
  integral kernels where present.

## Current Implementation Deviations

The active H2 PQS route no longer applies the forbidden global core/shell
Lowdin cleanup. The public `CartesianIDAHamiltonian` type and minimal
versioned writer/reader now exist. The active H2 residual-GTO route is still
H2/q5-specific producer work rather than a general Cr2-ready diatomic
constructor.
