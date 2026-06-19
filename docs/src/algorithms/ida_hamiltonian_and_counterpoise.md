# IDA Hamiltonian and Counterpoise

This page defines the public base Hamiltonian convention for Cartesian IDA
routes. The base all-electron Hamiltonian has one localized IDA working basis.

## Spaces and Dimensions

For a localized IDA basis with `n` functions:

- `K` is `n x n`.
- each unit-nuclear attraction matrix `U_A` is `n x n`.
- `Vee` is `n x n`.
- nuclear charges have length `ncenter`.
- nuclear positions are `ncenter x 3`, one row per center.

## Inputs

- localized IDA basis;
- kinetic matrix `K`;
- separated uncharged attraction matrices `{U_A}`;
- two-index IDA electron-electron matrix `Vee`;
- nuclear charges `{Z_A}`;
- nuclear positions `{R_A}`;
- spin counts `nup`, `ndn`.

## Outputs

- base Hamiltonian data:
  - `K`;
  - `{U_A}`;
  - `Vee`;
  - `nup`, `ndn`;
  - nuclear charges and positions;
  - nuclear repulsion `E_nn`;
- caller-assembled one-body matrices for normal and counterpoise branches.

## Pseudocode

1. Construct the localized IDA basis by the PQS shell and optional
   residual-Gaussian algorithms.

2. Build kinetic in that basis:

   ```math
   K_{ij} = \langle \chi_i | -\tfrac{1}{2}\nabla^2 | \chi_j \rangle.
   ```

3. Build each uncharged nuclear attraction matrix separately:

   ```math
   (U_A)_{ij} = \langle \chi_i | -1/r_A | \chi_j \rangle.
   ```

4. Build the two-index IDA electron-electron matrix `Vee` in the same
   localized basis.

5. Store nuclear charges and positions separately.

6. Let consumers assemble the charged one-body Hamiltonian:

   ```math
   H_1 = K + \sum_A Z_A U_A.
   ```

7. For monomer or counterpoise branches, keep the same basis and `Vee`, but
   choose the electron/spin counts and center weights `w_A`. A branch one-body
   matrix is assembled as:

   ```math
   H_1(w)=K+\sum_A w_A Z_A U_A.
   ```

   The branch nuclear repulsion uses the same weights:

   ```math
   E_{nn}(w)=\sum_{A<B}\frac{(w_AZ_A)(w_BZ_B)}{|R_A-R_B|}.
   ```

## Linear Algebra

The IDA electron-electron approximation in the localized basis is:

```math
(ij|kl)_{\rm IDA} \approx \delta_{ij}\delta_{kl} Vee_{ik}.
```

Under a deliberate orbital rotation `C`, the corresponding orbital-basis
contraction is:

```math
(pq|rs)_{\rm IDA}
=
(C_{\cdot p}\odot C_{\cdot q})^T
Vee
(C_{\cdot r}\odot C_{\cdot s}).
```

This is the intended IDA model. It is not a placeholder for dense four-index
ERI storage.

## Allowed Orthogonalizations

None at Hamiltonian assembly. Orthogonalization belongs to shell construction
and residual-Gaussian construction. Later orbital transformations are deliberate
consumer transformations, not Hamiltonian-construction repair.

## Forbidden Operations

- Do not charge and sum nuclear attraction matrices before the public
  Hamiltonian boundary if counterpoise is needed.
- Do not store only `H_1` when the consumer needs `{U_A}`.
- Do not introduce a public all-electron density-transform map for the base
  localized IDA basis.
- Do not replace IDA with dense four-index ERIs as the public Cartesian route.
- Do not use a global Lowdin cleanup as part of Hamiltonian assembly.

## Numerical Invariants

- `K`, each `U_A`, and `Vee` are finite and symmetric.
- `E_nn` is finite for nonzero internuclear separations.
- `H_1 = K + sum_A Z_A U_A` reproduces the normal charged one-body matrix.
- Counterpoise branches use the same basis and `Vee` as the full system while
  allowing branch-specific electron/spin counts and center weights.

## Operator and Gauge Conventions

For a two-center system `A`, `B`:

```math
H_1^{AB}=K+Z_AU_A+Z_BU_B,
```

```math
H_1^{A+\mathrm{ghost},B}=K+Z_AU_A,
```

```math
H_1^{\mathrm{ghost},A+B}=K+Z_BU_B.
```

The basis and `Vee` remain fixed for these branches.

Equivalently, these are center-weight choices `w=(1,1)`, `w=(1,0)`,
and `w=(0,1)` using the stored physical charges.

## Public API

`CartesianIDAHamiltonian` is the public in-memory object for this all-electron
one-basis contract. Use `one_body_hamiltonian(ham; center_weights=...)` and
`nuclear_repulsion(ham; center_weights=...)` for full and counterpoise
branches.

`write_cartesian_ida_hamiltonian` and `read_cartesian_ida_hamiltonian` provide
the minimal versioned JLD2 artifact. The artifact stores only `K`, `{U_A}` as
an `n x n x ncenter` tensor, `Vee`, charges, `ncenter x 3` positions, and spin
counts. It does not store nuclear repulsion or route diagnostics.

## Code Map

- `src/pqs_multilayer_complete_core_shell_h1.jl` builds the common complete
  core/shell H1 path and separated center contributions.
- `src/cartesian_final_basis_realization/pqs_complete_core_shell_final_basis.jl`
  builds IDA density interaction helpers in the completed localized basis.
- `src/cartesian_ida_hamiltonian.jl` contains the public one-basis IDA
  Hamiltonian object and minimal artifact reader/writer.
- `src/pqs_h2_residual_gto_handoff.jl` contains the current private H2
  residual-GTO producer/readback scaffolding.
- Future public code should expose a one-basis IDA Hamiltonian object rather
  than the current private H2 sidecar artifact.

## Current Implementation Deviations

The public `CartesianIDAHamiltonian` type and minimal versioned writer/reader
exist. The private H2 residual-GTO route still writes its older sidecar
artifact until the next replacement/deletion pass.
