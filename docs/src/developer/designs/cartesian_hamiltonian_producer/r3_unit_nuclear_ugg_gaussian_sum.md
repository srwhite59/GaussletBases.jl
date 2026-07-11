# R3 Unit-Nuclear U_GG Gaussian Sum

Status: implemented internal optimization and completed validation contract.

This page is the canonical contract for:

- `HP-R3UN-FN-01`;
- `HP-R3UN-TEST-01`.

The implementation is owned by `CartesianFinalBasisRealization`. It constructs
the exact terminal final-basis, uncharged nuclear-attraction block for each
nuclear center. Neutral Gaussian `G-A` and `A-A` nuclear blocks remain owned by
[Cartesian Gaussian raw blocks](cartesian_gaussian_raw_blocks_nuclear.md).

## Source And Entry Points

Implemented source:

```text
src/cartesian_final_basis_realization/pqs_terminal_one_body.jl
src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl
```

The principal helpers are `_accumulate_terminal_gaussian_sum!` and its
function-local workspace form. Supporting helpers validate factor families,
fill one term-first action, and add terminal blocks.
`pqs_terminal_residual_gto_augmented_unit_nuclear(...)` supplies center-specific
factor families and transforms the resulting `G-G`, `G-A`, and `A-A` blocks
into the augmented basis.

Commit `79e5cd474` implemented the allocation reduction. Manager-log Passes
089, 089A, and 089B preserve the accepted parity and allocation evidence.

## Numerical Contract

For center `A`, the terminal kernel assembles the Gaussian expansion

```text
U_GG[A] = -sum_t c_t tensor(Fx[t], Fy[t], Fz[t])
```

blockwise in the realized terminal basis. The minus sign belongs to the
uncharged attractive unit operator `U_A = -1/r_A`. Physical nuclear charge is
applied later by Hamiltonian assembly and must not be folded into this block.

The implementation is term-first and respects direct identity and compact
support-local terminal blocks. It evaluates one triangle and mirrors
off-diagonal blocks. Before accumulation it requires:

- a square destination of size `basis.final_dimension`;
- equal coefficient and x/y/z factor-family term counts;
- factor matrices large enough for all referenced support states;
- finite factor matrices symmetric within `1e-10` in infinity norm.

When no trusted base blocks are supplied, the augmented path computes one
matrix per center. It reuses one function-local set of action, tile, and block
workspaces across centers. A center already represented by a PGDG axis may
reuse that axis's factor family; translated centers are evaluated through the
same mapped ordinary one-body owner and the carried Coulomb expansion.

The producer-wide expansion and PGDG exponent sequence must agree before this
path runs. The number of nuclear locations must match the number of physical
charges, although charges are not applied inside `U_GG`.

Trusted same-construction `U_GG[A]` reuse and its live recomputation fallback
belong to [same-construction base reuse](r3_same_construction_base_reuse.md).

## Validation

The maintained focused gate is:

```text
test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl
```

It exercises the fallback path, exact augmented nuclear blocks, finiteness,
symmetry, base-block parity, and the H2 endpoint. Accepted Be2 and ignored Cr2
replay/allocation evidence remains in manager-log Passes 089-089B.

`HP-R3UN-TEST-01` is completed validation and maintenance permission. It does
not authorize a separate Gaussian-sum framework or development-only fixture.

## Boundaries

This contract does not own or change:

- neutral nuclear `G-A` or `A-A` kernels;
- terminal kinetic or moment products;
- residual selection, orientation, or exact transforms;
- IDA, MWG, routes, parent construction, or terminal realization;
- physical charge application;
- persistent caches or workspace objects;
- metadata, reports, artifacts, public APIs, solvers, or Cr2 workflows.

Workspace reuse is function-local. No persistent factor cache or stage object
is part of the implemented contract.
