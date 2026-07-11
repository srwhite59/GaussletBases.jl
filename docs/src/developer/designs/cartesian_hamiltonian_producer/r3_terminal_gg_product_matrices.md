# R3 Terminal G-G Product Matrices

Status: implemented internal optimization and completed validation contract.

This page is the canonical contract for:

- `HP-R3GG-FN-01`;
- `HP-R3GG-TEST-01`.

The implementation is owned by `CartesianFinalBasisRealization`. It constructs
the terminal final-basis `G-G` blocks needed by Residual Gaussian exact
augmented operators. Neutral Gaussian `G-A` and `A-A` blocks remain owned by
[Cartesian Gaussian raw blocks](cartesian_gaussian_raw_blocks_non_nuclear.md).

## Source And Entry Points

Implemented source:

```text
src/cartesian_final_basis_realization/pqs_terminal_one_body.jl
src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl
```

The reusable terminal kernel is `assemble_terminal_product_operator!`, with
the internal workspace form `_assemble_terminal_product_operator!`.
`pqs_terminal_residual_gto_augmented_products(...)` uses that kernel for:

- kinetic `K_GG`;
- first moments `x_GG`, `y_GG`, and `z_GG`;
- second moments `x2_GG`, `y2_GG`, and `z2_GG`.

Commits `fb9b0414a` and `5cd9e15a6` implemented destination and workspace
reuse. Manager-log Passes 087, 087A, and 087B preserve the accepted numerical
and allocation evidence.

## Numerical Contract

The kernel receives three one-dimensional factors and assembles their tensor
product blockwise in the realized terminal basis. It respects both terminal
representations:

- direct blocks have implicit identity coefficients on their owned support;
- compact blocks apply their support-local coefficient matrices.

Only one triangle is evaluated. Off-diagonal terminal blocks are mirrored into
the other triangle. The destination is accumulated into rather than silently
cleared, which permits the three Cartesian kinetic contributions to be added
to one matrix.

Before assembly, the kernel requires:

- a square destination of size `basis.final_dimension`;
- each one-dimensional factor to cover every referenced support state;
- finite factor entries;
- factor symmetry within `1e-10` in infinity norm.

The caller owns zeroing when a new product is required. In the augmented
operator path, one `G-G` scratch matrix and one function-local set of action,
tile, and block workspaces are reused. Each first or second moment is
transformed into the augmented basis immediately before the scratch matrix is
reused.

If a trusted base kinetic matrix is supplied, the augmented path may use it
instead of recomputing `K_GG`. That trust and fallback behavior belong to
[same-construction base reuse](r3_same_construction_base_reuse.md).

## Validation

The maintained focused gate is:

```text
test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl
```

It covers exact augmented operators, base `G-G` block parity, finiteness,
symmetry, and the H2 endpoint. Accepted Be2 and ignored Cr2 replay/allocation
evidence is retained in manager-log Passes 087-087B rather than repeated here.

`HP-R3GG-TEST-01` is completed validation and maintenance permission. It does
not authorize a separate product framework or development-only fixture.

## Boundaries

This contract does not own or change:

- neutral `G-A` or `A-A` raw blocks;
- unit-nuclear Gaussian sums;
- terminal basis realization;
- residual selection, orientation, or exact transforms;
- IDA or MWG interaction conventions;
- Qiu-White semantics, route construction, or parent construction;
- persistent caches or workspace objects;
- metadata, reports, artifacts, public APIs, solvers, or Cr2 workflows.

Workspace reuse is function-local. No persistent carrier or cache is part of
the implemented contract.
