# White-Lindsey Low-Order Density-Density Builder Contract

## Purpose

This note defines the missing private builder needed before the Cartesian
nesting route driver can write a real White-Lindsey low-order Ham bundle.

Contract label:

```julia
:white_lindsey_low_order_fixed_block_density_density_builder
```

This is a private development contract. It is not a public API, not a new Ham
schema, and not an adoption of supplement-required QW paths as White-Lindsey
benchmark semantics.

## Current Input

The builder should consume one of:

- the materialized White-Lindsey low-order seed report from
  `_white_lindsey_low_order_materialized_seed_report()`;
- or the report's materialized `_NestedFixedBlock3D` directly.

The current seed fixed block represents the final retained low-order basis as
disjoint retained pieces:

- direct core sites;
- face interiors as two-dimensional products of one-dimensional retained side
  functions;
- edges as one-dimensional retained side functions;
- corners as direct single sites.

The final retained basis is the fixed block's retained basis. The builder must
not reinterpret shell/support rows as a new algorithm and must not introduce a
Gaussian supplement as the White-Lindsey benchmark route.

## Available Data

The fixed block already supplies the one-body and basis facts needed by the
route-driver materialization report:

- `fixed_block.overlap`;
- `fixed_block.position_x`, `position_y`, and `position_z`;
- `fixed_block.x2_x`, `x2_y`, and `x2_z`;
- `fixed_block.kinetic`;
- `fixed_block.weights`.

`fixed_block.weights` are final IDA weights for this route: unsquared
integrals of the final retained basis functions. They are not positive
quadrature weights for arbitrary retained-column integration, and they are not
PQS retained-column diagnostic weights.

## Required Missing Output

The missing builder must produce a real retained-space two-index
density-density interaction matrix:

```julia
interaction_matrix::Matrix{Float64}
```

The matrix must be compatible with the existing Cartesian Ham bundle
conventions:

- square size equal to the final retained dimension;
- finite entries;
- symmetric to the tolerance chosen by the builder contract;
- density-density / IDA semantics, not a four-index Galerkin Coulomb tensor;
- paired with the same final retained basis and final IDA weights carried by
  the fixed block.

The eventual full payload target is an existing writer-compatible object or
private conversion payload containing:

- overlap;
- one-body Hamiltonian;
- interaction matrix;
- final retained-basis integral weights;
- labels and centers needed by the bundle writer;
- route metadata identifying the private White-Lindsey low-order benchmark
  path.

`cartesian_basis_bundle_payload` and `write_cartesian_basis_bundle_jld2` should
remain the preferred export contract once a real operator payload exists.

## Unsettled Decisions

The next implementation pass must not silently choose these details:

- the low-order density-density contraction formula;
- scaling and normalization ownership for the retained density-density matrix;
- whether the builder returns `OrdinaryCartesianOperators3D`, a new private
  payload, or a conversion object accepted by the existing writer;
- how to avoid dense parent four-index work at useful sizes;
- validation fixtures beyond the current tiny one-center seed;
- performance targets and allocation budget for a representative retained
  dimension.

Until those decisions are explicit, Ham export should remain blocked with:

```julia
:missing_pure_low_order_fixed_block_density_density_interaction_builder
```

## Reuse Candidates

Useful existing pieces include:

- `src/white_lindsey_materialized_seed.jl` for the fixed block, retained
  pieces, final weights, and one-body inventory;
- `cartesian_basis_bundle_payload` and `write_cartesian_basis_bundle_jld2` for
  the basis/Ham export shape;
- `OrdinaryCartesianOperators3D` and `OrdinaryCartesianIDAOperators` for
  existing payload conventions;
- ordinary Cartesian QW/IDA density-density helpers as convention references;
- existing Coulomb Gaussian expansion utilities;
- old low-order/nesting implementation material, if it can be identified
  locally without importing unreviewed route semantics.

Reuse candidates are not adoption decisions. In particular, supplement-required
atomic paths remain diagnostic inventory only for this benchmark route unless
a later manager-scoped pass changes that boundary explicitly.
