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

## Design Inventory: 2026-06-04

The strongest current reuse candidate is the existing nested fixed-block
pair-sum path:

- `src/cartesian_nested_faces.jl`
  - `_nested_factorized_weight_aware_pair_terms(...)`;
  - `_nested_weight_aware_pair_terms(...)`;
  - `_nested_shell_packet(...)`;
- `src/ordinary_qw_raw_blocks.jl`
  - `_qwrg_fixed_block_interaction_matrix(fixed_block, expansion)`.

The current White-Lindsey seed fixed block already carries `pair_sum`. A
focused probe on the tiny `7`/`5` seed found:

- `fixed_block.pair_sum !== nothing`;
- `_qwrg_fixed_block_interaction_matrix(fixed_block, expansion)` returns a
  finite symmetric `223 x 223` matrix;
- the symmetry error was `0.0` for the probed seed.

This suggests the likely builder formula is not a new Coulomb tensor
construction. The likely first builder should consume the retained fixed
block's existing density-density pair sum:

```julia
interaction_matrix =
    _qwrg_fixed_block_interaction_matrix(fixed_block, expansion)
```

Conceptually, `fixed_block.pair_sum` is the retained two-index
density-density matrix assembled during nested packet construction. It is the
Coulomb-Gaussian expansion sum of density-normalized retained pair factors,
not a four-index Galerkin Coulomb tensor and not a placeholder quadrature
matrix.

### Scaling And Normalization Ownership

Current ownership appears to be:

- Coulomb Gaussian exponents and coefficients are selected at sequence build
  time. For the one-center atomic seed,
  `_one_center_atomic_term_coefficients(...)` verifies that the PGDG bundle
  exponents match the supplied `CoulombGaussianExpansion`, then uses the
  expansion coefficients as `term_coefficients`.
- Raw pair-factor normalization is owned by the nested packet pair-term
  helpers, not by the eventual Ham writer.
- In the factorized packet path,
  `_nested_factorized_weight_aware_pair_terms(...)` projects axis integral
  weights, builds normalized axis pair-term tables, and contracts the
  coefficient-weighted three-axis term products.
- In the support-reference packet path,
  `_nested_weight_aware_pair_terms(...)` computes final retained weights,
  divides support coefficients by those weights, and contracts raw pair terms.
  This path is a small-fixture authority/debug comparison, not the intended
  scalable algorithm.
- `fixed_block.weights` are the final retained-basis integral weights. The
  builder should carry them through to the payload and must not divide by them
  again when consuming `fixed_block.pair_sum`.
- `_qwrg_fixed_block_interaction_matrix(...)` owns the final symmetrization of
  the stored pair sum.

The remaining scaling decision is whether the future builder should trust
`fixed_block.pair_sum` as already normalized and coefficient-summed, or expose
a more explicit receipt proving those facts before a Ham bundle can be written.

### Proposed Smallest Validation Fixture

The smallest useful fixture is the current one-center White-Lindsey seed:

```julia
_white_lindsey_low_order_materialized_seed_fixture(
    packet_kernel = :factorized_direct,
)
```

Validation should compare it against the same seed with:

```julia
packet_kernel = :support_reference
```

The first design probe found:

- factorized interaction size: `(223, 223)`;
- support-reference interaction size: `(223, 223)`;
- both outputs finite;
- both symmetry errors `0.0`;
- `norm(factorized - support_reference, Inf) =
  2.5757174171303632e-14`.

This fixture should check:

- `fixed_block.pair_sum !== nothing`;
- `length(fixed_block.weights) == retained_dimension`;
- all final weights are finite; for this seed, `minimum(fixed_block.weights) >
  0.0` is a useful regression check, but positivity must not be promoted to a
  generic quadrature-weight semantic;
- interaction matrix size is `(retained_dimension, retained_dimension)`;
- finite output;
- symmetry error;
- factorized-direct versus support-reference agreement.

Dense parent four-index construction is not needed as the intended algorithm.
It may be considered only as an additional tiny debug authority if a future
manager-scoped pass asks for it.

### Remaining Ambiguities

The next manager/user decisions are:

- whether `_qwrg_fixed_block_interaction_matrix(...)` is acceptable as the
  interaction source for the private benchmark Ham builder, or whether a new
  wrapper must first audit `pair_sum` provenance;
- whether the full payload should be an `OrdinaryCartesianOperators3D`
  instance, a new private fixed-block Ham payload, or a conversion object
  consumed by the existing bundle writer;
- how labels and centers should be populated for the current one-center
  fixed-block retained pieces;
- whether the one-body Hamiltonian should initially use only kinetic plus
  stored `gaussian_sum` nuclear attraction, and how by-center nuclear sidecars
  should be represented for this private route;
- what retained dimension and packet kernel should define the first
  performance target beyond the tiny seed.
