# White-Lindsey Terminal Basis Realization

Status: implemented under `HP-WLTERM-*`, `HP-COMP-WLDIAT-*`,
`HP-WLDIAT-COMPACT-*`, and `HP-WLDIAT-PARITY-*`. The registry owns exact ID
lifecycle and file permission. This page owns the current numerical boundary
for realizing native White-Lindsey retained units as the terminal basis used
by the shared Hamiltonian producer.

## Input And Output Boundary

The existing `:white_lindsey_low_order` route consumes construction-native
terminal supports, retained units, transform contracts, and parent-axis
bundles. It returns the same downstream object as PQS:

```julia
CartesianTerminalBasisRealization
```

containing ordered `CartesianTerminalBasisBlock` entries on authoritative
owned support rows. White-Lindsey does not own a parallel Hamiltonian object,
route result, report payload, or materialization path.

Both construction families first consume the common geometry contract in
`common_terminal_shell_decomposition.md`. They differ only afterward:

- PQS realizes a complete shell from a full source box by boundary-mode
  selection, owned-row restriction, and shell-local Lowdin;
- White-Lindsey realizes the same common shell support as compact facet,
  edge, corner, and boundary-stratum products of one-dimensional
  contractions;
- midpoint, outer-mismatch, and angular-extension thin slabs use the common
  compact face-stack realization under both families;
- true direct/core regions remain identity sectors.

## Realization Contract

Each retained unit must match its transform contract's source CPBs and use the
transform path appropriate to its unit kind.

Direct units use identity coefficients only on their owned support. A
`:white_lindsey_boundary_stratum_retained_unit` is not an identity function:
it realizes compact coefficients on its owned stratum support. Facets use the
neutral terminal face-product helper, edges use the corresponding nested edge
product, and corners contribute their local corner product. Compact thin-slab
units reuse the neutral face-stack helper.

All non-direct coefficients are support-local. Realization must preserve the
native retained-unit and transform-contract order, and each block's
`column_range` follows that order. Blocks may not acquire rows from an earlier
terminal region.

For every compact block, the local Gram matrix must satisfy the terminal
identity check at `identity_atol`, whose default is `1.0e-8`. Owned supports
must be pairwise disjoint, support indices and states must agree, and the
final object reports zero structural cross-block overlap. Failure is an error;
the implementation must not repair it by dropping support rows, relabeling a
unit as identity, or applying a global Lowdin transform.

## Retention Semantics

Public `ns` remains the common comparison input. The WL route-local size is

```text
q = ns - 2
```

while the direct nucleus-centered core side remains the common public-`ns`
value

```text
isodd(ns) ? ns : ns + 1
```

Odd-side enforcement applies only to that direct nucleus-centered core. It
must not be inherited by boundary strata or their one-dimensional
contractions. Boundary products use the requested WL retained count with
`enforce_symmetric_odd = false`.

For the canonical one-layer WL boundary examples:

```text
ns = 4, q = 2:  4^3 - 2^3 = 56 boundary columns
ns = 5, q = 3:  5^3 - 3^3 = 98 boundary columns
```

The old `26`-column result for `ns = 4` is not valid retained-count policy.
Full-support identity retention for WL boundary strata is likewise invalid.

White-Lindsey z-axis diatomics reject normalized `ns < 4` before route
construction. At supported sizes, physical parent extent and contact geometry
may saturate retained support or produce a final dimension different from
PQS. That is not permission to alter common shell ownership, fake compactness,
or promise equal PQS/WL dimensions.

## Downstream Seam

After terminal realization, WL uses the shared final-basis product,
unit-nuclear, one-body, IDA, residual-GTO/MWG, Hamiltonian, artifact, and
driver machinery. The bounded terminal inventory and due-diligence report
consume the same native support, retained, transform, and block facts. They do
not define WL coefficients or retained counts.

The removed route-global White-Lindsey H1/H1+J stack remains historical only.
It must not be revived through adapters, reports, route stages, or alternate
Hamiltonian construction.

## Source Ownership

Current ownership is limited to:

- `src/cartesian_final_basis_realization/white_lindsey_terminal_basis_realization.jl`
  for WL block realization and validation;
- `src/cartesian_final_basis_realization/terminal_face_product_blocks.jl` for
  the neutral facet and thin-slab coefficient primitive;
- `src/cartesian_final_basis_realization/CartesianFinalBasisRealization.jl`
  for module wiring;
- the established terminal-lowering, retained-unit, and transform-contract
  owners for native WL strata and compact slabs;
- `src/pqs_source_box_route_driver_helpers.jl` for narrow route-to-terminal
  wiring.

Common shell geometry belongs to its separate contract. The registry remains
authoritative for exact source surfaces.

## Guardrails

This contract does not change route skeletons, shellification ownership,
public `ns` normalization, direct-core centering, PQS source dimensions,
complete-shell or slab retained policies, mapped-COMX defaults, artifacts,
public APIs, reports, raw blocks, RG/MWG/IDA semantics, solvers, ECP, or Cr2
workflow. It does not choose a longitudinal `L` or an angular-resolution
scale.

If a WL unit cannot be realized from its native support and transform facts,
construction must report the missing fact. It must not restore old WL
materialization, add a persistent route object, or infer a replacement policy.
