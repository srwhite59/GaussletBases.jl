# White-Lindsey Terminal Basis Realization

Status: the four White-Lindsey source families are implemented. Matched
aspect-aware shared-shell counts are approved pending under
`HP-PQS-ASPECTSHELL-FN-01`; current source still uses one scalar contraction
count for those shells. `HP-WLTERM-TEST-01` retains narrow inherited committed
test maintenance, and the aspect test ID separately owns the approved matched
H2 regression. Machine authority owns exact lifecycle and file permission.
This page owns the numerical boundary for realizing native White-Lindsey
retained units as the terminal basis used by the shared producer.

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

For an eligible shared z-axis complete shell, the matched aspect policy in
`pqs_complete_shell_aspect_source_modes.md` supplies one common outer shape:

```text
D = (ns, ns, L)
```

PQS consumes the outer boundary quotient. White-Lindsey uses the
axis-specific inner one-dimensional counts:

```text
K = D .- 2 = (ns - 2, ns - 2, L - 2)
```

Each facet or edge chooses the component of `K` for its free axis. Summed over
all native facets, edges, and corners, the shell count must be
`prod(D)-prod(K)`, equal to the matched PQS shell count. This requirement is
dimension parity, not coefficient or subspace identity.

For the canonical one-layer WL boundary examples:

```text
ns = 4, q = 2:  4^3 - 2^3 = 56 boundary columns
ns = 5, q = 3:  5^3 - 3^3 = 98 boundary columns
```

These examples are cubic outer shapes. The old `26`-column result for
`ns = 4` is not valid retained-count policy. Full-support identity retention
for WL boundary strata is likewise invalid. Nonshared and one-center shells
outside the matched aspect policy continue to use the existing scalar
route-local count.

White-Lindsey z-axis diatomics reject normalized `ns < 4` before route
construction. Eligible standard matched shared-shell constructions require
equal PQS/WL direct-core, complete-shell, slab, and bare final dimensions.
Other supported constructions may still differ when their region eligibility
or retention policy differs. Neither case permits altered shell ownership or
fake compactness.

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
public `ns` normalization, direct-core centering, PQS source dimensions, slab
retention, mapped-COMX defaults, artifacts, public APIs, raw blocks,
RG/MWG/IDA semantics, solvers, ECP, or Cr2 workflow. White-Lindsey consumes
the matched longitudinal `L`; it does not choose `L` or an angular-resolution
scale.

If a WL unit cannot be realized from its native support and transform facts,
construction must report the missing fact. It must not restore old WL
materialization, add a persistent route object, or infer a replacement policy.
