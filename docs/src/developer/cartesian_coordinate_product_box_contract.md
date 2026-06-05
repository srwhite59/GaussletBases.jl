# Cartesian Coordinate Product Box Contract

This developer note defines the Coordinate Product Box (CPB) vocabulary for Cartesian nesting work. It is intended to keep shellification geometry separate from recipe-specific lowering.

A CPB is an axis-aligned product of coordinate intervals. Singleton intervals are allowed, so a CPB can describe a filled rectangular box, a facet or face, an edge, a corner, or a slab. A shell region is not itself a CPB; it is an owned support region described by an outer box, an inner exclusion box, or a union of CPB pieces.

## Shared vocabulary

- `ShellificationRegion`: a disjoint owned support region from atom-growth or another shellification policy.
- `CoordinateProductBox` or `CPB`: a product of coordinate intervals.
- `LoweringSource`: a CPB plus a recipe-specific interpretation.
- `IntermediateRetainedSpace`: the contracted source-domain space where source operators are first built.
- `ShellRealization`: an optional map from intermediate space to final owned support.
- `FinalRetainedUnit`: the column-owning unit used by pair planning and Hamiltonian assembly.

## Low-order White--Lindsey

The low-order White--Lindsey route lowers a shellification region by decomposing its shell support into small disjoint CPB boundary strata: facets, edges, and corners. Facets use one fixed axis and two doside axes. Edges use one doside axis and two fixed axes. Corners are direct parent sites.

Thus low-order White--Lindsey is a boundary-stratum CPB lowering.

## PQS

Projected q-shell lowering uses the same shellification-region idea but starts from a filled source CPB. The filled source CPB is not the shell; it is the coordinate-product source domain for an intermediate contracted mode space.

The intended PQS sequence is:

```text
filled source CPB
-> one-dimensional COMX/source transforms
-> product-box modes
-> boundary COMX-product mode selection
-> intermediate retained source space
-> source-space operator blocks from one-dimensional factors
-> optional shell projection and Lowdin cleanup
-> final shell-realized retained unit
```

The source-space operators are built before shell realization. Shell projection and Lowdin cleanup belong to the final realization stage.

## Implementation guidance

Future contracts should distinguish owned support, source CPB, intermediate retained space, shell realization, and final retained unit. Low-order face, edge, and corner strata are one lowering recipe, not the universal model for pair blocks. PQS should enter the driver as a lowering recipe through a filled source CPB rather than as a parallel route that bypasses the shared geometry and realization vocabulary.
