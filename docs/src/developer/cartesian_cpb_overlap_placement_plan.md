# CPB Overlap Placement Plan

Date: 2026-06-09

This note defines the next retained/global placement boundary for local CPB
overlap data. It is a design plan only. It does not implement placement, make
route-global overlap available, or change private overlap behavior.

Current state:

- local overlap provider data exists through
  `parent_overlap_axis_factor_packet -> cpb_interval_pair ->
  cpb_overlap_axis_blocks -> cpb_overlap_dense_block ->
  cpb_local_overlap_block_record -> cpb_local_overlap_block_collection`;
- a real probe-enabled report can build one collection from the
  `(:product, :product)` source pair with dense local block shape `(25, 25)`;
- the private overlap side recognizes that collection as structured local
  overlap source data;
- placement requirement and placement candidate fingerprints track missing and
  placeholder placement facts;
- global overlap remains blocked, and even all placeholder facts still block on
  `:placement_not_implemented`.

## Input Object

Retained/global placement must consume a
`CartesianCPBBlockProviders.CPBLocalOverlapBlockCollection`.

It must not consume:

- loose dense matrices;
- scalar report-field clouds;
- old support-row or fixed-block data as route authority;
- report aliases that duplicate local provider metadata.

Each collection record should reference its local provider output through
`source_block`. The source may be a `CPBOverlapDenseBlock` when local dense
materialization has been requested, or a `CPBOverlapAxisBlockSet` when a later
placement path can transform axis blocks directly. Summaries stay compact and
must not duplicate dense matrix payloads.

## Placement Source Record

The placement layer should treat each local overlap record as a structured
source record. The facts that matter are:

- `block_key` or pair key;
- `term = :overlap`;
- left and right CPB summaries;
- interval-pair summary;
- dense-block summary or axis-block summary;
- source kind, for example `:cpb_overlap_dense_block`;
- local ordering;
- factor-space, factor-convention, normalization, and index-domain labels;
- placement status;
- retained-transform status;
- route/global nonclaim flags.

For the current overlap provider, the local ordering is
`:parent_compatible_x_slowest_z_fastest`, and overlap factors carry the
`factor_space = :parent_axis_bundle_pgdg_intermediate`,
`factor_convention = :axis_bundle_one_body_overlap`, and
`index_domain = :parent_axis_indices` contract. Placement must preserve these
labels in compact summaries so later route code does not infer semantics from
field names or matrix sizes.

## Required Placement Facts

Before any retained/global overlap placement can occur, the adapter must have
explicit facts for:

- left retained transform;
- right retained transform;
- left retained column range;
- right retained column range;
- global or retained dimension;
- placement plan;
- accumulation rule.

These facts must be structured enough to be validated independently. They
should not be represented as unrelated scalar fields copied through route
reports. A later implementation can introduce a compact placement-plan object
or placement-candidate object that owns these statuses.

## Schematic Operation

For a local dense CPB overlap block `O_cpb`, the intended retained block shape
is:

```text
O_retained = T_left' * O_cpb * T_right
```

where `T_left` maps left retained columns into the left local CPB product-space
ordering, and `T_right` maps right retained columns into the right local CPB
product-space ordering.

The retained block is then accumulated into the route/global overlap matrix:

```text
global_overlap[left_column_range, right_column_range] += O_retained
```

This is schematic only. No placement code should be written until the
transform, range, dimension, and accumulation contracts are explicit and
covered by focused tests.

## Transform Contract

Each retained transform must state what it maps and which ordering it expects.
For the current local dense overlap block, the source ordering is local CPB
product-space ordering with x slowest and z fastest:

```text
:parent_compatible_x_slowest_z_fastest
```

The transform contract must define:

- source local support count and shape;
- target retained column count;
- source ordering;
- target retained-column ordering;
- whether it is a direct local-to-retained transform or a retained-unit
  transform with additional metadata;
- any normalization convention relevant to overlap placement.

If a transform uses a different local ordering, it must carry an explicit
permutation or conversion. Placement must not silently assume that matrix rows
or columns match the CPB provider ordering.

## Global Dimension Source

Placement must record a `global_dimension_source`.

Preferred sources are structured retained layout facts such as retained-unit
column ranges, terminal retained units, or a reviewed placement plan that owns
the full retained column inventory.

`retained_dimension` is a compatibility fallback only. It is acceptable for
fingerprints and transition tests when no stronger structured retained layout
object exists, but the result must label that source explicitly, for example:

```text
global_dimension_source = :retained_dimension_compatibility
```

A future route-ready placement path should prefer a source such as:

```text
global_dimension_source = :retained_layout_column_ranges
```

or another label owned by the reviewed placement object.

## Accumulation Rule

Before numerical placement, the implementation must choose and record an
accumulation rule. The rule must state:

- whether local retained blocks are added into global ranges;
- how diagonal pairs are handled;
- how off-diagonal pairs are handled;
- whether transpose or symmetry fill is inferred or explicit;
- whether duplicate records for the same range are allowed and summed;
- how blockers are reported when two records conflict.

The simplest first implementation should require explicit left and right
column ranges and should add only the provided block into those ranges. Any
transpose or symmetry fill should be explicit in the placement plan rather than
inferred from `block_key` alone.

## Blockers And Statuses

The current blocker vocabulary remains valid for the planning boundary:

- `:missing_retained_transform`;
- `:missing_left_column_range`;
- `:missing_right_column_range`;
- `:missing_global_dimension`;
- `:missing_placement_plan`;
- `:missing_accumulation_rule`;
- `:missing_placement_or_retained_transform`;
- `:placement_not_implemented`.

Potential refined blockers for a future implementation include:

- `:retained_transform_ordering_mismatch`;
- `:left_column_range_dimension_mismatch`;
- `:right_column_range_dimension_mismatch`;
- `:local_dense_block_shape_mismatch`;
- `:unsupported_overlap_accumulation_rule`;
- `:conflicting_overlap_placement_records`.

This pass does not rename existing statuses or change code. A future
implementation should add refined blockers only when a concrete validation
branch needs them.

## Nonclaims

This plan does not implement:

- placement code;
- retained transforms;
- retained column-range validation;
- route/global overlap availability;
- route-global overlap stage adoption;
- kinetic, position, x2, or Coulomb placement;
- Hamiltonian assembly;
- IDA/MWG semantics;
- PQS Lowdin or projection;
- exports or artifacts.

Local CPB overlap collection availability is not global overlap availability.
It means only that structured local CPB product-space overlap data exists and
is ready for a reviewed placement layer to consume.

## Next Code Step

The next likely implementation is a metadata-only placement plan skeleton, not
numerical placement. That skeleton should carry the collection reference,
record-level placement facts, transform and range statuses, global dimension
source, accumulation-rule status, compact summaries, and blockers. It should
continue to block global overlap until a reviewed numerical placement engine is
implemented and tested.

The private placement plan skeleton now exists as a metadata-only status
carrier. It records local collection, record-placement, transform, range,
dimension-source, placement-plan, and accumulation-rule statuses, but it is not
a placement engine. It must remain blocked until numerical placement is
reviewed and implemented separately.

The current real probe-enabled report audit finds a structured dimension source
through `report.retained_units`, so `report.retained_dimension` is not the only
dimension source in that fixture. It still does not carry structured retained
transforms, source-pair retained column ranges for `(:product, :product)`, a
reviewed overlap placement plan, or an accumulation rule. The next code step
should design those structured carry objects before numerical placement.
