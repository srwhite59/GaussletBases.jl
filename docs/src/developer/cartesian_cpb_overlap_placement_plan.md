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

## Current Metadata Status

Implemented so far:

- private placement plan skeleton;
- `CPBRetainedTransformCarry`;
- `CPBSourcePairPlacementRange`;
- `CPBOverlapPlacementFacts`;
- `CPBReviewedOverlapPlacementPlan`;
- private skeleton consumption of `CPBOverlapPlacementFacts`.

Still missing:

- real retained transforms;
- source-pair retained column ranges;
- numerical transform application;
- global overlap accumulation.

The private placement plan skeleton is a metadata-only status carrier. It
records local collection, record-placement, transform, range, dimension-source,
placement-plan, and accumulation-rule statuses, but it is not a placement
engine. It must remain blocked until numerical placement is reviewed and
implemented separately.

The current real probe-enabled report audit finds a structured dimension source
through `report.retained_units`, so `report.retained_dimension` is not the only
dimension source in that fixture. It still does not carry structured retained
transforms, source-pair retained column ranges for `(:product, :product)`, a
reviewed overlap placement plan, or an accumulation rule.

`CPBOverlapPlacementFacts` is implemented as a metadata-only coherence bundle,
and the private placement skeleton can consume it directly. The current real
probe report negative fingerprint builds facts from the real local overlap
collection without placeholders. It reports:

- `available_requirements = (:local_cpb_overlap_collection,)`;
- `missing_requirements = (:missing_retained_transform,
  :missing_left_column_range, :missing_right_column_range,
  :missing_global_dimension, :missing_placement_plan,
  :missing_accumulation_rule)`.

Those facts preserve left and right CPB summaries into the private placement
skeleton. Missing placement ranges also report `:missing_global_dimension`
because global dimension is an independent required placement fact. Non-matrix
retained-transform references block with
`:unsupported_retained_transform_reference`.

There is still no transform application, placement, global overlap
accumulation, route adoption, kinetic, position, x2, Coulomb, Hamiltonian,
IDA/MWG, PQS Lowdin/projection, export, or artifact work.

## Reviewed Placement Plan Object

The reviewed overlap placement plan object is metadata-only. It owns:

- placement plan kind;
- accumulation rule;
- symmetry or transpose policy;
- duplicate record policy;
- accepted block keys and record inventory;
- required global dimension source;
- status and blocker;
- route/global nonclaim flags.

It should not apply transforms or assemble any matrix. Its purpose is to
replace placeholder `placement_plan` and `accumulation_rule` values with a
compact reviewed contract before numerical placement code exists.

The reviewed placement plan object now exists as metadata-only contract data.
`CPBOverlapPlacementFacts` and the private placement skeleton can use it instead
of placeholder `placement_plan` and `accumulation_rule` values. A real-report
fingerprint can supply this reviewed plan, making `:placement_plan` and
`:accumulation_rule` available while real retained transforms, source-pair
column ranges, and placement-derived global dimension remain missing. This is
still not numerical placement and does not make route-global overlap available.

`CPBOverlapPlacementFacts` now checks the reviewed plan's accepted record
inventory against the local overlap collection block keys. The check records
compact tuples for provided, accepted, rejected, and duplicate block keys. It
blocks metadata coherence when an available local collection record is not
accepted by the reviewed plan, or when duplicate block keys violate a
`:reject_duplicate_block_keys` policy. It also validates that available local
records match the reviewed plan's local ordering contract, and that available
placement ranges use the reviewed plan's required global dimension source. If
placement ranges are missing, the dimension-source gate is not checked and the
existing missing range/dimension requirements remain the authority. This
remains metadata-only validation: it does not apply transforms, place blocks,
or assemble route/global overlap.

## Structured Carry Objects For Placement

The next implementation boundary should introduce compact carry objects before
any numerical placement code. These objects should be private or internal at
first and should be consumed by the placement skeleton. They should not make
global overlap available by themselves.

### CPBRetainedTransformCarry

`CPBRetainedTransformCarry` should describe one side of one local CPB overlap
record. It should be keyed by side and by local source identity, not by a flat
report field. A transform can be shared by key across records when the same
source CPB and retained range are reused, but the placement skeleton should
resolve that sharing into per-record left/right carry summaries.
The current implementation is metadata-only: it validates shape, count, and
ordering metadata, but does not apply transforms.

Required fields:

- `side = :left` or `:right`;
- `block_key` or pair key;
- source CPB summary or source CPB key;
- source local CPB product-space shape;
- source local ordering, currently
  `:parent_compatible_x_slowest_z_fastest`;
- target retained column count;
- target retained column range;
- transform object;
- transform convention and provenance;
- status;
- blocker.

The transform carry should answer only whether a transform is available and
contractually compatible with the local source shape and retained target. It
should not apply the transform. Validation should check:

- source shape matches the local CPB block side shape;
- source ordering matches, or an explicit permutation/conversion is carried;
- target retained column count matches `length(target_retained_column_range)`;
- transform matrix dimensions match source support count and target retained
  count.

The current metadata contract accepts only `AbstractMatrix` transform objects.
Opaque transform references block with
`:unsupported_retained_transform_reference` until reference-aware metadata is
reviewed.

Expected blockers include:

- `:missing_retained_transform`;
- `:retained_transform_source_shape_mismatch`;
- `:retained_transform_target_count_mismatch`;
- `:retained_transform_ordering_mismatch`;
- `:unsupported_retained_transform_reference`.

### CPBSourcePairPlacementRange

`CPBSourcePairPlacementRange` should describe the retained column ranges for a
local source-pair record. Ranges should remain separate from transform carry
objects. A transform may know its own target range, but the source-pair range
object is the placement authority that says where the left and right retained
blocks would land in the global overlap matrix.
The current implementation is metadata-only: it validates range, dimension,
and transform target-count metadata, but does not place matrices.

Required fields:

- `block_key` or pair key;
- `left_column_range`;
- `right_column_range`;
- global dimension;
- global dimension source;
- range source and provenance;
- status;
- blocker.

Validation should check:

- left and right ranges are present;
- range lengths match the left and right transform target column counts;
- ranges lie inside the global dimension;
- range source is a reviewed retained layout or placement plan, not old
  support-row or fixed-block authority.

Preferred global dimension source is retained-unit column ranges or another
structured retained layout. `retained_dimension` remains a compatibility
fallback only and must be labeled as such.

Expected blockers include:

- `:missing_left_column_range`;
- `:missing_right_column_range`;
- `:left_column_range_dimension_mismatch`;
- `:right_column_range_dimension_mismatch`;
- `:missing_global_dimension`.

### CPBOverlapPlacementFacts

`CPBOverlapPlacementFacts` should bundle the available carry objects for one
local overlap collection. It is the object that a future placement adapter
should consume instead of loose transform, range, dimension, and rule
arguments.

Required fields:

- local CPB overlap collection reference or compact summary;
- left and right `CPBRetainedTransformCarry` summaries by record;
- `CPBSourcePairPlacementRange` summaries by record;
- global dimension and source;
- placement-plan status;
- accumulation-rule status;
- available requirements;
- missing requirements;
- status;
- blocker;
- route/global nonclaim flags.

Transforms are logically per side of each placement record. They may be stored
or cached by CPB key, retained-unit key, or another reviewed source key, but
the facts bundle should expose record-level left/right carry summaries so
shape and range validation is local and auditable.

The facts bundle should remain blocked unless all required carries are
available and a reviewed placement plan plus accumulation rule are present. It
should still not materialize a global matrix; it should only decide whether the
inputs are coherent enough for a later placement engine.
`CPBOverlapPlacementFacts` now exists as a metadata-only coherence bundle and
remains blocked until a placement engine is reviewed.
The private placement plan skeleton can consume this facts bundle as its
structured status source, but it remains a blocked fingerprint until a
placement engine exists.

Nonclaims for these carry objects:

- no placement code;
- no transform application;
- no matrix assembly;
- no global overlap availability;
- no kinetic, position, x2, or Coulomb work;
- no Hamiltonian assembly;
- no IDA/MWG semantics;
- no PQS Lowdin or projection;
- no exports or artifacts.
