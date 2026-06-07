# PQS Source-Box Operator Framework

This note is the framework contract for the projected q-shell (PQS) operator
lane. It exists to prevent architecture drift across multi-pass work,
especially after handoffs, reviews, or context compaction.

This is not a public API contract and not a claim that the production route
already consumes this machinery. It is the current design frame for private
PQS operator work and retained-unit unification.

## Purpose

PQS should provide a compact, source-box-first way to build high-order
Cartesian retained spaces and operator blocks.

The practical goal is:

```text
small raw product boxes
-> one-dimensional operator factors
-> retained rules
-> source-box pair operator plans
-> optional shell-realization or support-row adapters
-> retained operator blocks
```

The lane should move PQS toward a usable route for molecular calculations
without turning shell-row support contractions into the algorithm. Shell-row
representations are useful for construction checks, current-route authority
comparisons, and debug oracles. They are not the desired computational
primitive for PQS operator assembly.

## 2026-06-07 Module-Spine Checkpoint

The current private Cartesian module spine is:

```text
CartesianCPB
-> CartesianShellification
-> CartesianTerminalLowering
-> CartesianRetainedUnits
-> CartesianRetainedUnitTransformContracts
-> CartesianUnitPairs
-> CartesianPairOperatorPlans
-> CartesianPairBlockMaterialization preflight/direct-direct, PQS source, and LW adapter checkpoints
```

`CartesianRawProductSources` is the side module for raw product source-box
facts used by the PQS source pilot. It owns source CPBs, source-mode
dimensions, source-mode ordering, and metadata-only axis transform facts. It
does not own retained rules, shell projection, Lowdin cleanup, final retained
units, or pair blocks.

`CartesianPairOperatorPlans` consumes retained-unit pairs plus
`CartesianRetainedUnitTransformContracts` for transform and realization paths.
It must not infer realization paths directly from retained-unit kinds.

`CartesianPairBlockMaterialization` is currently a preflight layer plus:

- local direct/direct one-body pair-block pilots for overlap, position, `x2`,
  and kinetic;
- PQS/PQS raw source-space safe one-body helpers for overlap, position, `x2`,
  and kinetic;
- metadata-only PQS source-to-shell-realization bridge summaries for single
  source results and batches;
- metadata-only final PQS pair-block readiness summaries over those bridge
  summaries;
- metadata-only White--Lindsey boundary-stratum adapter preflight
  classification.

The PQS helpers materialize source-space blocks only. They consume ready
`:pqs_source_pair_preflight` records and caller-supplied 1D source factors,
following the old source-mode product pattern from
`_pqs_raw_product_box_pair_mode_matrix(...)` and the safe-term descriptors from
`_source_box_separable_term_factor_kinds(...)`. They do not apply shell
projection, Lowdin cleanup, final retained PQS transforms, Coulomb/IDA
materialization, full operator assembly, Ham export, or artifact writing.
The bridge summaries only record the planned handoff from source-space blocks
to later shell realization. They record source block term/status, source-mode
dims/counts and ordering, transform/source contract keys, realization paths,
compact status and blocker counts for batches, and materialization flags. They
do not build shell projection, Lowdin data, final retained PQS pair blocks,
Hamiltonians, exports, artifacts, IDA/MWG data, or Coulomb blocks. The private
PQS source safe-term descriptor helper is local selector metadata for the
supported source safe terms. It is not public API and not a route-adoption
surface.

`pqs_source_pair_final_block_readiness_summary(bridge_summary)` consumes either
a single PQS source shell-realization bridge summary or a bridge batch summary.
It reports whether a future final retained PQS pair block could be attempted.
At the current checkpoint the expected status remains blocked by
`:shell_realization_not_materialized`; blocked bridge summaries propagate their
own blockers. This helper is metadata-only and does not build shell projection,
Lowdin data, final retained PQS pair blocks, Hamiltonians, exports, artifacts,
IDA/MWG data, or Coulomb blocks.

For White--Lindsey boundary-stratum pairs,
`CartesianPairBlockMaterialization` now recognizes pair-operator records whose
source-operator path is `:white_lindsey_boundary_stratum_adapter_path`. The
metadata-only preflight materialization path is
`:white_lindsey_boundary_stratum_adapter_preflight`; current readiness remains
blocked by `:white_lindsey_boundary_stratum_pair_block_adapter_not_materialized`.
This is adapter-boundary metadata only. It does not build LW numerical blocks,
coefficient maps, doside transforms, Hamiltonians, exports, artifacts, IDA/MWG
data, or Coulomb.

`white_lindsey_boundary_stratum_adapter_summary(record)` is the matching
metadata-only reuse summary for those preflight records. It records only symbol
targets for future adapters: facet/face -> `_nested_face_product`, edge ->
`_nested_edge_product`, corner -> `_nested_corner_piece`, and facet/edge side
helper -> `_nested_doside_1d`. For batch or plan-level inputs,
`reuse_metadata_available_count` and `reuse_metadata_blocked_count` count only
old-kernel reuse metadata availability, not numerical adapter readiness or
pair-block materialization. It does not call those old kernels.

The PQS guardrails are unchanged: support-row or shell-row contraction is an
oracle/debug path, not the PQS algorithm; shell projection and Lowdin cleanup
belong to shell realization, not raw product-box operator construction;
retained-column diagnostic weights are not final IDA weights; final PQS IDA
weights are first-class unsquared integrals of final retained functions.

## Source Documents

Future work in this lane should read this page first, then the detailed policy
notes:

- `docs/src/developer/raw_product_source_retained_transform_policy.md`
- `docs/src/developer/projected_q_shell_policy.md`
- `docs/src/developer/performance_review_contracts.md`
- external agent background memo:
  `/Users/srw/Dropbox/chatarchive/handoff/external_analysis/pqs_source_box_manuscript_packet_2026-06-01/agent_oriented_background_2026-06-01.md`

The geometry/lowering vocabulary is defined in
[`cartesian_coordinate_product_box_contract.md`](cartesian_coordinate_product_box_contract.md).
That note introduces **Coordinate Product Box** (CPB) as the common
coordinate-product source/support primitive and separates shellification owned
support from filled source boxes, intermediate retained spaces, and
shell-realization maps. Read it before extending this framework toward shared
LW/PQS driver stages.

If those documents and this one disagree, stop and update the framework or
policy explicitly before continuing implementation. Do not resolve the conflict
silently in code.

The external memo is drift-control background, not a public API contract. It
is especially useful for keeping the high-level PQS picture straight: PQS is a
single product-shell object built from a full raw source box and boundary
product modes, not a stitched panel/endcap route; edge and corner modes are
part of the structured shell; shell realization/support-local contraction are
compatibility/oracle paths; structured PQS rank loss is a construction
diagnostic rather than silent compression.

## Weight Vocabulary

PQS weights are first-class objects when they are final IDA weights, but the
word "weight" also appears in provisional and diagnostic contexts. These names
must stay distinct:

- `raw_source_weight`: a provisional source/product/parent weight or integral
  attached before the final retained basis is fixed. Raw-weighted pair-factor
  conversion may use these only at the source/raw object that owns the positive
  normalization convention.
- `final_ida_weight`: the unsquared integral of the final retained basis
  function. This is a first-class object and is required for IDA. A PQS route
  may use IDA weights only after they have been recomputed or propagated
  through the actual final retained transform.
- `retained_column_diagnostic_weight`: a debug/reference value attached to a
  retained column in a private fixture or adapter. It is not an IDA weight and
  must not be used for IDA division unless a reviewed contract explicitly
  promotes it to `final_ida_weight`.

The warning in this framework is therefore not "PQS can never have weights."
The warning is that current retained-column diagnostics, raw source weights,
and shell-realization metadata are not substitutes for final IDA weights.

## Key Spaces

The lane uses four spaces that must remain distinct.

### Parent Lattice

The parent lattice is the full Cartesian gausslet grid. It owns the physical
axis points, parent indexing, PGDG operator data, and any parent-level
coordinate/order contract.

Parent rows are not the preferred loop space for PQS operator blocks when a
smaller raw product source exists.

### Raw Product Source Box

A raw product source box is a small local product domain:

```text
axis intervals + source_mode_dims = (nx, ny, nz)
-> 1D source transforms
-> source product modes
```

Once the axis intervals and total source-mode dimensions are chosen, the 1D
source transforms should be deterministic under the selected rule. For the
current PQS lane, these are COMX/source transforms over the full source box.

Operator factors belong here first. For kinetic:

```text
K_raw =
    Kx x Sy x Sz
  + Sx x Ky x Sz
  + Sx x Sy x Kz
```

where each factor is a 1D cross-axis factor between the left and right source
boxes.

### Retained Rule

A retained rule defines what source-box subspace or source-box columns are kept
before any compatibility adapter is considered. In the compact mode-selected
PQS reference path, this can be boundary COMX-product mode selection.

For product/doside units, this is the product/doside retained transform or an
identity selector over a source box.

The retained rule is an algorithmic object. It is not derived by asking the
current shell-realized fixture what transform it happens to expose.

When a retained rule has an explicit source-mode transform, it can be used in
the operator formula:

```text
O_retained = T_left' * O_raw_product * T_right
```

### Shell Realization

Shell realization is separate from raw product-box operator construction. For
the current shell-supported PQS fixture, realization is:

```text
selected product-box modes
-> project to shell rows
-> full-rank symmetric Lowdin cleanup
-> isometric shell-supported retained representation
```

If the active current-route object is shell-realized, the realization transform
is a compatibility/validation object unless a separate framework update makes
it an algorithmic source-box object. It must not force the operator algorithm
to loop over shell rows except as an oracle.

If the current fixture lacks a compact source-space transform, that does not by
itself mean one should be forced out of the fixture. It may mean the current
fixture is not the algorithmic object. In that case, define the algorithmic
raw-box plan, retained rule, and pair operator plan first, then add an adapter
that compares or realizes that object against shell rows.

### Support-Row Adapter

A support-row adapter maps between retained/source-box facts and explicit
parent or shell support rows for compatibility, diagnostics, and current-route
authority checks.

The adapter owns support indices, support-local coefficients, and shell-row
oracle behavior. It must not define the source-box-first algorithm.

## Intended Data Flow

The intended operator data flow is:

```text
left raw product source box
right raw product source box
-> 1D cross-axis factors
-> raw product-box pair operator
-> left/right retained-rule transforms
-> optional realization adapters only after explicit framework review
-> retained operator block
```

For any pair of retained units that have raw product source boxes:

```text
O_block = T_left' * O_raw(left_source, right_source) * T_right
```

That formula is applied only after the participating retained rules have been
defined as source-box objects. Do not start from a shell-row fixture and try to
derive an algorithmic `T` merely because the formula wants one.

This is the organizing principle for:

- PQS/PQS blocks;
- PQS/product blocks;
- product/product blocks;
- future local/Gaussian one-body source-box blocks;
- future GTO/source-box transfer helpers where applicable.

Support-local contraction may still be used to verify or adapt the result:

```text
support rows + explicit retained coefficients -> oracle block
```

but this support-row path must be labeled as oracle/debug/reference. It is not
the algorithmic target for PQS.

## Required Objects

Future implementation should keep these objects conceptually separate even if
the initial code uses private named tuples or adapters.

The object contract is primary. Define these objects first, then implement
helpers that consume them. Current shell-realized fixtures are validation and
compatibility evidence; they are not the authority for inventing the
source-box algorithm.

### Raw Product Box Plan

Owns:

- parent axis intervals;
- total source-mode dimensions;
- 1D source transforms and local coefficients;
- source-mode ordering;
- PGDG/no-fallback provenance diagnostics where available;
- no retained rule.

Must not own:

- shell-row support coefficients;
- Lowdin cleanup;
- packet construction;
- retained weight semantics.

### Retained Rule

Owns:

- source-box identity;
- retained dimension;
- rule kind, such as boundary selector, product/doside transform, identity
  selector, or direct support adapter;
- source-mode selector/transform when the rule is algorithmic in source-box
  space;
- diagnostics describing whether the transform is intended algorithmic input
  or oracle-only.

Must not pretend retained-column diagnostic weights are positive quadrature
weights or final IDA weights.

Must not use the current shell-realized fixture as the algorithmic retained
rule merely because it carries support-local coefficients.

### Retained Unit

Owns:

- column range;
- role and route provenance;
- support ownership metadata;
- retained transform metadata;
- capability labels for operator families.

The unit may carry support-local coefficients for debug, current-route
compatibility, or authority comparison. Those coefficients do not by
themselves define the intended PQS operator algorithm.

### Source Box Pair Operator Plan

Owns:

- left/right retained units;
- left/right source boxes, when available;
- supported terms;
- 1D cross-factor construction plan;
- selected algorithm path;
- oracle path, if present;
- cost and provenance diagnostics.

The pair plan is the place to say whether a pair has a source-box algorithm
available, is oracle-only, or is unsupported.

### Shell-Realization Or Support-Row Adapter

Owns:

- shell support rows and states;
- shell projection and Lowdin cleanup metadata, when present;
- support-local coefficients for compatibility;
- oracle comparison hooks;
- current-route authority comparison metadata.

Must not own:

- source-box operator policy;
- packet adoption;
- retained PQS weight semantics;
- the decision that a shell-row fixture is the algorithmic PQS object.

## Object Contract Sketch

This section makes the private object contract concrete enough that future
implementation passes can be mechanical. Names may remain private named
tuples or helper-local records at first, but the field ownership and
diagnostic boundaries should not drift.

### `RawProductBoxPlan`

Required fields:

- `object_kind`, such as `:cartesian_raw_product_box_plan_3d`;
- `source_box` and `axis_intervals`;
- `source_mode_dims` as total source-mode lengths, not interior counts;
- `source_mode_count = prod(source_mode_dims)`;
- `axis_transform_plan` and per-axis local coefficient matrices;
- `source_mode_indices` and source-mode ordering;
- optional per-axis centers or interval-center views when pair plans need
  coordinate diagnostics;
- provenance fields: construction helper, metric/backend source, and whether
  numerical fallback was invoked by this helper;
- diagnostics: deterministic under box/dims, PGDG/no-fallback provenance when
  checked, `retained_rule_attached=false`, `packet_adoption=false`.

Invariants:

- The plan is deterministic for a fixed box, source-mode dimensions, and
  transform rule.
- `source_mode_dims` are total source-mode lengths.
- Axis coefficient matrix row counts match the physical interval lengths.
- Axis coefficient matrix column counts match `source_mode_dims`.
- The plan owns no retained rule and no shell-row support coefficients.

Must not own:

- boundary selector semantics as a retained rule;
- shell projection or Lowdin cleanup;
- support-local coefficients;
- retained weights or IDA division semantics;
- packet, fixed-block, QW/Hamiltonian, or public/default adoption state.

### `RetainedRule`

Required fields:

- `rule_kind`;
- `source_box_id` or direct reference to the source-box plan;
- retained dimension and column ordering;
- transform/selector metadata in source-mode space when algorithmic;
- compatibility status when the rule is not algorithmic;
- supported operator families or explicit unsupported fields;
- retained-weight semantics, always explicit;
- diagnostics separating algorithmic input from oracle/debug metadata.

Current and anticipated variants:

| Rule variant | Algorithmic status | Required transform facts | Notes |
|---|---:|---|---|
| Boundary COMX-product mode selection | Algorithmic for raw-box PQS | Boundary mode indices, boundary column indices, retained count, source-mode ordering | This is the compact mode-selected PQS rule. No shell projection or Lowdin is part of this rule. |
| Product/doside retained transform | Algorithmic for product/doside units | Axis intervals, axis coefficient matrices, axis function indices, retained count | Used by current private product/product and PQS/product source-box helpers. |
| Identity/direct source selector | Algorithmic only when the source box itself is the retained space | Source-mode identity or selection, retained count | Useful for simple product slabs or debug fixtures; not a license to reinterpret support-dense rows as product modes. |
| Support-dense direct rows | Compatibility/fallback, not product-box algorithmic | Support indices/states, support-local coefficients | Existing atom-box/support-dense fallback path. It can be an oracle or fallback, not a source-box rule. |
| Shell projection plus Lowdin realization | Adapter/compatibility until reviewed otherwise | Shell support rows, projection matrix shape, Lowdin cleanup shape/diagnostics, isometry diagnostics | Current shell-realized PQS fixture exposes this as a transform fact. It reports `source_box_operator_application_ready=false` because no exact compact source-space transform is available yet. |

Anti-invariants:

- Do not use Lowdin cleanup alone as a full raw-to-retained transform.
- Do not infer an algorithmic retained rule from support-local coefficients.
- Do not divide by retained-column diagnostic weights or mark them positive
  quadrature weights/final IDA weights.
- Do not use the current shell-realized fixture as the algorithmic PQS rule
  unless a future framework update explicitly defines that route.

### `SourceBoxPairOperatorPlan`

Required fields:

- left/right source-box plans;
- left/right retained rules;
- left/right retained counts and output block shape;
- pair kind and pair policy;
- supported terms, initially safe terms only: `:overlap`,
  `:position_x/y/z`, `:x2_x/y/z`, and `:kinetic`;
- term-to-factor mapping, for example kinetic as `(K,S,S) + (S,K,S) + (S,S,K)`;
- per-axis cross factors or enough data to build them from caller-supplied
  axis operators;
- factor provenance: explicit metric/operator data, PGDG provenance when
  checked, and whether numerical fallback was invoked by this helper;
- storage/cost diagnostics: materialized raw pair matrix or streamed factors,
  retained-block materialization, expected scaling category;
- oracle path metadata, if an oracle comparison exists;
- authority/adoption flags: packet/fixed-block/QW/Hamiltonian/public/default
  adoption must be false for private plans.

Invariants:

- Pair operators are built in raw source-box spaces first.
- Retained rules are applied after 1D factor construction.
- Shell-row support-local contraction may be present only as oracle/debug
  metadata.
- Unsupported terms reject rather than silently falling back.

Pair policy labels:

- `:source_box_algorithm_available`: a real source-box rule exists for both
  sides and the term is supported.
- `:oracle_only`: the pair can be checked through support rows but does not
  have an algorithmic source-box retained rule.
- `:support_local_fallback`: existing support-dense fallback path, not PQS
  source-box algorithmic policy.
- `:unsupported`: required source, retained rule, or factor data is missing.

### Shell-Realization / Support-Row Adapter

Required fields:

- representation stage, such as `:shell_realized_pqs_fixture`;
- shell support indices/states and support count;
- source-box facts, if known, but not as algorithmic authority;
- boundary selection metadata when present;
- shell projection stage label and matrix shape;
- Lowdin cleanup stage label, method, shape, rank/cutoff/eigenvalue
  diagnostics;
- support-local coefficient shape and optional coefficient equality checks;
- isometry diagnostics when metric data is supplied;
- compact transform availability and missing reason;
- oracle/compatibility flags: `support_local_oracle_used`,
  `shell_row_oracle_only`, `source_box_operator_application_ready`.

Roles:

- compatibility with current shell-supported fixtures;
- debug/authority comparison against fixed-block or packet fields;
- isometry and coefficient sanity checks;
- future adapter boundary if an explicit source-space realization transform is
  defined.

Anti-patterns:

- Calling the adapter the active PQS algorithm;
- optimizing shell-row support contraction as the source-box strategy;
- deriving a compact transform from the fixture just because the operator
  formula wants one;
- using adapter-retained weights for IDA or positive quadrature semantics.

## Existing Helper Migration Map

| Existing helper/fact | Object role | Status |
|---|---|---|
| `_cartesian_raw_product_box_plan(...)` | `RawProductBoxPlan` | Shared private source-box plan; owns intervals, source dims, 1D transforms, ordering, and provenance. |
| `_pqs_raw_product_box_plan(...)` | PQS raw-box plan plus boundary selector metadata | Algorithmic for mode-selected raw-box PQS references; no shell projection or Lowdin in raw-box operators. |
| `_product_doside_retained_unit_plan(...)` | Product/doside `RetainedRule` | Algorithmic product retained rule used by private product/product and PQS/product source-box helpers. |
| `_product_doside_source_box_pair_plan(...)` | `SourceBoxPairOperatorPlan` for product/product | Private/shadow source-box pair plan for safe terms; current product-staged helpers remain authoritative where applicable. |
| `_pqs_product_source_box_pair_plan(...)` | `SourceBoxPairOperatorPlan` for PQS/product | Private source-box mixed plan using PQS boundary selection and product/doside retained transform. |
| `_pqs_pqs_source_box_pair_plan(...)` | `SourceBoxPairOperatorPlan` for compatible raw-box PQS/PQS fixtures | Private reference/shadow plan for same-box and compatible cross-box raw-box fixtures, not current shell-realized route adoption. |
| `_pqs_shell_realization_plan(...)` | Shell-realization adapter | Projects selected modes to shell rows and applies Lowdin for isometry checks; not raw-box operator construction. |
| `_pqs_product_box_realization_plan(...)` | Setup bundle containing raw-box plan plus shell-realization adapter | Useful setup grouping only if the stages remain separate. |
| `_pqs_current_route_shell_realization_transform_fact(...)` | Shell-realization/support-row adapter fact | Current-route metadata precursor; reports `compact_source_space_transform.available=false` and `source_box_operator_application_ready=false`. |
| `_pqs_current_route_retained_unit_inventory(...)` | Retained-unit inventory/compatibility map | Diagnostic current-route inventory; not an algorithmic retained-rule authority. |
| `_pqs_current_route_retained_pair_inventory(...)` | Pair inventory/diagnostic policy map | Labels shell-realized PQS pairs as oracle-only and product/product as source-box-available. |
| `_pqs_current_route_safe_term_matrices(...)` | Support-local oracle matrix builder | Debug/authority comparison path; must not be optimized as the PQS operator algorithm. |
| `_pqs_current_route_safe_term_authority_comparison(...)` | Current authority comparison | Confirms agreement with existing fixed-block/packet fields where available; no adoption implied. |
| `_pqs_atom_box_support_dense_units(...)` | Support-dense adapter/fallback | Direct support-row retained units; not product/doside and not source-box PQS. |
| `_pqs_contact_cap_product_doside_unit(...)` | Product/doside retained rule adapter for contact cap | Source-box-capable product/doside bridge, still private route vocabulary. |

## Next Implementation Gate

Before implementing any current-route PQS/PQS operator block, the following
must be true in object-contract terms:

- Each side has a `RawProductBoxPlan` with validated source-mode ordering and
  axis factor provenance.
- Each side has an algorithmic `RetainedRule` in source-box space. For current
  shell-realized fixtures, the existing shell-realization fact is not enough
  because it reports no exact compact source-space transform.
- The `SourceBoxPairOperatorPlan` states the supported safe terms, pair policy,
  factor construction, output shape, cost model, and oracle path.
- Any shell-realization/support-row adapter is explicitly separate from the
  retained rule used in `T_left' * O_raw_product * T_right`.
- The validation plan compares source-box output to support-row contraction
  only as an oracle and does not route packet/fixed-block construction through
  the new helper.
- The pass has a clear answer to whether it targets mode-selected raw-box PQS
  fixtures, a new explicit shell-realization retained rule, or current-route
  compatibility only.

If any of those fields are ambiguous, stop for design review instead of
filling the gap with shell-row support-local contraction.

## Current PQS/PQS Source-Block Checkpoint

As of 2026-06-07, the narrow PQS/PQS source-space pilot has this
private-module surface inside `CartesianPairBlockMaterialization`:

```text
pqs_source_pair_overlap_block(record; overlap_1d)
pqs_source_pair_position_block(record; axis, overlap_1d, position_1d)
pqs_source_pair_x2_block(record; axis, overlap_1d, x2_1d)
pqs_source_pair_kinetic_block(record; overlap_1d, kinetic_1d)
pqs_source_pair_one_body_block(record, term; ...)
pqs_source_pair_shell_realization_bridge_summary(result_or_batch)
pqs_source_pair_final_block_readiness_summary(bridge_summary)
```

with matching plan-level plural helpers. The accepted one-body terms are:

```text
:overlap
:position_x, :position_y, :position_z
:x2_x, :x2_y, :x2_z
:kinetic
```

The returned blocks are source-space blocks in raw product source-mode
ordering. Their status flags intentionally distinguish source materialization
from final operator materialization:

```text
source_operator_blocks_materialized = true
final_pair_blocks_materialized = false
shell_realization_materialized = false
operator_blocks_materialized = false
hamiltonian_data_materialized = false
artifacts_materialized = false
```

Do not describe these as final PQS pair blocks. They are source-side inputs to
later shell realization and final retained-block assembly.

The bridge summaries are metadata-only checkpoints over those source-side
inputs. A single-result bridge records source block term/status, source-mode
dims/counts and ordering, transform/source contract keys, realization paths,
and materialization flags. Batch bridge summaries add compact term counts,
status counts, blocker counts, skipped-record counts, and source-mode-ordering
uniformity. The private safe-term descriptor helper is local implementation
metadata for selecting supported safe terms and preserving source-term,
required-factor, batch-path, and blocker facts in one place; it is not public
API and must not be treated as a retained rule or realization object.
Bridge summaries do not build shell projection, Lowdin data, final retained
PQS pair blocks, Hamiltonians, exports, artifacts, IDA/MWG data, or Coulomb
blocks.

The final-block readiness summary is the next metadata-only checkpoint. It
consumes a single bridge summary or a bridge batch summary and reports whether
a future final retained PQS pair block could be attempted. Current source
bridges correctly block with `:shell_realization_not_materialized`, and blocked
bridge summaries propagate their blockers. The readiness summary builds no
shell projection, Lowdin data, final retained PQS pair block, Hamiltonian,
export, artifact, IDA/MWG data, or Coulomb block.

## First PQS/PQS Implementation Target

The first PQS/PQS implementation target was the private mode-selected raw-box
PQS/PQS safe-term **source-space** block. It intentionally does not target the
current-route shell-realized fixture. The reason is object-contract clarity:
mode-selected raw-box PQS has an algorithmic source-space retained rule on both
sides, while the shell-realized current-route fixture is still a
shell-realization/support-row adapter with
`source_box_operator_application_ready=false`.

Target object setup:

- Left and right `RawProductBoxPlan` objects must provide axis intervals,
  total source-mode dimensions, source-mode ordering, per-axis source
  transforms, and factor provenance.
- Left and right `RetainedRule` objects must be boundary COMX-product mode
  selection rules over their respective source boxes.
- The retained rule must provide boundary mode indices, boundary column
  indices, retained count, and source-mode ordering.
- No shell projection, Lowdin cleanup, support-local coefficient matrix, or
  retained-column diagnostic weight may participate in the algorithmic retained
  rule.
- Any support-row data used for comparison must be attached through a separate
  adapter/oracle object.

In-scope fixtures:

- Same-box cubic raw-box PQS/PQS fixture.
- Same-box rectangular `q x q x L` fixture.
- Compatible shifted/cross-box fixture only if existing raw-box metadata
  already proves interval compatibility and source-mode ordering.
- Small private fixtures only; no current-route Be2 shell-realized adoption in
  this first implementation.

In-scope safe terms:

- `:overlap`;
- `:position_x`, `:position_y`, `:position_z`;
- `:x2_x`, `:x2_y`, `:x2_z`;
- `:kinetic`.

These safe terms now exist as source-space helpers. The next PQS/PQS target is
not another source helper; it is the explicit bridge from source-space blocks
through shell realization to final retained PQS pair blocks.

Raw product-box pair construction:

- Build per-axis cross factors from caller-supplied axis operator data.
- Use overlap factors for `:overlap`.
- Use one position or `x2` factor and overlap factors on the other axes for
  position and `x2` terms.
- Use `(K,S,S) + (S,K,S) + (S,S,K)` for kinetic.
- Prefer streaming separable factor products directly to the retained block
  for the implementation path. A fully materialized raw product-box matrix is
  acceptable only as an explicit small-fixture oracle.

Retained block formula:

```text
O_boundary = B_left' * O_raw_product(left_box, right_box) * B_right
```

where `B_left` and `B_right` are boundary mode selection rules. This is the
mode-selected raw-box specialization of:

```text
O_retained = T_left' * O_raw_product * T_right
```

Validation ladder:

- Same-box cubic fixture: compare overlap to identity and all safe terms to an
  explicit product-box boundary-column selection reference.
- Same-box rectangular fixture: repeat the safe-term comparison with
  `q x q x L` where `L != q`.
- Compatible shifted/cross-box fixture: compare inside the private helper to
  explicit raw product-box pair boundary-column selection only if the existing
  raw-box metadata already supports the cross intervals.
- Unsupported terms must reject.
- Support-row contraction may be used only as an optional debug/oracle
  comparison through a separate adapter, never as the reference that defines
  the source-box algorithm.

Explicit exclusions:

- No shell-realized fixture compact-transform extraction.
- No shell projection or Lowdin in the raw-box operator path.
- No support-local shell-row contraction as algorithm.
- No support-local fallback optimization.
- No retained-column diagnostic weights, positive quadrature semantics, or IDA
  division.
- No packet/fixed-block adoption.
- No QW/Hamiltonian, local/ECP/Gaussian/MWG/interaction, public/default, CR2,
  or science-status change.

Future implementation stop conditions:

- Stop if the helper needs a current-route shell-realized compact transform.
- Stop if the helper needs shell rows or support-local coefficients to define
  the retained rule.
- Stop if a private source-box helper needs final IDA weights for an operator
  term before the route has a reviewed final-weight contract; do not substitute
  retained-column diagnostics.
- Stop if the fixture would require packet/fixed-block/QW/Hamiltonian/public
  behavior changes.
- Stop if cross-box compatibility is ambiguous; keep the first pass to
  same-box fixtures.

## Invariants

These invariants are lane-wide.

- Operators for product-structured retained units should be built in raw
  source-box spaces first.
- Shell-row projection and Lowdin cleanup belong to realization, not raw-box
  operator construction.
- The support-local shell-row path is an oracle/debug/reference path unless a
  framework update explicitly says otherwise.
- Retained-column diagnostic weights are debug/reference metadata only. They
  must not become IDA division weights or positive quadrature-weight carriers.
  Final PQS IDA weights are different: they are first-class unsquared integrals
  of the final retained basis functions and must be recomputed or propagated
  through the final transform before any IDA path can use them.
- PGDG analytic provenance and "no numerical fallback invoked by this helper"
  are different claims and should be reported separately.
- Current fixed-block or packet authority remains authoritative until a private
  source-box path has exact checks on the relevant fixtures.
- Private source-box helpers do not imply public/default route adoption.
- A full retained matrix diagnostic is allowed for small fixtures, but it is
  not the intended scaling pattern.
- Object contracts come first. If a fixture does not expose the object needed
  by the source-box algorithm, define the object contract rather than forcing
  the fixture to act as that object.

## Must-Not-Become Patterns

The lane should stop if it starts drifting into any of these patterns.

- "Optimize support-local fallback" as the PQS operator strategy.
- Labeling shell-row contractions as the active PQS algorithm.
- Treating Lowdin cleanup alone as the full source-to-retained transform.
- Treating a current shell-realized fixture as the source-box algorithmic
  object.
- Forcing a compact source-space transform out of a compatibility fixture when
  the object contract has not been defined.
- Mixing shell-row realization into raw product-box operator construction.
- Claiming PQS/product or PQS/PQS production support from oracle-only evidence.
- Treating retained-column diagnostic weights as IDA-safe final weights.
- Adding local/ECP/Gaussian/MWG/interaction work before the safe source-box
  operator path is structurally clean.
- Installing sidecars, packets, QW/Hamiltonian hooks, or public/default routes
  from a diagnostic helper.
- Keeping stale docs or helper names that encode the wrong contract.

## Diagnostic Versus Algorithmic Paths

Use these labels consistently.

Algorithmic source-box path:

- builds 1D source-box factors;
- forms or streams raw product-box pair operators;
- applies retained transforms;
- avoids shell-row loops except for verification.

Shell-realization/support-row adapter path:

- may record source-box facts;
- projects selected modes to shell rows;
- applies Lowdin cleanup;
- returns or checks shell-supported retained columns;
- can validate a source-box retained rule;
- is not automatically the retained rule used in
  `T_left' * O_raw * T_right`.

Support-local oracle path:

- uses explicit support-row coefficients;
- contracts rows directly;
- validates current-route semantics;
- may be slow;
- must not be described as the route algorithm.

Current-route authority comparison:

- compares private helper output to fixed block or sequence packet fields where
  those fields exist;
- proves agreement for those fields and fixtures only;
- does not imply adoption.

The recent `_pqs_current_route_*` helpers should be read through this lens.
In particular, `_pqs_current_route_retained_pair_policy(...)`,
`_pqs_current_route_safe_term_matrices(...)`, and
`_pqs_current_route_safe_term_authority_comparison(...)` are validation and
oracle scaffolding for the current shell-realized route. They are not the
algorithmic policy for PQS shell pairs, and their support-local shell-row
contracts must not be copied into source-box implementation as the target
algorithm.

## Current Evidence

The current evidence supports the framework direction, but it is still private
and partial.

What is established:

- Raw product-box PQS self references exist for safe terms.
- PQS/product source-box mixed blocks exist for overlap, position, `x2`, and
  kinetic, including nontrivial product retained transforms.
- Product/product source-box shadow blocks exist, showing the vocabulary is not
  PQS-only.
- Cross-PQS source-box references exist for compatible raw-box fixtures, and
  same-box plus compatible cross-box blocks are now validated inside the
  private helper against explicit raw product-box boundary-column selection
  references.
- The private all-electron local-Gaussian one-body lane now covers the
  positive Coulomb `gaussian_sum` component used for nuclear attraction.
  Product/product, PQS/product, and PQS/PQS pair families have explicit
  term-table helpers and centered analytic wrappers. The centered wrappers
  generate per-axis factors with the existing `CoulombGaussianExpansion` and
  `gaussian_factor_matrices(...)` machinery, require the analytic primitive
  backend, and then feed the explicit source-box helpers. Nuclear charge and
  attraction sign application remain outside these helpers.
- Commit `770b7be` adds the first private route-shaped safe-term consumer.
  It composes PQS/PQS, PQS/product, and product/product source-box blocks for
  overlap, `position_x/y/z`, `x2_x/y/z`, and kinetic. Every route pair is
  labeled `:source_box_algorithm_available`. Product/product blocks go through
  `_product_doside_source_box_reference_block(...)`, which still compares to
  the existing product-staged retained helpers as authority.
- Commits `95d7b11` and `804bdd9` add the first private raw-box route
  producer checkpoint. Explicit fixture facts now produce the same route
  descriptor through `RawProductBoxPlan -> RetainedRule -> route descriptor`.
  The producer uses left/right mode-selected raw-box PQS retained rules and an
  identity product/doside slab retained rule, then feeds the produced
  descriptor into `_pqs_pqs_product_route_shaped_safe_term_consumer(...)`.
  Sampled validation covers the shifted cubic `q5/L5` fixture and a
  rectangular `q5/L7` fixture with `L != q`; consumer output matches the
  source-box shadow or hand-built route path to roundoff. Timing and
  allocation summaries are captured as diagnostic evidence only, not as
  performance thresholds.
- Commit `047af1d` adds the first private geometry/recipe facts producer for
  that route lane. The helper turns small homonuclear-style fixture facts into
  explicit source-box inputs for the raw-box route producer: `parent_dims`,
  `bond_axis`, `q`/`L` or explicit `source_mode_dims`, `left_start`,
  `right_shift`, and a product-slab fixed index/rule become left/right PQS
  source boxes, a product/doside slab source box, source-mode dimensions,
  metadata, provenance, and diagnostics. This is private fixture
  infrastructure, not a general diatomic route geometry policy or public
  builder. The shifted `q5/L5` and rectangular `q5/L7` samples match the
  explicit-fixture producer and safe-term consumer path to roundoff.
- Commit `17dd86d` validates that this geometry facts producer is mechanically
  axis-general over fixture bond-axis labels. A non-`:z` `:x` fixture emits
  source boxes, source-mode dimensions, product slab fixed-axis metadata,
  retained dimension, pair count, and safe-term consumer output matching the
  explicit route-producer path to roundoff. The same checkpoint adds focused
  guards for invalid bond axes, missing `q` when `source_mode_dims` is absent,
  malformed source-mode dimensions, and source-mode axis lengths below two.
  This is still explicit fixture/recipe infrastructure, not an atom-centered
  or CR2 geometry builder.
- A private q4 current-route retained-unit inventory covers all retained
  columns and has a route-wide safe-term authority comparison.
- A Be2-like strict PQS q5 inventory-shape check pins the 8-unit, 36-pair,
  multi-shared-shell route shape.
- A blockwise timing probe shows that shell-realized PQS/PQS support-local
  oracle contraction is the dominant cost center.
- The object-contract sketch above defines the intended private roles for
  raw product boxes, retained rules, source-box pair operator plans, and
  shell-realization/support-row adapters.
- The first PQS/PQS implementation target is narrowed to private
  mode-selected raw-box PQS/PQS safe-term fixtures.

What is not established:

- Production packet or fixed-block adoption of the source-box algorithm.
- Current-route shell-realized PQS/PQS source-box adoption.
- Nuclear-attraction assembly adoption: the private source-box helpers produce
  positive local-Gaussian `gaussian_sum` blocks only, with no nuclear
  charge/sign application and no route integration.
- ECP or other local-potential production paths.
- MWG/IDA interaction support through retained PQS source-box transforms.
- CR2 validation or public/default route readiness.
- That dense retained or raw product-box pair validation matrices scale to
  larger routes.

## Performance Expectations

The desired cost model is controlled by source-box dimensions and retained
dimensions, not by shell support-row counts.

For PQS/PQS safe terms, the expected algorithmic shape is:

```text
1D source-box factors
-> small source-box operator action
-> small retained transforms
-> retained block
```

The observed Be2-like timing shows why this matters. A sampled shell-realized
PQS/PQS support-local oracle block of shape `(98, 98)` took tens of seconds
for overlap/kinetic even though allocation was small. That is evidence that
the oracle path is CPU-expensive, not evidence that support-row contraction
should become the algorithm.

Future performance reports should separate:

- inventory construction time;
- source-box factor time;
- retained-transform application time;
- oracle comparison time;
- full retained matrix materialization time, when used;
- allocation and dense matrix materialization claims.

## Open Questions

These questions are intentionally unresolved.

- How much route-scale validation should materialize dense raw source-box pair
  oracle matrices versus relying on sampled checks? The private helper streams
  separable factors for the algorithmic block and materializes dense raw
  source-box pair matrices only as small-fixture validation oracles.
- What is the smallest Be2-like sampled validation that proves same-shell and
  cross-shell PQS/PQS blocks without running a full route-wide oracle?
- How should support-dense atom-box fallback enter the broader retained-unit
  vocabulary without pretending to be product/doside?
- What is the next semantic lane after positive local-Gaussian one-body
  factors: nuclear charge/sign assembly, current-route authority comparison,
  or electron-electron design? Electron-electron should be scoped separately.
- What evidence is required before any packet/fixed-block adoption discussion?

## Framework Update Rule

Update this document before implementation when any of these changes are
proposed:

- changing the meaning of PQS retained transforms;
- treating shell-row support-local contraction as more than an oracle;
- adding a new retained-unit kind to the PQS route vocabulary;
- adding a new operator family such as local/Gaussian/MWG/IDA;
- moving from diagnostic helpers to construction, packet, fixed-block, QW, or
  public/default adoption;
- changing weight semantics;
- changing performance expectations or validation gates.

Every PQS baton loop should cite this document in `RUN.md` and in each blurb.
Each architecture-sensitive PQS baton loop should also create or refresh a
loop-local `INVARIANTS.md` that quotes the governing invariants from this
framework. Before each pass and before each commit, reread this framework and
the loop-local invariants file. If implementation pressure suggests the
framework is wrong or incomplete, stop and make the framework update explicit.

## Next Intended Correction

The private route-shaped raw-box safe-term consumer checkpoint is commit
`770b7be`. It is a route-shaped consumer, not route adoption: it takes an
already-built private descriptor with left PQS raw plan, right PQS raw plan,
and product/doside unit, then delegates numerical blocks to the source-box
helpers. PQS/PQS uses `_pqs_pqs_source_box_reference_blocks(...)` with
helper-internal explicit raw product-box boundary-selection validation.
PQS/product uses `_pqs_product_source_box_reference_blocks(...)`.
Product/product uses `_product_doside_source_box_reference_block(...)`, so it
is labeled as source-box vocabulary while still checking against the existing
product-staged retained helpers.

The algorithmic path remains source-box first. Dense raw source-box pair
matrices are validation-only. The consumer does not use shell projection,
Lowdin cleanup, support-local fallback, support coefficient matrices, retained
PQS weight semantics, or IDA division, and it does not change packet,
fixed-block, QW/Hamiltonian, public/default, local/ECP/Gaussian/MWG/
interaction, IDA/MWG, or CR2 behavior.

The private mode-selected raw-box PQS/PQS safe-term target now has helper
coverage for same-box cubic, same-box rectangular, and compatible shifted
cross-box fixtures. The helper streams 1D factor products into retained
blocks; dense raw product-box pair matrices are materialized only for explicit
small-fixture validation oracles. A current-route shell-realized PQS/PQS block
should still wait until an explicit source-space realization rule is defined
or the pass is scoped as compatibility/oracle-only.

The private raw-box route producer checkpoint is commits `95d7b11` and
`804bdd9`. It is the first producer-side complement to the route-shaped
consumer. It starts from explicit fixture facts, builds left/right
`RawProductBoxPlan` objects, attaches boundary COMX-product mode-selection
`RetainedRule` facts for the two PQS units, creates an identity
product/doside slab retained rule for the middle unit, and emits the existing
route descriptor shape consumed by
`_pqs_pqs_product_route_shaped_safe_term_consumer(...)`.

This producer is still private/shadow-only. The sampled validation matrix
covers the shifted cubic `q5/L5` fixture and a rectangular `q5/L7` fixture
with `L != q`. The produced descriptor output matches the source-box shadow
and hand-built descriptor paths to roundoff. Producer and consumer timing plus
allocation summaries are recorded only as diagnostics; they are not readiness
thresholds. Dense raw source-box matrices remain validation-only and are not
the producer algorithm.

Commit `047af1d` records the first geometry/recipe facts producer checkpoint
for this lane. `_pqs_pqs_product_raw_box_homonuclear_geometry_facts(...)`
accepts `parent_dims`, `bond_axis`, `q`/`L` or explicit `source_mode_dims`,
`left_start`, `right_shift`, and a product slab fixed index/rule. It emits the
explicit left/right PQS source boxes, product/doside slab source box,
source-mode dimensions, metadata, provenance, and diagnostics consumed by
`_pqs_pqs_product_raw_box_route_from_geometry_facts(...)`, which then calls
the existing raw-box route producer. This is a private fixture geometry
producer only: it is not a broad diatomic geometry policy, public builder, or
operator-algebra authority.

The geometry checkpoint validates the shifted `q5/L5` and rectangular
`q5/L7` samples against the explicit-fixture route producer and safe-term
consumer path to roundoff. It adds no shell projection, Lowdin cleanup,
support-local fallback as an algorithm, support coefficient matrices,
retained-column diagnostic weights, IDA division, packet or fixed-block
adoption, QW/Hamiltonian routing, public/default behavior,
local/ECP/Gaussian/MWG/interaction terms, IDA/MWG change, or CR2 science
claim.

Commit `17dd86d` records the focused axis-general validation checkpoint for
the same helper. The helper is now validated for `:z` and `:x` fixture bond
axes: the `:x` sample emits the expected left/right PQS source boxes,
product/doside slab source box with fixed axis `1`, source-mode dimensions,
retained dimension, pair count, pair policy, and safe-term consumer output
matching the explicit route-producer path to roundoff. This axis check is
mechanical; it does not add center inference, atom boxes, shell-realized
current-route PQS, or broad diatomic route geometry policy.

The same checkpoint adds high-level guard coverage for invalid `bond_axis`,
missing `q` when `source_mode_dims` is absent, malformed source-mode
dimensions, and source-mode axis lengths below two. Those guards keep the
private geometry facts lane explicit:
`geometry/recipe facts -> explicit source boxes/source-mode dimensions ->
RawProductBoxPlan -> RetainedRule -> route descriptor -> source-box safe-term
consumer`.

Commit `28c3dbc` records the private input-gate checkpoint for this producer.
The gates are misuse protection for private fixture work, not a public route
contract. Source boxes must be nonempty and inside `parent_dims`; `parent_dims`
must be a positive 3D integer tuple; PQS source-mode dimensions are total
source-mode lengths and must have at least two modes per axis for boundary
selection; the identity product/doside slab requires exactly one fixed axis;
and unsupported safe terms reject before source-box construction.

These producer gates do not change the source-box-first boundary. They add no
shell projection, Lowdin cleanup, support-local PQS oracle, support
coefficient matrix use, retained-column diagnostic weights, IDA division,
packet or fixed-block adoption, QW/Hamiltonian routing, public/default
behavior, local/ECP/Gaussian/MWG/interaction terms, IDA/MWG change, or CR2
claim.

The private source-box local-Gaussian one-body checkpoint is commits
`b12f1b4`, `7608264`, `d948d51`, `c262d40`, `5dc0642`, and `dea7feb`.
It covers the positive Coulomb `gaussian_sum` component used by all-electron
non-ECP nuclear attraction:

```text
sum_t c_t * Gx_t * Gy_t * Gz_t
```

The pair families covered are product/product, PQS/product, and PQS/PQS. Each
family has a caller-supplied explicit term-table helper, and each has a
centered analytic wrapper where the per-axis term tables are generated from
the existing `CoulombGaussianExpansion` and
`gaussian_factor_matrices(...)` machinery. These wrappers require the
analytic primitive backend and then delegate to the explicit source-box
helpers.

This checkpoint is still private/source-box/reference infrastructure. The
helpers listed above build positive `gaussian_sum` blocks only. Commit
`549ae2f` adds the separate physical all-electron non-ECP wrapper layer:

```text
V_nuc,A = -Z_A * gaussian_sum(center_A)
```

The physical wrapper preserves one block per nucleus/center as the primary
result because counterpoise workflows need to keep center contributions
separate. A summed total block is allowed only as a derived convenience. The
positive source-box helper output remains available as component metadata, so
the sign/charge wrapper does not obscure which analytic Gaussian terms were
used.

Focused validation uses small fixtures and ignored probes covering
product/product, PQS/product, same-box PQS/PQS, rectangular PQS/PQS, and
shifted/cross-box PQS/PQS comparisons against explicit source-box references
or explicit-table helper output. The nuclear wrapper checks include sign,
charge scaling, by-center preservation, and rejected malformed center/charge
inputs.

The checkpoint adds no ECP terms, electron-electron terms, MWG/IDA change,
retained-column diagnostic positive-weight semantics, retained-column
diagnostic-weight IDA division, shell-row support-local contraction as an
algorithm, packet or fixed-block adoption, QW/Hamiltonian routing,
public/default behavior, or CR2 science claim. Starting electron-electron
interactions is a separate semantic lane and should wait for explicit review.

### Electron-Electron Source-Box Design Lane

Electron-electron source-box work is a separate semantic lane, not a
continuation of the one-body nuclear sign/charge wrapper.

The existing ordinary/MWG/IDA path has two distinct pair-factor objects:

- raw pair-factor terms, such as `pair_factor_terms_raw`, which live at the
  auxiliary/raw quadrature level;
- density-normalized pair factors, such as `pair_factor_terms`, produced by
  dividing raw pair terms by raw auxiliary weights in the existing IDA/MWG
  convention.

That weight division belongs to the raw/source or support quadrature object
that owns the positive weights. It must not be moved onto retained PQS
columns. Retained PQS columns are transformed source-box columns, not positive
quadrature carriers, and retained-column diagnostic weights must continue to
report semantics as not IDA-safe. This does not demote final PQS IDA weights:
those are separate first-class unsquared integrals of the final retained basis
functions, as defined above.

The first electron-electron object contract should describe:

- `RawPairFactorData`: per-axis pair-factor term tensors, raw weights, density-
  normalized pair factors if already produced by the intended IDA/MWG
  convention, expansion coefficients, backend/provenance, and a flag saying no
  numerical fallback was invoked by this helper;
- left/right `RawProductBoxPlan` and `RetainedRule` objects, as used by the
  safe-term source-box lane;
- a retained-pair density object or matrix-block object whose output is in the
  repo's two-index density-density interaction convention, not a four-index
  Galerkin Coulomb tensor;
- explicit diagnostics for where source weights live and for
  `retained_weight_division_allowed = false`.

The smallest first fixture was product/product, because product/doside units
already have source-box retained transforms and current product-staged helpers
can provide comparison data. PQS/product and PQS/PQS now have corresponding
private density-normalized source-box blocks. Current support-local or
fixed-block interaction paths may be used only as validation/oracle
comparisons until an explicit source-box pair object is promoted into a real
route. No promotion has happened in this lane.

Commit `9bed286` adds that first tiny product/product fixture. It consumes
caller-supplied density-normalized pair factors, requires raw/source weights
as provenance, does not divide by weights again, and returns a two-index
density-density block. Raw-weighted inputs intentionally rejected at that
checkpoint rather than guessing the convention.

The raw-weighted conversion boundary is now explicit for product/product:
per-axis raw pair-factor terms are converted to the density-normalized
convention by dividing each term matrix by the raw source weight outer product
on that axis:

```text
F_density[t][i,j] = F_raw[t][i,j] / (w[i] * w[j])
```

The conversion is axis-wise and term-wise, matching the existing ordinary
IDA/MWG construction of `pair_factor_terms` from `pair_factor_terms_raw`.
After conversion, the density-normalized helper remains the low-risk core.
Diagnostics distinguish `pair_factor_normalization = :raw_weighted`,
`source_weight_division_owner = :source_box_raw_weights`, and
`source_weight_division_applied_by_helper = true`.

Commit `27d15dd` adds the private PQS/product density-normalized source-box
density-density block. It uses mode-selected raw product-box PQS on the PQS
side and the product/doside retained transform on the product side. Inputs are
caller-supplied density-normalized per-axis pair factors; raw/source weights
are checked and carried as provenance only. The helper does not divide by
weights and does not assign retained PQS columns positive quadrature or
IDA-safe weight semantics.

Commit `b1ee2a5` adds the matching private PQS/product raw-weighted
conversion wrapper. It accepts raw-weighted per-axis pair factors plus
explicit positive source weights, divides by the source-weight outer product,
then delegates to the PQS/product density-normalized core. This is raw/source
quadrature-weight normalization only; retained PQS columns still do not carry
positive quadrature weights, and retained-column diagnostics are not
IDA-division weights.

Commit `653d35d` adds the private PQS/PQS density-normalized source-box
density-density block. It uses mode-selected raw product-box PQS on both
sides, applies boundary COMX-product mode selection on both retained rules,
and returns a retained two-index density-density block. Inputs are
caller-supplied density-normalized per-axis pair factors; explicit positive
source weights are required only as provenance and positivity checks. The
helper does not use shell projection, Lowdin cleanup, support coefficient
matrices, support-local/shell-row contraction as the algorithm, retained PQS
diagnostic weights, or retained-column diagnostic-weight/IDA division.

Commit `3b64e91` adds the matching private PQS/PQS raw-weighted conversion
wrapper. It accepts raw-weighted per-axis pair factors plus explicit positive
source weights, divides each raw term matrix by the source-weight outer
product on that axis, and delegates to the PQS/PQS density-normalized core.
This closes the private pair-family symmetry for product/product,
PQS/product, and PQS/PQS without assigning positive quadrature or IDA-division
semantics to retained PQS columns.

Commit `93a9af8` adds the private route-shaped source-box density-density
consumer:

```text
left mode-selected raw-box PQS
right mode-selected raw-box PQS
middle product/doside slab
-> six upper-triangular source-box pair blocks
-> complete retained two-index density-density matrix
```

This consumer is private/reference infrastructure only. It composes the
already-proven pair-family helpers: three PQS/PQS pairs, two PQS/product
pairs, and one product/product pair. Product/PQS and other lower-triangular
cross blocks are filled by transpose only after the synthetic/caller-supplied
pair-factor terms are checked as symmetric. The output remains the repo's
retained two-index density-density convention, not a four-index Galerkin
Coulomb tensor. Inputs are still explicit synthetic or caller-supplied
per-axis pair-factor terms; real repo MWG/IDA pair-factor provenance has not
been adapted.

The route consumer supports both density-normalized input and a raw-weighted
route mode. In raw-weighted mode, normalization is still source-box
normalization only: raw per-axis pair factors are divided by explicit
raw/source quadrature-weight outer products through the existing raw-weighted
wrappers, then the density-normalized cores are used. Retained PQS columns
still do not carry positive quadrature weights, retained-column diagnostic
weights are not used, and retained-column diagnostic-weight/IDA division
remains forbidden. The consumer adds no shell projection, Lowdin cleanup,
support-local/shell-row algorithm, support coefficient matrix,
packet/fixed-block adoption, QW/Hamiltonian routing, MWG/IDA semantic change,
public/default route, ECP behavior, or CR2 science claim. It records
timing/allocation diagnostics for complete retained-matrix assembly only as
private route-consumer diagnostics.

Commit `e868f49` adds the matching private density-density route producer
wrapper, `_pqs_pqs_product_raw_box_density_density_route_producer(...)`. It
starts from explicit source-box fixture facts for the same
left-PQS/right-PQS/product-slab layout, builds the raw-box route descriptor
through the existing private route producer machinery, then delegates to the
route-shaped density-density consumer above. The producer returns both the
descriptor and the consumer result so tests and future review passes can
compare the descriptor-only consumer path with the explicit-fixture producer
path. This remains private/reference infrastructure; it is not production
interaction adoption.

The producer supports the same synthetic or caller-supplied pair-factor modes
as the consumer. Density-normalized factors are accepted as already
normalized. Raw-weighted factors are normalized only by explicit raw/source
quadrature-weight outer products in the existing raw-weighted wrappers before
delegating to the density-normalized helpers. Retained PQS columns are still
not positive quadrature weights and are not used for retained-column
diagnostic-weight or IDA division. Every pair remains on the source-box
algorithmic path; shell projection, Lowdin cleanup, support-local/shell-row
contraction, support coefficient matrices, packet/fixed-block adoption,
QW/Hamiltonian routing, real MWG/IDA pair-factor provenance, public/default
behavior, ECP behavior, and CR2 science claims remain out of scope. Dense raw
product-box matrices, if mentioned in diagnostics, are validation-only.

Commit `3028556` adds a private diagnostic PGDG IDA source-factor provenance
extractor, `_pqs_source_box_ida_factor_provenance(...)`. It records the real
repo IDA gausslet/source-box factor data from PGDG intermediates: density-
normalized `pair_factor_terms`, raw `pair_factor_terms_raw`, explicit
source/raw quadrature `weights`, centers, backend metadata, and shape/term-
count diagnostics. The extractor adapts only the IDA gausslet/source-box
layer. MWG supplement/residual coupling remains separate and is explicitly not
consumed. Retained PQS columns are not introduced as positive quadrature
weights or IDA-division weights.

Commit `5de13b1` adds the private diagnostic route adapter
`_pqs_pqs_product_raw_box_density_density_route_producer_from_ida_provenance(...)`.
It feeds the IDA provenance object into the existing explicit-input route
producer for the same left-PQS/right-PQS/product-slab source-box layout. The
adapter delegates to `_pqs_pqs_product_raw_box_density_density_route_producer(...)`
and does not duplicate route pair assembly. Density-normalized mode passes
the provenance `pair_factor_terms`; raw-weighted mode passes
`pair_factor_terms_raw` and uses only the explicit PGDG source/raw quadrature
weights through the existing raw-weighted wrappers. The output remains a
retained two-index density-density matrix, not a four-index Galerkin Coulomb
tensor. Diagnostics record
`input_pair_factor_data = :ida_gausslet_source_box_provenance`,
`interaction_path = :ida_gausslet_source_box`,
`real_ida_gausslet_source_box_provenance_adapted = true`, and
`mwg_supplement_residual_path = false`. The generic mixed
MWG/IDA-adapted flag stays false to preserve the IDA/MWG split.

Commit `c443af2` adds a private validation-only dense-parent IDA authority
comparison for the same small route. `_pqs_pqs_product_route_parent_coefficient_matrix(...)`
builds a retained-to-parent coefficient map from the source-box route facts:
left/right PQS columns are embedded from mode-selected raw product-box axis
coefficients and the product slab uses its product/doside retained
coefficients. `_pqs_pqs_product_dense_parent_ida_authority_comparison(...)`
then compares the source-box IDA adapter block with

```text
C_route' * V_parent_ida * C_route
```

where `V_parent_ida` is the existing dense parent IDA matrix, such as the
matrix from `_qwrg_diatomic_interaction_matrix(...)`. The focused `q5/q5/q7`
fixture agrees to roundoff, with dense-parent projection max error about
`1.8e-15`. This is an authority comparison only. Dense parent projection is
not the source-box algorithm, and it does not adopt packet/fixed-block/
QW/Hamiltonian/public/default behavior. MWG supplement/residual coupling
remains separate and unadapted, and retained PQS columns remain non-
quadrature retained columns with no retained-column diagnostic-weight/IDA
division.

Commits `33a76ed` and `30fd052` add the private component route smoke helper
`_pqs_pqs_product_source_box_component_route_smoke(...)` for the same small
source-box layout:

```text
left mode-selected raw-box PQS
middle product/doside slab
right mode-selected raw-box PQS
```

The component smoke combines the already-private one-body and two-body
source-box pieces without becoming a Hamiltonian route. Nuclear attraction is
reported as one retained matrix per nucleus/center, preserving by-center
identity for counterpoise workflows, and the total nuclear-attraction matrix
is only an explicit sum of those per-center pieces. The electron-electron
component remains an IDA gausslet/source-box retained two-index
density-density block, not a four-index Galerkin Coulomb tensor.

The component smoke supports both `:density_normalized` and `:raw_weighted`
electron-electron factor modes. Raw-weighted mode is mechanical: it delegates
to the existing raw-weighted IDA route producer path, where explicit
raw/source quadrature weights own the normalization before the
density-normalized cores are used. Dense parent IDA authority is still
validation-only and, at this checkpoint, density-normalized-only. MWG
supplement/residual coupling remains separate and unadapted. The helper adds
no retained-column diagnostic weights, retained-column diagnostic-weight/IDA
division, shell/support-local algorithm, packet/fixed-block/QW/Hamiltonian
adoption, public/default route, ECP behavior, or CR2 science claim.

The ignored private report probe
`tmp/work/pqs_component_route_smoke_report.jl` records this component route
checkpoint in generated artifacts under
`tmp/work/pqs_component_route_smoke_report_outputs/`: `report.txt` and
`summary.tsv`. The first three route-size rows are deliberately small and are
only a private scaling signal, not a production benchmark:

| route variant | retained dim | modes | no-go status | dense authority | density-row timing/allocation signal |
|---|---:|---|---|---|---|
| `q5_L5_parent5x5x7_slab_z4` | `221` | density-normalized and raw-weighted | clear for both rows | density row available; raw row skips with `density_normalized_authority_only` | about `0.87 s` / `394 MB` nuclear and `0.0016 s` / `1.6 MB` electron-electron |
| `q5_L7_parent5x5x9_slab_z5` | `285` | density-normalized and raw-weighted | clear for both rows | density row available; raw row skips with `density_normalized_authority_only` | about `0.79 s` / `416 MB` nuclear and `0.0025 s` / `2.1 MB` electron-electron |
| `q5_L9_parent5x5x11_slab_z6` | `349` | density-normalized and raw-weighted | clear for both rows | density row available; raw row skips with `density_normalized_authority_only` | about `0.80 s` / `421 MB` nuclear and `0.0028 s` / `2.8 MB` electron-electron |

All three variants report six nuclear pairs, six electron-electron pairs, IDA
term count `45`, finite output, and source-box-first no-go diagnostics. In
the same run the largest raw-weighted electron-electron row is the `349`
dimension variant at about `0.0045 s` and `3.8 MB`. The probe does not adopt
packet/fixed-block/QW/Hamiltonian behavior, public/default routing, MWG
supplement/residual coupling, ECP behavior, or CR2 science status.

The private component-route report adapter checkpoint extracts the Be2/PQS
summary shape from script-local report assembly into
`_pqs_pqs_product_component_route_smoke_report_adapter(...)` and
`_write_pqs_pqs_product_component_route_smoke_report(...)`. The ignored
Be2/PQS smoke probe now feeds existing source-box/PQS component summaries and
caller-supplied final-residual MWG component facts into that structured report
object. This remains reporting infrastructure only: it records the
source-box/PQS IDA lane and the ordinary final-residual MWG supplement lane
side by side without merging them into one algorithm. The focused tracked
test `test/nested/pqs_component_route_report_adapter_runtests.jl` uses only
synthetic summary/fact objects, so it avoids `/Users/srw/BasisSets` and other
machine-local basis files. Its strict-mode checks require explicit residual
owner metadata, zero final-residual MWG authority error, no owner inference
from `raw_to_final`, no raw GTO/GTO or fixed/raw-GTO MWG blocks, and no
retained/source-box/final-residual weight or IDA division.

The follow-up private CR2 sidecar schema checkpoint adds
`_pqs_pqs_product_component_route_smoke_cr2_sidecar_schema(...)` as
sidecar/reporting infrastructure. The ignored Be2/PQS smoke probe writes the
generated sidecar summary at
`tmp/work/be2_pqs_component_route_smoke_cr2_sidecar_schema_report.txt`.
The sidecar collects adjacent facts for CR2-facing consumers; it does not
define one unified Hamiltonian route and does not merge the source-box/PQS
IDA lane with the ordinary final-residual MWG supplement lane. The sidecar
makes explicit the source-unit labels and retained ranges from route
descriptor unit keys, fixed/residual/final dimensions and ranges, explicit
residual owners, component block shapes, authority errors, and provenance. It
also records the absences: no shell-label reconstruction from centers, no
nearest-grid or nearest-center label heuristic, no raw GTO/GTO or fixed/raw-
GTO MWG blocks, no owner inference from `raw_to_final`, no retained/source-
box/final-residual weight or IDA division, and no packet/fixed-block/
QW/Hamiltonian/public/default/ECP/SCF/HF/CR2 science adoption. Shell labels
for full Be2 q5/CR2 artifacts are still `:unavailable` in this checkpoint;
full shell/source-unit labels need explicit repo export before shell-start
contraction becomes a contract.

The tracked private driver `bin/cartesian_ham_builder.jl` is an
old-style development spine for the second migration. It follows the simple
`nestpgg.jl` pattern: plain `using` statements, top-level defaults, optional
config include, simple `name=value` overrides, and straight-line printed
diagnostics. Its first checkpoint is metadata/reporting only. It names the
intended stages from system metadata and recipe metadata through parent
description, retained-unit partitioning, pair inventory, and final retained
linear-algebra plan, but it does not build a Hamiltonian or adopt a public
route. The default dry run writes no artifact; optional JLD2/TSV output is a
compact report of the same metadata.

The driver is deliberately private scaffolding. It records the source-box-first
route shape, retained ranges, pair-family helper names, density-normalized or
raw-weighted pair-factor mode, and no-go diagnostics. It must not be treated as
packet/fixed-block/QW/Hamiltonian adoption, MWG/IDA semantic change,
retained-column diagnostic-weight division, repo-side `ray_id`, ECP behavior,
CR2 science status, or a shell/support-local algorithm.

The repo-owned private writer
`_write_pqs_pqs_product_component_route_smoke_cr2_sidecar_schema_report(...)`
now emits that CR2 sidecar schema text report. The ignored Be2/PQS smoke
probe calls this writer instead of carrying a duplicate local text format.
Focused synthetic coverage checks the writer sections and key fields,
including source-unit labels, unavailable shell labels, disabled center/
nearest-grid label heuristics, residual owners, fixed/residual/final
dimensions and ranges, component block shapes, and lane-boundary no-go
records. This writer is still private/reporting-only and adds no unified
Hamiltonian route, shell-label reconstruction, raw GTO/GTO or fixed/raw-GTO
MWG blocks, retained-column diagnostic-weight/IDA division, or packet/fixed-block/
QW/Hamiltonian/public/default adoption.

The private fixed-side retained-unit metadata checkpoint adds
`_pqs_current_route_fixed_side_retained_unit_metadata(...)` as a reporting-only
view of `_pqs_current_route_retained_unit_inventory(...)`. It normalizes the
explicit inventory unit keys, roles, categories, retained ranges, retained
counts, support counts, representation kinds, available source-mode
dimensions, and shell-realized PQS oracle flags into sidecar-adjacent metadata
records. It does not construct a route, install a sidecar, thread data into
the CR2 sidecar schema by itself, or change packet/fixed-block/QW/Hamiltonian
behavior.

The private CR2 sidecar schema can optionally carry this fixed-side metadata
when callers pass `fixed_side_retained_unit_metadata` explicitly. Omitted
metadata preserves existing caller behavior and leaves the text report without
the optional section. When present, the repo-owned writer emits a
`[fixed_side_retained_unit_metadata]` section with the same reporting facts and
the same no-go records: shell labels remain unavailable, center/grid label
reconstruction and nearest-grid heuristics remain disabled, and no
retained-column diagnostic-weight/IDA division or route/adoption behavior is
introduced.

The source-shell/source-mode label lane is provenance metadata, not operator
data. A label such as `(shell_id, ix, iy, iz)` identifies a
construction-native function or retained column in a source shell. It does not
assert support locality, does not define a ray, does not imply a positive
quadrature weight, and is not enough to reconstruct parent-to-final or
source-to-final coefficient maps. Consumers such as CR2 may build ray/cone or
tail-contraction groupings from exported labels, dimensions, and representative
centers, but repo exports facts rather than choosing the grouping policy.
Needed operator matrices should continue to be built in repo as small-small
blocks; labels are not a substitute for exporting a large parent-basis
coefficient map.

The `source_modes.local_axis_*` fields are normalized identity labels in the
row's source shell. They run from `1:nx`, `1:ny`, and `1:nz` for that shell's
contracted dimensions. Parent/support-state coordinates are separate optional
metadata, exported only under explicit parent-lattice field names when the
construction has native facts. They must not be read back as local labels or
as repo-defined rays. The private `source_axis_*` columns are compatibility
native-tuple fields; diagnostics should treat `local_axis_*` and
`parent_lattice_axis_*` as the clarified contract fields.

The ignored full Be2 q5 report
`tmp/work/be2_pqs_fixed_side_retained_unit_metadata/be2_pqs_fixed_side_retained_unit_metadata_report.txt`
exercises that extractor on the current-route inventory. It records fixed
dimension `1483`, retained-unit count `8`, complete coverage `1:1483`,
product/doside count `3`, support-dense/direct-support count `2`,
shell-realized PQS count `3`, and source-box-ready shell-realized PQS count
`0`. The explicit inventory labels are
`:outer_mismatch_z_low_slab`, `:outer_mismatch_z_high_slab`, `:left_atom_box`,
`:right_atom_box`, `:contact_cap_slab`,
`:regular_shared_molecular_shell_1`, `:regular_shared_molecular_shell_2`,
and `:regular_shared_molecular_shell_3`. The three shared PQS retained counts
are `(98, 98, 114)`, support counts are `(1738, 1346, 1002)`, and retained
ranges are `(1174:1271, 1272:1369, 1370:1483)`.

This report makes fixed-side retained-unit labels explicit. It does not make
true shell-start labels available and does not make shell-start contraction a
stable contract. `shell_label_status` remains `:unavailable`; center/grid
label reconstruction and nearest-grid heuristics remain disabled. The
shell-realized PQS fixtures remain metadata/oracle-only and not
source-box-operator-ready. Coarse timings for the ignored run were about
`32.6 s` for construction, `14.1 s` for current-route inventory, and `3.4 s`
for metadata extraction.

The current electron-electron source-box checkpoint is therefore:

- product/product accepts caller-supplied density-normalized factors and has
  the `ad74d3c` raw-weighted conversion wrapper;
- PQS/product accepts caller-supplied density-normalized factors and has the
  `b1ee2a5` raw-weighted conversion wrapper;
- PQS/PQS accepts caller-supplied density-normalized factors and has the
  `3b64e91` raw-weighted conversion wrapper;
- the private route-shaped consumer at `93a9af8` assembles the
  left-PQS/right-PQS/product-slab retained matrix from the existing
  source-box pair-family helpers;
- the private route producer at `e868f49` builds the same route descriptor
  from explicit source-box fixture facts and returns the descriptor plus the
  density-density consumer result;
- the private IDA provenance extractor at `3028556` records PGDG
  gausslet/source-box density-normalized factors, raw factors, source/raw
  weights, centers, and diagnostics without consuming MWG supplement/residual
  coupling;
- the private adapter at `5de13b1` feeds that IDA provenance object into the
  existing explicit route producer while preserving the same private
  source-box route and retained two-index output convention;
- the private dense-parent authority check at `c443af2` validates the small
  route against `C_route' * V_parent_ida * C_route` to roundoff while keeping
  dense parent projection validation-only;
- the private component smoke at `33a76ed` and `30fd052` combines by-center
  nuclear-attraction component matrices with the IDA source-box retained
  two-index density-density route, including density-normalized and
  raw-weighted electron-electron modes;
- the private component-route report adapter consumes those source-box/PQS
  summaries plus final-residual MWG supplement facts as adjacent report lanes,
  with tracked synthetic coverage for owner, authority-error, raw/fixed-GTO
  MWG-block, and retained-column diagnostic-weight/IDA-division gates;
- all current outputs are retained two-index density-density blocks, not
  four-index Galerkin Coulomb tensors;
- explicit route calls may still use synthetic or caller-supplied data, and
  the private adapter can now use IDA gausslet/source-box provenance from PGDG
  intermediates; MWG supplement/residual provenance is not connected;
- source weights are provenance/positivity checks for density-normalized
  input, or the owner of raw-weight normalization inside the raw-weighted
  wrappers;
- raw-weighted route input divides only by explicit raw/source weight outer
  products and then delegates to density-normalized cores;
- dense parent IDA projection is a diagnostic authority comparison only and
  does not replace the source-box-first route;
- PQS/product and PQS/PQS use no shell projection, Lowdin cleanup, support
  coefficient matrix, or support-local oracle as the algorithm, and the route
  consumer preserves the same boundary;
- retained PQS columns have no positive quadrature-weight or IDA-division
  semantics.

## Private Route-Adapter Migration Contract

The next repo-side migration object should be a private Be/Be2-shaped route
adapter/report, not a public route and not a packet or Hamiltonian builder.
Its job is to collect the adjacent component facts that already exist, check
their boundaries, and state whether the private route is ready for the next
implementation pass.

The intended inputs are:

- explicit source-box retained-unit facts for the fixed/PQS lane: route
  descriptor unit keys, retained ranges, source boxes, source-mode dimensions,
  retained-rule kinds, and supported source-box pair families;
- by-center nuclear-attraction data: source-box nuclear axis layers, Gaussian
  expansion, center coordinates, nuclear charges, and stable center labels;
- IDA electron-electron data: `_pqs_source_box_ida_factor_provenance(...)`,
  term coefficients, pair-factor normalization mode, raw/source weights for
  raw-weighted input, and the source-box route descriptor;
- optional dense-parent IDA authority data for validation only, currently for
  density-normalized rows;
- ordinary final-residual MWG supplement facts from
  `_qwrg_final_residual_mwg_component_blocks(...)` or a report carrying the
  same explicit residual-owner metadata, component shapes, and authority
  errors;
- optional reporting sidecars such as fixed-side retained-unit metadata,
  source-shell/source-mode metadata, and fixed-column/source-relation labels.

The intended output is a private route-adapter record plus a text or TSV
report. It should include:

- route shape, retained dimension, retained-unit records, retained ranges, and
  source-mode dimensions;
- one by-center nuclear retained matrix per nucleus/center, plus a total that
  is explicitly the sum of those center pieces;
- the source-box IDA retained two-index density-density block, with pair
  counts, pair-family counts, helper names, term count, finite/symmetry
  checks, and pair-factor normalization mode;
- final-residual MWG supplement component facts as a separate lane: fixed,
  fixed-residual, residual-residual, and final shapes, explicit residual
  owners, authority error, and raw-GTO/fixed-raw-GTO block absence flags;
- authority-comparison facts and skip reasons;
- timing/allocation facts for component assembly and report adaptation;
- no-go diagnostics for public/default behavior, packet/fixed-block/
  QW/Hamiltonian adoption, MWG/IDA semantic changes,
  retained-column diagnostic-weight/IDA division, owner/shell/ray inference,
  ECP behavior, and CR2 science status.

The fixed/PQS source-box lane and the final-residual MWG lane stay adjacent,
not merged, until a later reviewed adapter explicitly composes them. By-center
nuclear attraction stays separately attributed by center so counterpoise and
center-level diagnostics can inspect the pieces. Source-box IDA output remains
the retained two-index density-density convention, not a four-index Galerkin
Coulomb tensor. Raw-weighted IDA input may be accepted only through the
existing raw/source quadrature-weight normalization wrappers; retained PQS,
source-box, or final-residual columns must not acquire positive quadrature
weight or IDA-division semantics.

Available authority comparisons today are limited:

- source-box safe-term and component smokes can compare against current
  fixed-block or sequence-packet retained fields where those fields exist;
- product/product low-order and kinetic source-box blocks can compare to
  existing product-staged retained kernels;
- density-normalized source-box IDA rows can compare to dense-parent IDA
  projection, `C_route' * V_parent_ida * C_route`, as validation only;
- final-residual MWG supplement facts can require zero reported authority
  error and explicit residual owner metadata from the ordinary MWG lane;
- shell-realized PQS support-local paths remain oracle/debug comparisons and
  are not the source-box algorithm.

Performance reporting is part of the private adapter contract. A complete
smoke should record retained dimensions, source/residual/final dimensions,
unit counts, upper-triangular pair counts, IDA term count, center count,
whether dense diagnostic matrices were used, wall time, allocations, and
GC time when practical. The performance category remains
`diagnostic/prototype` until representative Be/Be2-sized timings show that the
adapter can run without accidental dense-parent or support-row scaling.

A complete private Be/Be2 component-route smoke requires all of the following:

- fixed/PQS source-box retained-unit records are present for the route under
  test;
- by-center nuclear source-box matrices are finite and their explicit sum
  matches the reported total;
- IDA source-box electron-electron rows run for density-normalized input, and
  raw-weighted rows run or report a clear unsupported reason;
- source-box IDA pair counts and pair-family counts cover the intended
  upper-triangular route pairs;
- available authority comparisons pass or record explicit skip reasons;
- final-residual MWG component facts have consistent fixed/residual/final
  dimensions and explicit residual owner metadata;
- no owner, shell, ray, or source label is inferred from centers,
  nearest-grid reconstruction, support order, support indices, or
  `raw_to_final`;
- raw GTO/GTO and fixed/raw-GTO MWG blocks remain absent;
- no retained/source-box/final-residual weight division is introduced;
- the report states that no public/default route, packet/fixed-block/
  QW/Hamiltonian path, ECP behavior, MWG/IDA semantics, or CR2 science status
  changed.

The explicit gate before public/default or packet/QW/Hamiltonian adoption is a
separate review. That review must have a private adapter report at
representative Be/Be2 scale, passing authority comparisons, documented
performance, complete no-go diagnostics, and a manager-approved decision about
how the source-box fixed/PQS lane composes with the ordinary final-residual
MWG supplement lane. Until then, the correct implementation seam is private
report/adaptation around `_pqs_pqs_product_source_box_component_route_smoke(...)`,
`_pqs_pqs_product_component_route_smoke_report_adapter(...)`, and the CR2
sidecar/report writers, not construction or public-route adoption.

Stop implementation if any of these are unclear:

- whether a candidate factor tensor is raw-weighted or density-normalized;
- which raw/source weights own the IDA/MWG division;
- whether the output is a two-index density-density interaction matrix or a
  different Coulomb representation;
- whether a PQS retained column would need to carry a positive quadrature
  weight;
- whether validation requires changing packet/fixed-block/QW/Hamiltonian,
  MWG/IDA semantics, or public/default routing.

Any next implementation blurb should state:

- which `RetainedRule` variant is used on each side;
- whether raw product-box pair operators are materialized or streamed through
  separable factors;
- which safe terms are supported;
- what support-row oracle is used only for validation;
- why no packet/fixed-block/QW/Hamiltonian/default behavior changes.

No packet adoption, fixed-block construction adoption, QW/Hamiltonian routing,
ECP, MWG/interaction implementation, public/default route, or CR2 claim is
part of that correction. Additional local/Gaussian work beyond the private
one-body checkpoint remains separately scoped.

Any next implementation step must still stop before real MWG/IDA factor
provenance beyond the private IDA gausslet/source-box adapter, MWG
supplement/residual adaptation, packet/fixed-block/Hamiltonian adoption,
public/default route behavior, ECP behavior, or CR2 claims.
