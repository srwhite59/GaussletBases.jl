# Raw Product Source And Retained Transform Policy

This note records the intended next abstraction for high-order Cartesian
construction. It is a design contract, not a claim that all routes already use
this machinery.

## Core Idea

Most retained Cartesian units should be understood as:

```text
raw product source space -> retained-space transform -> retained unit
```

The raw product source space is a small local product block built from the
parent gausslet axes. The retained-space transform then performs the local
selection, projection, cleanup, or identity map needed for the construction.

The first reusable sizing object should be the source-box dimension plan:

```text
box intervals + angular-spacing policy + constraints
-> source_mode_dims = (nx, ny, nz)
-> deterministic 1D source transforms
-> retained rule
```

The `source_mode_dims` are total source-mode lengths. They are not interior
lengths and should not be described as "retain plus two" in PQS-facing
contracts. Route-local selector counts may still exist as diagnostics, but the
source-box contract should make the total dimensions primary. For bond-aligned
diatomic boxes, `L` is simply the source-mode dimension on the bond axis after
the side-dimension policy has been applied.

Operator blocks should be assembled in the raw product source spaces first,
then finished by small dense transforms:

```text
O_retained = T_left' * O_raw_product * T_right
```

This is the preferred framing for ordinary product slabs, midpoint caps,
rectangular shells, and projected q-shells. It should reduce the number of
special mixed-block cases because each side presents the same shape: a raw
product source plus a retained transform.

## Examples

A midpoint slab `q x q x 1` is the simple case. Its raw source is the product
space

```text
X_q * Y_q * fixed_z
```

and its retained transform is usually identity or close to identity. This is a
special case of a product-source retained unit, not a separate Hamiltonian
idea.

A projected q-shell, `PQS(q, L)`, starts from the full local product block

```text
X_q * Y_q * Z_L
```

The shell is not defined by subtracting a contracted inner block. There are two
separate stages.

The raw product-box stage is:

1. build the full local product transform;
2. select boundary COMX-product modes whose local product-mode index is first
   or last on at least one axis;
3. build raw-box operator references by selecting those same product-box
   columns on the left and right.

The shell-realization stage is:

1. project the selected product-box modes onto shell rows;
2. apply full-rank symmetric Lowdin cleanup;
3. use the resulting isometric shell-supported representation.

For a rectangular shell this gives the counting picture
`q*q*L - (q-2)*(q-2)*(L-2)`, but that count is not the construction rule. The
construction rule starts with boundary product-mode selection from the full
product block. Shell-row projection and cleanup are a later realization step.
This avoids any dependence on how an interior cube would have been contracted.

## Operator Consequence

Mixed blocks such as PQS/product should not be treated as strange special
cases. Both sides should expose:

- raw product source metadata;
- a retained transform;
- support and boundary ownership diagnostics;
- metric/operator capability labels.

Then overlap, position, kinetic, `x2`, and Gaussian local one-body blocks can
be formed in the raw product source spaces and transformed down. For
factorized two-electron or MWG-style terms, the same philosophy may apply, but
the contraction must go through the factor terms carefully. It must not assume
positive quadrature weights for functions that are not quadrature-carrying.

## Raw Source Pair Operator Packets

A real retained basis has several retained units: atom cores, midpoint slabs,
PQS shells, contact caps, mismatch shells, supplement blocks, and fallback
pieces. The final operator matrix is assembled from blocks between every pair
of retained units.

The missing intermediate object is therefore a pair packet for raw product
source spaces:

```text
source_i, source_j
-> RawProductSourcePairOperatorPacket3D
-> T_i' * O_raw(i,j) * T_j
-> retained operator block
```

Only the upper triangle needs to be built for symmetric operators. The packet
should be lazy or cacheable so large constructions do not materialize all raw
pairs too early.

A useful vocabulary is:

```text
RawProductSource3D
    parent axes
    local box
    axis transforms
    source dimension
    support map / boundary map
    raw source quadrature weights

RetainedTransform3D
    source id
    retained dimension
    transform T
    selection / projection / cleanup diagnostics
    quadrature weight role for retained columns

RawProductSourcePairOperatorPacket3D
    left source id
    right source id
    operator kind
    raw overlap / position / kinetic / gaussian factors / pair factors
    backend provenance
    supported term list
    symmetry status
```

This is the structure that makes a PQS/product mixed block ordinary. The PQS
side has a full `q x q x L` product source and a raw-box retained rule that
selects boundary COMX-product modes. The product or midpoint-slab side has its
own product source and usually a simple retained transform. The operator block
is built between the two raw sources, then transformed on both sides.

Shell realization is separate from this raw-box rule. When a PQS object must be
realized on shell rows, the selected boundary modes are projected to shell rows
and followed by Lowdin cleanup. PQS/product source-box operators use the
raw-box boundary-mode rule, not that shell-realization map.

## Private Source-Box Migration Plan

The next implementation line should make the current private helpers converge
on three explicit records. These are private planning/reference records first,
not public API.

`CartesianRawProductBox3D` should own the source-box facts:

- parent axis intervals;
- total source-mode dimensions `(nx, ny, nz)`;
- 1D COMX/source transforms for each axis;
- source-mode ordering and product-mode count;
- optional 1D operator factors for overlap, position, `x2`, and kinetic;
- coordinate/backend/provenance diagnostics;
- raw-source quadrature-weight role, when relevant.

It should not own retained column selection, shell projection, Lowdin cleanup,
packet fields, or QW/Hamiltonian routing. Given intervals and
`source_mode_dims`, raw product-box construction should be deterministic.

`CartesianRetainedRule3D` or `CartesianRetainedUnit3D` should own the mapping
from a raw source box to retained columns:

- source-box id or embedded raw-box plan;
- retained rule kind such as `:identity`, `:product_doside`,
  `:pqs_boundary_mode_selection`, `:shell_projection_lowdin`, or
  `:support_dense_fallback`;
- retained dimension and column range, once assigned;
- retained transform metadata, including whether the transform is separable,
  selected columns, dense/materialized, or factored;
- cleanup/isometry diagnostics when the rule uses projection plus Lowdin;
- retained-column weight semantics and IDA/MWG division permissions.

For the corrected mode-selected PQS source-box path, the retained rule is
boundary COMX-product mode selection. Shell-row projection plus Lowdin is a
separate realization rule, not part of raw-box operator construction.

`CartesianSourceBoxPairOperatorPlan3D` should own pair-level operator facts:

- left/right raw source boxes and retained rules;
- axis interval compatibility and 1D cross-axis factors;
- supported terms and symmetry status;
- raw pair operator provenance and backend diagnostics;
- whether the retained block can be formed by
  `T_left' * O_raw_box_pair * T_right`;
- reference/shadow/adoption status.

This pair plan is the place where product/product, PQS/product,
support-dense fallback, and future GTO cross-overlap consumers should meet.
It should not imply all-pairs production assembly or packet adoption.

Current helpers can migrate as adapters:

- `_pqs_raw_product_box_plan(...)` is the present PQS-shaped
  `CartesianRawProductBox3D` adapter.
- `_pqs_shell_realization_plan(...)` is the present shell-projection/Lowdin
  retained-rule adapter.
- `_pqs_product_box_realization_plan(...)` should remain only a wrapper that
  returns both plans together for diagnostics and compatibility.
- `_pqs_raw_product_box_reference_block(...)` should consume the raw-box plan
  directly for PQS/PQS self references.
- `_pqs_product_source_box_pair_plan(...)` and
  `_pqs_product_source_box_reference_block(...)` should consume raw-box and
  retained-rule facts, not descriptor shell rows.
- `_product_doside_retained_low_order_block(...)`,
  `_product_doside_retained_kinetic_block(...)`, and product-staged metric
  helpers are existing product-side adapters into the same retained-unit
  vocabulary.
- `_pqs_source_box_gto_cross_overlap_shadow(...)` should stay a raw-box
  consumer and not become the final-basis handoff authority.

The shared private helper `_cartesian_raw_product_box_plan(...)` is now the
current `CartesianRawProductBox3D`-shaped source-box fact owner. It wraps the
1D source-axis transform planning and records source-box intervals, total
source-mode dimensions, z-fast source-mode ordering, axis transforms,
axis-local coefficients, and retained-rule-free diagnostics. PQS raw plans
wrap and validate this shared plan when it is supplied. Descriptor-only PQS
plans remain structural adapters and explicitly report the shared plan
unavailable, because descriptors do not carry axis bundles and should not fake
bundle-derived source transforms.

Private PQS/product and PQS/GTO source-box shadows now exercise shared-backed
PQS raw plans where those shared plans are available. This proves consumer
compatibility with the shared source-box seam while preserving the old
authority boundary: the shared raw-box plan does not drive packet construction,
fixed-block construction, QW/Hamiltonian assembly, IDA/MWG semantics,
shell-realization adoption, or any public/default route.

The latest cleanup moved private source-box consumers closer to that boundary.
PQS raw-box operator helpers now assemble through raw-plan methods; descriptor
level PQS raw-box reference remains only a convenience adapter that builds a
raw plan and delegates. The PQS/product shadow layout is raw-plan-first, with
the descriptor method kept only as a compatibility wrapper that checks
descriptor/raw-plan consistency. The PQS/GTO source-box shadow is also
raw-plan-first, with the descriptor method acting only as an adapter. The
final-basis GTO handoff remains separate and authoritative.

The product side now has a matching private retained-unit metadata adapter.
`_product_doside_retained_unit_plan(...)` records product/source axis
intervals, physical/source axis lengths, retained axis counts, column range and
retained count, axis coefficient matrices, and `axis_function_indices` for an
existing product/doside staged unit. It is metadata-only: it does not rebuild
coefficients, does not change product/doside retained-block math, does not
adopt packet construction, and does not change IDA or retained-weight
semantics. `_pqs_product_source_box_pair_plan(...)` now obtains its
product-side metadata through this adapter, so the private pair plan reads as
raw product-box source plus retained-rule facts on both sides. This is still
only migration/shadow infrastructure, not a generic retained-unit framework.

Product/product source-box unification is now a completed private/shadow
checkpoint, not a future design choice. Product/product now uses the same
private source-box pair vocabulary.
`_product_doside_source_box_pair_plan(...)` builds the product/product pair
metadata from two retained-unit plans, records 1D cross factors, and keeps
existing product-staged helpers authoritative. `_product_doside_source_box_reference_block(...)`
supports overlap, position, `x2`, and kinetic; it compares each retained block
against `_product_doside_retained_low_order_block(...)` or
`_product_doside_retained_kinetic_block(...)`. `_product_doside_source_box_shadow_blocks(...)`
builds a small two-block product/product shadow layout from those reference
blocks and checks transpose consistency for symmetric real terms. The shared
separable term descriptor is now
`_source_box_separable_term_factor_kinds(...)`, used by both product/product
and PQS/product source-box blocks. Existing product-staged helpers remain the
numerical authority; this checkpoint does not adopt packet construction,
fixed-block construction, QW/Hamiltonian assembly, IDA/MWG, or public/default
routes.

A focused provenance checkpoint now keeps fallback reporting separate from
input-operator provenance. `_cartesian_raw_product_box_plan(...)` records
`integration_contract = :pgdg_exact` and `numerical_reference_fallback = false`
for PGDG-backed test fixtures. PQS raw-box self references and PQS/product
source-box references forward that raw-box fallback status. Product/product
source-box references report `operator_factor_source = :explicit_metric_operator_data`,
`input_metric_operator_data = :caller_supplied_explicit_data`,
`input_metric_operator_data_pgdg_checked = false`, and
`numerical_reference_fallback = false`. For product/product, the fallback flag
means only that the helper did not invoke a numerical-reference fallback; it is
not a claim that the caller-supplied matrices have been proven PGDG analytic
inside the helper. The PQS/GTO source-box shadow likewise records that it does
not use a numerical-reference fallback. This is provenance reporting only; it
does not change backend selection or quadrature policy.

After this raw product-box migration loop, the private state is:

- the shared raw product-box plan exists;
- PQS raw plans wrap and validate the shared plan when axis bundles are
  available;
- PQS self, PQS/product, shadow-layout, and PQS/GTO source-box shadows are
  raw-plan first;
- descriptor-level methods are compatibility adapters;
- product/doside retained metadata is exposed through a retained-unit adapter;
- product/product source-box pair, reference-block, and two-block shadow
  helpers exist privately;
- no packet, fixed-block, QW/Hamiltonian, IDA/MWG, public, or default-route
  adoption has happened.

Likely next design choices are deliberately separate: clean up
shell-realization consumers that use projection plus Lowdin, optimize
PQS/product retained blocks, or extend the vocabulary to support-dense
fallback. None of those should be folded into this closeout checkpoint.

During migration, these paths remain authoritative:

- `_nested_shell_packet(...)` remains authoritative for active packet matrix
  construction.
- Existing fixed-block, QW, Hamiltonian, IDA/MWG, supplement, and public route
  consumers remain authoritative.
- Existing final-basis GTO handoff remains authoritative for handoff use.
- The new raw-box plans and pair plans are private reference/shadow evidence
  until an explicit adoption pass changes that.

Validation should be staged and mechanical:

- PQS self: raw plan plus boundary selector reproduces overlap, position,
  `x2`, and kinetic source-box references for cubic and rectangular boxes.
- PQS/product: raw-box pair factors plus product retained transform reproduce
  explicit source-box references, including a non-identity product transform.
- Product/product: source-box reference and shadow helpers compare against
  existing product/doside retained low-order and kinetic helpers, which
  continue to match product-staged metric and kinetic references.
- Support-dense fallback: mixed product/support and support/support paths
  remain explicit fallback references, not silently optimized kernels.
- GTO cross overlap: raw-box/GTO shadow matches dense source-box reference and
  remains separate from final-basis handoff.
- Diagnostics must state which stage is active: raw-box reference,
  shell-realization reference, shadow packet/layout, or construction adoption.

Non-goals for this migration stage:

- no public API or default-route promotion;
- no QW/Hamiltonian construction adoption;
- no all-pairs production packet builder;
- no local/ECP/Gaussian/MWG/interaction implementation;
- no retained PQS weight division or positive retained-weight IDA claim;
- no assumption that shell-row projection plus Lowdin is part of raw-box
  operator construction;
- no broad generic framework layer unless it removes a concrete duplicated
  helper path.

## Weight Contract

Retained units need an explicit weight role.

The safest positive quadrature object is the raw product source. Raw
gausslet-derived product sources should have positive source-point weights in
non-pathological cases:

- product slabs;
- midpoint caps treated as product slabs;
- atom-local cubes and shells;
- rectangular product shells;
- projected q-shells.

For these raw source spaces, negative or zero source weights in normal geometry
are a construction problem, conditioning problem, or gauge problem.

Retained-column weights need a stricter label. If the retained transform is
identity-like, retained weights may remain positive. If the retained transform
includes projection and Lowdin cleanup, as in PQS, induced per-column weights
can become small or signed depending on the gauge and conditioning. The
intended PQS construction should induce positive retained weights in ordinary
cases, but code must not silently rely on that fact. IDA/MWG-style paths that
divide by retained weights should either operate at the raw source level or
require an explicit retained-column `:positive_required` role plus a finite
positive check.

This label is independent of whether an operator object contains the intended
gausslet IDA/MWG electron-electron interaction. For PQS, the active fixed/final
interaction may still be the normal nested pair-sum IDA/MWG route, while the
private PQS retained-transform payload reports retained-column weights as
`:debug_reference_only` with `ida_weight_division_allowed = false`. That
combination is intentional: the interaction path remains active, but retained
PQS columns are not themselves promoted to positive quadrature carriers.

Supplement or final-basis functions are different. Angular GTO supplements
such as p and d functions, residual directions, and other non-quadrature
additions do not carry positive quadrature weights. Their integral or weight
can be zero, signed, or not meaningful as a quadrature mass. Nothing in the
construction should divide by those weights, and IDA-style positive-weight
paths must be forbidden for those functions.

A useful payload field would be:

```text
quadrature_weight_role =
    :raw_source_positive
    :positive_required
    :not_quadrature_weight
    :debug_reference_only
```

The intended meaning is:

- `:raw_source_positive`: finite positive weights are attached to raw product
  source points and can be used by raw-source quadrature contractions.
- `:positive_required`: finite positive retained-column weights are part of
  the construction contract before IDA/MWG-like use.
- `:not_quadrature_weight`: weights are not valid quadrature masses; reject
  any path that divides by them.
- `:debug_reference_only`: the payload is a reference/debug object and should
  not be promoted to production quadrature behavior without a separate review.

## Fit With The Current Framework

The existing framework direction is:

```text
Cartesian parent
-> shell/region descriptor
-> contraction rule
-> resolved payload
-> metric/operator helper
```

This note refines what an execution-ready payload should contain. A useful
payload should not merely say "PQS" or "product". It should expose the raw
product source and retained transform that let operator code build the raw
block first and finish with small linear algebra.

For pairwise operator construction, the next useful object is a
`RawProductSourcePairOperatorPacket3D`-style packet. That packet should be the
place where raw product-source operator factors, backend provenance, symmetry,
and supported-term lists are recorded before retained transforms are applied.

Current private PQS sidecar fixtures are still low-order/reference-scoped.
They prove payload discovery and low-order checks, but they do not yet install
this full raw-product-source contract into production metric, QW, Hamiltonian,
or CR2 routes.

## Private Metadata Checkpoint

The current private implementation now has the first raw-source planning
scaffold:

- raw product source metadata;
- retained transform metadata;
- an upper-triangular raw source pair plan;
- a resolved raw source pair audit;
- a private raw overlap packet;
- a private toy raw `:axis_index_x` packet for the identity product/slab
  fixture;
- private physical raw `:position_x`, `:position_y`, and `:position_z` packets
  for product/slab fixtures, consuming explicit 1D axis metric data with
  non-integer positions;
- private retained low-order blocks for materialized product/slab transforms;
- a nonidentity materialized product/slab retained transform;
- product/product cross-pair physical position packets for same-support
  identity/nonidentity fixtures;
- product/product cross-pair physical position packets for distinct shifted
  supports with nontrivial explicit axis metrics;
- private factorized product/doside raw packets for `:overlap` and physical
  `:position_x`, `:position_y`, and `:position_z`, built from explicit 1D axis
  factors over staged local intervals;
- private retained-block comparisons against the existing product-staged metric
  block kernel for overlap and physical position terms.

The retained-block calculations are intentionally tiny product/slab fixtures.
The baseline commits for the raw-packet line include `8ba7ab4 Rename toy
product packet axis index term`, `d0fb995 Add private physical product position
packets`, `3ef0fa1 Generalize private physical product packets`, `9589f9a Test
nonidentity product retained position blocks`, `26f88d6 Add private product
cross position packets`, and `c6d4f41 Test product cross packets on distinct
supports`. The factorized packet checkpoint continues through `107e68c Add
private factorized product overlap packets`, `3c1eb95 Add private factorized
product position packets`, and `1c28f8f Compare factorized product packets to
staged metrics`.

```text
source dimension 4 -> retained dimension 4
retained transform kind = :product_axis_transform
T' * I_raw * T == I_4
T' * axis_index_x_raw * T == diag(1, 1, 2, 2)
T' * position_x_raw * T == diag(0.25, 0.25, 1.75, 1.75)
T' * position_y_raw * T == diag(-0.5, 0.5, -0.5, 0.5)
T' * position_z_raw * T == diag(3.25, 3.25, 3.25, 3.25)
```

The `:axis_index_x` packet is a toy axis-index diagnostic for the fixture. It
is not a physical `position_x` operator and should not be extended as if it
were one.

The private physical position packets are separate. They build raw matrices
over raw product source support rows from explicit 1D axis metric data and do
not build dense full-parent 3D matrices. Retained blocks are checked only
where the retained transforms are materialized. The current private checks
cover self-pair retained overlap/position, nonidentity retained transforms,
same-support product/product cross pairs, and distinct-support product/product
cross pairs. All of these remain fixture/reference checks, not production
metric execution.

The private factorized product/doside packets are the next validation layer.
They build materialized raw support-row matrices for these small fixtures, but
the entries are assembled from explicit 1D overlap/position factors restricted
to staged local intervals. Retained blocks are then checked with the same
`T_left' * O_raw * T_right` algebra and compared against
`_fill_product_staged_metric_blocks!(...)`, the existing product-staged
retained metric block kernel. This proves the low-order product/product
factorization for the private fixtures; it is still not a large-region fast
operator assembly path.

The follow-up retained-block checkpoint removes that raw support-row matrix
from the private product/doside comparison path. The helper
`CartesianContractedParentMetrics._product_doside_retained_low_order_block(...)`
directly returns retained product/product blocks from staged 1D factors for
`:overlap` and physical `:position_x`, `:position_y`, and `:position_z`. It
matches both `_fill_product_staged_metric_blocks!(...)` and the private
raw-packet retained path on the focused fixtures. Guard tests pin that
unsupported terms reject, non-`:product_doside` units reject, and axis metadata
or axis metric mismatches reject. This checkpoint is recorded by
`bde3dac Add private product retained block helper` and
`0470f27 Test product retained block helper guards`.

The product low-order adoption checkpoint then uses these private helpers only
inside the existing product/doside low-order metric plumbing. As of
`250b3d7 Adopt product retained block fill helper`,
`_fill_product_staged_metric_blocks!(...)` delegates product/doside
`:overlap` and physical `:position_x`, `:position_y`, and `:position_z`
retained block fills to `_product_doside_retained_low_order_block(...)`. As of
`57738aa Factor product linear vector helper`,
`_staged_unit_linear_vectors(...)` delegates product/doside retained weights
and first moments to `_product_doside_retained_linear_vectors(...)`.
Support-dense and generic fallback paths remain unchanged. Pair blocks and
one-unit linear vectors remain distinct concepts; weights and first moments
are not pair packets.

The helper seam is still private low-order infrastructure. Its adoption only
removes duplicate product/doside internals in the existing metric path. It does
not change public/default routes, QW/Hamiltonian builders, backend/default
policy, PGDG/quadrature policy, IDA or positive-weight assumptions, CR2 paths,
or science validation. It also does not add PQS/product mixed blocks,
support-dense mixed packet generalization, all-pairs production matrices,
kinetic, nuclear/local, Gaussian, MWG, or interaction terms. GTO supplement
functions remain outside retained-column positive quadrature assumptions.

The first private non-low-order operator checkpoint is kinetic. As of
`04b7095 Add private product kinetic block helper`,
`_product_doside_retained_separable_sum_block(left_unit, right_unit,
axis_factor_terms)` builds retained product/doside blocks from explicit
separable axis-factor triples, and
`_product_doside_retained_kinetic_block(left_unit, right_unit, axis_ops)` is a
thin wrapper for

```text
(Kx, Sy, Sz) + (Sx, Ky, Sz) + (Sx, Sy, Kz)
```

The focused tests compare product/doside self and cross-unit kinetic blocks
against support-local retained references formed from `_staged_unit_entries(...)`
and `_contract_pair_block(...)`. This is private signed operator-block
infrastructure only. Kinetic is not a weight object, not a positive quadrature
claim, and not an IDA division path.

The real-block shadow checkpoint is recorded by
`ecad56a Test product kinetic helper on real blocks`. The existing q4
bond-aligned endcap-panel shared-shell source-policy fixture exposes six
`:product_doside` units through its staged sidecar. The test compares
`_product_doside_retained_kinetic_block(...)` against the already-built real
route matrix block
`fixed_block.kinetic[left.column_range, right.column_range]` for all 21
upper-triangular product/product pairs. Cross-pair mirror blocks are checked
against the transposed kinetic block, and the observed maximum block error was
about `8.9e-16`. This remains private shadow validation only, not construction
adoption.

The product-only shadow-matrix checkpoint is recorded by
`8cd2f17 Add product-only kinetic shadow matrix`. The private helper
`_product_doside_retained_kinetic_shadow_matrix(...)` fills only the
product/product kinetic subblocks using
`_product_doside_retained_kinetic_block(...)`. In the same q4 bond-aligned
endcap-panel fixture, it covers six product units and 21 product/product
blocks, matching `fixed_block.kinetic` product subblocks to about `8.9e-16`.
Rows and columns outside the product units remain explicitly zero/absent, so
the helper does not claim support-dense or mixed kinetic coverage.

The full private kinetic shadow-matrix checkpoint is recorded by
`59c2595 Add full kinetic shadow matrix`. The private helper
`_staged_retained_kinetic_shadow_matrix(...)` fills product/product blocks
with `_product_doside_retained_kinetic_block(...)` and fills every other
staged-unit pair with support-local kinetic fallback through
`_staged_unit_entries(...)` and `_contract_pair_block(...)`. In the same q4
bond-aligned endcap-panel fixture, it checks the complete shadow matrix
against `fixed_block.kinetic`: 21 product/product blocks, 7 fallback blocks,
28 total upper-triangular blocks, and maximum full-matrix error about
`6.6e-14`. This is private shadow/reference infrastructure only, not kinetic
construction adoption.

The packet-build source plan checkpoint is recorded by
`c78e9e4 Add packet build source plan metadata` and
`6beebad Tighten packet build source metadata`. The private records
`_CartesianPacketBuildSource3D` and `_CartesianPacketBuildPlan3D`, with helpers
such as `_cartesian_packet_build_source(...)` and
`_cartesian_packet_build_plan(...)`, describe the resolved payload layout that a
future packet builder might consume. The tightened contract is deliberately
metadata-only: the fields are candidate packet fields, not executable operator
fields; the readiness flag is column-layout readiness, not operator readiness;
operator data is explicitly unavailable and unchecked; no numerical packet
matrices are built; and `_nested_shell_packet(...)` remains the authoritative
current packet builder.

The q4 bond-aligned endcap-panel fixture is the current structural evidence for
that plan shape. It records parent dimension `539`, contracted dimension `313`,
six product/doside payloads, support-dense fallback payloads, and complete
column coverage. Its candidate packet fields are `:overlap`, `:position_x`,
`:position_y`, `:position_z`, `:weights`, `:first_moments`, and `:kinetic`.
Missing or not-implemented packet fields remain explicit: `:x2_x`, `:x2_y`,
`:x2_z`, `:nuclear_one_body`, `:local_coulomb_one_body`,
`:local_ecp_one_body`, `:gaussian_local_terms`, `:gaussian_sum`, `:pair_sum`,
`:mwg_interaction`, and `:interaction`. The support union summary is
informational only; parent support completeness is not required, and overlapping
payload support is allowed.

This plan does not drive construction. It changes no fixed-block construction,
contracted-parent metric packet execution, QW/Hamiltonian path,
backend/default policy, PGDG/quadrature policy, IDA or positive-weight
semantics, CR2 path, or science behavior. Before more construction-boundary
coding, repo-auditor/user review should decide whether and how a future packet
builder should consume this source shape upstream, rather than reconstructing
it after `_nested_shell_packet(...)` has already built the authoritative
packet.

The first safe-field shadow checkpoint is recorded by
`b662efa Add packet build source safe-field shadow`. The private helper
`CartesianContractedParentMetrics._cartesian_packet_build_source_safe_field_shadow(...)`
consumes a `_CartesianPacketBuildSource3D` plus explicit axis data and
reconstructs the safe field set only:

```text
overlap
position_x / position_y / position_z
weights
first_moments
kinetic
```

On the q4 bond-aligned endcap-panel fixture, the shadow fields match the
authoritative/reference data with maximum errors:

```text
overlap       1.71e-15
position_x    5.33e-15
position_y    5.33e-15
position_z    1.00e-14
weights       3.51e-14
first_moments 0.0
kinetic       6.57e-14
```

`first_moments` are compared against the contracted-parent metric packet
reference, not `_CartesianNestedShellPacket3D`, because the shell packet does
not carry first moments. Product/product blocks use the retained
product/doside helpers; support/support and product/support blocks use the
support-local fallback.

This helper is still source-driven shadow infrastructure, not packet
construction adoption. `_nested_shell_packet(...)` remains the current
numerical packet authority, and no fixed-block construction, metric-packet
public behavior, QW/Hamiltonian path, backend/quadrature policy, IDA or
positive-weight semantics, CR2 path, or science behavior changed. The next
architecture question is whether and how to produce `_CartesianPacketBuildSource3D`
upstream before `_nested_shell_packet(...)`, rather than reconstructing it from
a sidecar after packet construction has already happened.

The first pre-packet source checkpoint is recorded by
`bdea72b Add endcap panel pre-packet source shadow`. The private helper
`_cartesian_endcap_panel_pre_packet_build_source(...)` consumes endcap/panel
`owned_units`, `coefficient_matrix`, `unit_column_ranges`, `support_indices`,
and parent `dims`, and produces a `_CartesianPacketBuildSource3D`-shaped
source before packet matrix construction. This is diagnostic/layout-only:
`_nested_shell_packet(...)` remains authoritative.

The helper diagnostics record:

```text
packet_construction_consumes_source = false
source_object_builds_packet_matrices = false
nested_shell_packet_remains_authoritative = true
```

The helper currently uses `coefficient_matrix` only for dimensions and column
coverage; coefficient values are not verified by this pre-packet source
checkpoint, so its effective contract is `coefficient_values_checked = false`.
The focused test calls the helper after the q4 endcap/panel `layer` object
exists, but it passes only ingredients that are already available before packet
construction: owned units, coefficient matrix, unit column ranges, support
indices, and dimensions.

On the q4 endcap/panel fixture, the pre-packet source matches the corresponding
post-sidecar source in dimensions, payload kinds and counts, column ranges,
support indices, and candidate/missing fields. Its safe-field shadow matches
the authoritative layer packet fields for overlap, position, weights, and
kinetic, while first moments are compared against the contracted-parent metric
packet reference. This is not sequence-level source composition and not packet
construction adoption.

The sequence-level pre-packet source checkpoint is recorded by
`10bd5fe Add sequence pre-packet source shadow`. The private helper
`_cartesian_nested_sequence_pre_packet_build_source(...)` builds a
sequence-level `_CartesianPacketBuildSource3D` from sequence ingredients before
packet construction, again in diagnostic/shadow mode only. Its explicit inputs
are parent `dims`, `core_indices`, `core_coefficients`, `core_column_range`,
`shell_layers`, `layer_column_ranges`, the full sequence `coefficient_matrix`,
and sequence `support_indices`.

The sequence representation follows existing staged-sidecar conventions rather
than inventing new coefficient semantics. The core is represented as one
`:support_dense` payload. Product endcap/panel shell layers expand their
`owned_units` into `:product_doside` payloads with sequence-global column
ranges. Non-product layers use support-dense fallback through
`coefficient_matrix[:, layer_range]`.

The helper diagnostics preserve the same authority boundary:

```text
packet_construction_consumes_source = false
source_object_builds_packet_matrices = false
nested_shell_packet_remains_authoritative = true
```

Two auditor caveats are part of the contract. First, `core_coefficients` are
accepted separately from `coefficient_matrix[:, core_column_range]`; the helper
checks shape and support but does not currently check coefficient-value
equality. Second, product/doside layer payloads are reconstructed from
`owned_units`, not from the layer coefficient slice. That matches the sidecar
logic, while the safe-field shadow comparison is the check that catches packet
or coefficient mismatches. This is acceptable diagnostic infrastructure, not
packet construction adoption.

On the q4 bond-aligned endcap/panel fixture, the pre-sequence source matches
the post-sidecar source in dimensions, payload counts and kinds, column ranges,
support indices, coverage summaries, and candidate/missing fields. Its
safe-field shadow matches the authoritative fixed-block overlap,
`position_x`, `position_y`, `position_z`, weights, and kinetic fields. First
moments compare against the contracted-parent metric packet reference.

The private PQS/product reference checkpoint is recorded by
`157c35f Add PQS product low-order reference block` and
`93ddb3b Add PQS product kinetic reference block`. Private reference coverage
now exists for PQS/product overlap, `position_x`, `position_y`, `position_z`,
and kinetic retained blocks. These helpers reconstruct PQS support
coefficients as

```text
seed * cleanup_transform
```

where `seed` is derived from the projected q-shell descriptor, then compare
the reconstructed coefficients to the stored support-local PQS coefficients.
The factored PQS/product block is then compared against a separately built
support-local oracle. The kinetic reference is a signed operator check using

```text
(K, S, S) + (S, K, S) + (S, S, K)
```

Reverse product/PQS blocks are transpose/reference behavior only for real
symmetric blocks. They should not be generalized to nonsymmetric
derivative-like operators without a separate contract and test.

This checkpoint is not an optimized PQS/product kernel, packet or metric
adoption, fixed-block sidecar installation, QW/Hamiltonian path,
public/default behavior, CR2 claim, or retained PQS positive-weight/IDA
semantics. Local Coulomb, local ECP, Gaussian/local terms, MWG/interaction,
and all-pairs production assembly remain out of scope.

This is a contract checkpoint, not production metric execution. The PQS raw
overlap packet exists only as raw packet/reference plumbing. PQS retained
transforms are applied only in the explicit private PQS/product reference
helpers described above, and the PQS Lowdin cleanup matrix is not treated as
the full raw-to-retained transform. Production/optimized/adopted PQS/product
mixed retained blocks remain unsupported. Production support-dense/mixed
packets, production product/product operator assembly, all-pairs matrix
construction, kinetic adoption outside the private reference/shadow helpers,
nuclear/local, Gaussian, MWG, interaction, QW/Hamiltonian, public/default,
backend/default, PGDG/quadrature, CR2, and science/energy behavior remain
unchanged and outside this private fixture line.

The older support-local PQS fixture still needs an explicit way to represent or
resolve its shell-realized support-row transform:

```text
raw_product_modes
-> raw_boundary_projection
-> full_rank_symmetric_lowdin_cleanup
-> retained_columns
```

That factored projection-plus-Lowdin path is the shell-realized/support-local
fixture contract. It is not the retained transform used by the newer raw-box
PQS/product source-box operators, where `T_PQS` is boundary COMX-product mode
selection only.

The private raw product-box self-block checkpoint corrects the earlier
contract boundary. `_pqs_raw_product_box_reference_block(...)` stays entirely
in the mode-selected `q x q x L` product-box space. It supports only overlap,
`position_x`, `position_y`, `position_z`, `x2_x`, `x2_y`, `x2_z`, and kinetic.
The reference oracle is explicit product-box column selection:

```text
O_boundary = P_boundary' * O_product_box * P_boundary
```

where `P_boundary` selects the boundary COMX-product mode columns. The helper
records `max_1d_source_overlap_error`, `max_product_overlap_error`, and
`selected_overlap_error`, and reports that shell projection is postponed. It
does not use shell-row `support_coefficient_matrix` data, does not apply
Lowdin at the raw-box stage, and does not assign positive quadrature or IDA
division semantics to retained PQS columns.

The private `_pqs_product_box_realization_plan(...)` helper now builds the raw
product-box plan and shell-realization plan at the same setup point while
keeping them as separate objects. Its raw plan carries 1D operator factors and
the boundary selector. Its shell plan carries the shell projection matrix,
Lowdin cleanup, and isometry diagnostics. Operator references still use the
raw product-box 1D factor path first; shell projection is not part of raw-box
operator construction.

The raw product-box plan has since been split into the explicit private helper
`_pqs_raw_product_box_plan(...)`. PQS self blocks, PQS/product source-box
blocks, and PQS/GTO cross-overlap shadows now take that raw plan directly.
`_pqs_shell_realization_plan(...)` carries only the later shell-row projection
and Lowdin isometry, while `_pqs_product_box_realization_plan(...)` remains a
wrapper for diagnostics and compatibility. This makes the separation concrete
in code rather than only in prose: source-box pair/operator contractions are
1D-factor raw-box operations, and shell realization is a distinct consumer.

The first PQS/product source-box mixed-block checkpoint now also covers a
nontrivial product retained transform. A focused private test uses a
non-identity product-axis coefficient transform and checks overlap, one
position block, one `x2` block, and kinetic against the explicit source-box
reference with that product transform applied. This keeps the contract at
`T_PQS' * O_raw_box_pair * T_product`; it does not use shell projection,
Lowdin, support-local PQS coefficients, retained PQS weights, or IDA division.

The private `_pqs_product_source_box_shadow_blocks(...)` checkpoint is the
first block-layout consumer of those references. It builds a small two-block
shadow layout containing one mode-selected PQS source-box unit and one
product/doside retained unit, then fills PQS/PQS, PQS/product, product/PQS by
transpose for symmetric real terms, and product/product blocks. The supported
terms are `:overlap`, `:position_x/y/z`, `:x2_x/y/z`, and `:kinetic`. Focused
tests cover a rectangular PQS source box and a non-identity product/doside
transform. This is still private/shadow-only layout evidence, not a packet
builder: it does not use shell-row projection, Lowdin,
`support_coefficient_matrix` as a PQS oracle, retained PQS weights, or IDA
division, and it changes no packet construction, QW/Hamiltonian,
public/default, CR2, local/ECP/Gaussian/MWG/interaction, or IDA/MWG behavior.

The private PQS source-box to GTO cross-overlap checkpoint adds
`_pqs_source_box_gto_cross_overlap_shadow(...)` and
`_pqs_source_box_gto_axis_projection(...)`. The helper keeps the
mode-selected raw product-box contract: it uses PQS source-box intervals,
source coefficients, and boundary COMX-product mode indices, projects existing
1D Cartesian/GTO primitive-axis overlap tables through the source-box axis
coefficients, and returns a boundary-mode by GTO overlap shadow. The
rectangular `5 x 5 x 7` test compares against an explicit dense source-box
reference built from parent-row/GTO overlap followed by
`transpose(C_source_box) * S_parent_gto` and boundary-column selection.

This is not final-basis handoff adoption. The existing final-basis GTO handoff
remains unchanged and authoritative for handoff use. The shadow helper does
not use shell-row projection, Lowdin, `support_coefficient_matrix` as a PQS
oracle, retained PQS weights, retained-weight IDA division, packet
construction, QW/Hamiltonian paths, public/default routes, CR2 artifacts,
local/ECP/Gaussian/MWG/interaction paths, or IDA/MWG behavior changes.

When shell rows are needed, the realization consumer should remain separate:

```text
selected_product_box_modes
-> project_to_shell_rows
-> Lowdin_cleanup
-> isometric_shell_realized_representation
```

That shell-realization helper should consume the raw-box result or equivalent
selected product-box modes by conjugation. It should not be folded into the
raw product-box reference path.

General physical raw source pair operator packets beyond these private
product/slab fixtures, kinetic adoption into existing real route consumers,
nuclear, Gaussian/local, interaction/MWG, optimized/adopted product/PQS
pairs, and product/support-dense pairs need separate design. The current
checkpoint changes no all-pairs matrix construction, QW or Hamiltonian path,
public/default route, backend/default policy, PGDG or quadrature policy, CR2
path, or science status.

## Non-Goals

This policy does not change:

- public/default construction routes;
- QW or Hamiltonian kernels;
- backend defaults;
- PGDG or numerical-reference policy;
- supplement handling;
- CR2 or science validation status.

The next decision should be a read-only adoption-risk audit before routing any
real construction through the private shadow helper. Full kinetic adoption
still requires a deliberate production contract for support-dense and mixed
blocks, so it should not be inferred from the private full-shadow result.
Local/Gaussian one-body terms, optimized/adopted PQS/product mixed blocks,
support-dense mixed packet generalization, and all-pairs production assembly
remain separate design questions. Any adoption pass should remain private until
it has explicit operator checks, signed-block diagnostics, and a clear
production-readiness decision.
