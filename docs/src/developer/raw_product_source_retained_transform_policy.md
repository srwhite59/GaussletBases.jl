# Raw Product Source And Retained Transform Policy

For PQS operator work, the active governing framework is
[`pqs_source_box_operator_framework.md`](pqs_source_box_operator_framework.md).
Read that framework before using this note for implementation. If the two
notes appear to disagree, stop and update the framework/policy explicitly
instead of resolving the conflict silently in code.

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

The PQS/product source-box path is already the direct 1D-factor retained-block
path. `_pqs_product_source_box_reference_block(...)` avoids dense 3D raw
source-box pair matrix materialization: it builds 1D cross-axis factors,
selects PQS boundary COMX-product modes, applies product/doside retained mode
metadata, and assembles the retained block directly. The supported terms are
overlap, `position_x/y/z`, `x2_x/y/z`, and kinetic, with kinetic represented
as `(K,S,S) + (S,K,S) + (S,S,K)`. The helper consumes explicit axis
metric/operator data and reports that it did not invoke a numerical fallback;
PGDG/exact provenance is only claimed when upstream raw-box or fixture
diagnostics actually check that provenance.

The private multi-term checkpoint adds
`_pqs_product_source_box_reference_blocks_from_pair_plan(...)` and
`_pqs_product_source_box_reference_blocks(...)`. These helpers reuse one
`_pqs_product_source_box_pair_plan(...)` across requested terms and call the
same direct retained-block assembly for each term. `_pqs_product_source_box_shadow_blocks(...)`
now uses that multi-term path for its PQS/product component blocks. The
PQS/PQS source-box seam is now also explicit:
`_pqs_pqs_source_box_pair_plan(...)`,
`_pqs_pqs_source_box_reference_blocks_from_pair_plan(...)`,
`_pqs_pqs_source_box_reference_blocks(...)`, and
`_pqs_pqs_source_box_reference_block(...)` build mode-selected PQS/PQS blocks
from 1D source-box factors and boundary COMX-product mode selectors. Self and
compatible cross blocks now use the helper-internal explicit raw product-box
boundary-column selection oracle; the older raw-box self helper remains a
reference/debug path when that explicit validation is disabled.
Cross-PQS blocks now accept distinct compatible raw product-box plans when
source-mode dimensions, source-mode ordering, and boundary selector structure
match. The cross helper builds 1D cross factors
`C_left' * M[left_interval, right_interval] * C_right`, then assembles the
retained block directly from boundary selectors. The first focused cross-PQS
fixture uses shifted cubic boxes, left `(1:5,1:5,1:5)` and right
`(3:7,1:5,1:5)`, with source dims `(5,5,5)` and retained count `98`. Its
explicit source-box oracle is now internal to the private helper and agrees to
about `5.7e-14`; reverse-orientation transpose consistency is about
`7.6e-15`. Dense raw product-box pair matrices are materialized only for that
validation oracle, not for the streamed 1D-factor algorithmic block. The
supported terms remain overlap, `position_x/y/z`, `x2_x/y/z`, and kinetic.

`_pqs_product_source_box_shadow_blocks(...)` now carries a tiny private
all-pairs inventory for its two retained units. The units are `:pqs` with
source family `:mode_selected_raw_product_box` and `:product` with source
family `:product_doside`. The upper-triangular entries are `(:pqs, :pqs)` via
`:_pqs_pqs_source_box_reference_blocks`, `(:pqs, :product)` via
`:_pqs_product_source_box_reference_blocks`, and
`(:product, :product)` via
`_product_doside_source_box_reference_block(...)`. That product/product path
uses source-box vocabulary while still comparing to the existing
product-staged retained helpers as authority. The
shadow layout now uses the PQS/PQS helper for its PQS/PQS component while
keeping the existing PQS/product and product/product paths. Diagnostics record
the cost and boundary explicitly:
`dense_raw_source_box_pair_matrix_materialized = false`,
`dense_raw_pair_storage_avoided = true`,
`retained_block_assembled_directly_from_1d_factors = true`,
`source_box_pair_storage_scaling = :one_dimensional_factors_plus_retained_block`,
`pair_plan_reused_for_terms = true` on multi-term paths, and
`all_pairs_inventory_private = true` for the two-unit shadow inventory.
Existing single-term helpers remain compatibility wrappers and no packet,
fixed-block, QW/Hamiltonian, IDA/MWG, retained-weight, public/default, or
science route is adopted. The inventory is not a production all-pairs packet
builder.

The next route-like private checkpoint adds
`_pqs_pqs_product_source_box_all_pairs_inventory(...)` and
`_pqs_pqs_product_source_box_shadow_blocks(...)`. This hard-coded three-unit
shadow has retained units `(:pqs_left, :pqs_right, :product)` and exactly six
upper-triangular pair entries: `(:pqs_left, :pqs_left)`,
`(:pqs_left, :pqs_right)`, `(:pqs_left, :product)`,
`(:pqs_right, :pqs_right)`, `(:pqs_right, :product)`, and
`(:product, :product)`. PQS/PQS entries use
`:_pqs_pqs_source_box_reference_blocks`, PQS/product entries use
`:_pqs_product_source_box_reference_blocks`, and the product/product entry
uses `_product_doside_source_box_reference_block(...)`, which still compares
to the existing product-staged retained helpers as authority. The focused
fixture uses the shifted cubic PQS pair plus a small product/doside slab; the
full retained dimension is `200`, pair count is `6`, and component max error
is about `5.7e-14`. Product/product, PQS/product, and PQS/PQS now participate
in the same private source-box vocabulary for overlap, position, `x2`, and
kinetic. This remains a private shadow/reference inventory, not a generic
route inventory framework and not packet construction adoption.

The route-shaped private consumer checkpoint now has a small descriptor
normalizer. `_pqs_pqs_product_safe_term_route_descriptor(...)` records an
already-built route consisting of a left PQS raw plan, a right PQS raw plan,
and one product/doside unit. It records roles, units, unit summaries, expected
ranges, retained dimension, retained unit count, expected pair count,
supported safe terms, metadata/provenance, and no-adoption diagnostics. This
descriptor is private/shadow-only route metadata: it does not construct route
geometry and is not a generic retained-unit framework.

`_pqs_pqs_product_route_shaped_safe_term_consumer(...)` accepts the descriptor
route kind `:pqs_pqs_product_source_box_safe_term_route` while preserving
compatibility with the older fixture route kind. The consumer still delegates
all numerical work to `_pqs_pqs_product_source_box_shadow_blocks(...)`; no new
operator algebra was added. `_pqs_pqs_product_supported_safe_terms(...)`
centralizes validation for the current safe set: overlap, `position_x/y/z`,
`x2_x/y/z`, and kinetic. Unsupported terms such as `:weights` reject before
being treated as route-supported.

Commit `770b7be` records this route-shaped raw-box consumer as a private
checkpoint. Every route pair is explicitly labeled
`:source_box_algorithm_available`. Cross-PQS/PQS uses the helper-internal
explicit raw product-box boundary-selection oracle for validation. Dense raw
source-box pair matrices remain validation-only, not the algorithmic path.
The route consumer remains source-box-first and avoids shell projection,
Lowdin cleanup, support-local fallback, support coefficient matrices,
retained PQS weight semantics, and IDA division.

The private raw-box route producer checkpoint is recorded by commits
`95d7b11` and `804bdd9`. It starts from explicit fixture facts and produces
the same descriptor through
`RawProductBoxPlan -> RetainedRule -> route descriptor`. The producer uses
left/right mode-selected raw-box PQS retained rules plus an identity
product/doside slab retained rule, then feeds the produced descriptor into
`_pqs_pqs_product_route_shaped_safe_term_consumer(...)`. Sampled validation
covers a shifted cubic `q5/L5` fixture and a rectangular `q5/L7` fixture with
`L != q`; consumer output matches the source-box shadow or hand-built route
path to roundoff. Timing and allocation summaries are captured as diagnostic
evidence only, not as performance thresholds. Dense raw source-box matrices
remain validation-only.

This producer checkpoint remains private/shadow-only. It adds no shell
projection, Lowdin cleanup, support-local PQS oracle, support coefficient
matrix use, retained PQS weights, IDA division, packet or fixed-block
adoption, QW/Hamiltonian routing, public/default behavior, local/ECP/Gaussian/
MWG/interaction terms, IDA/MWG change, or CR2 science claim.

Commit `047af1d` adds the private geometry/recipe facts producer checkpoint
for this same lane. The geometry facts helper accepts `parent_dims`,
`bond_axis`, `q`/`L` or explicit `source_mode_dims`, `left_start`,
`right_shift`, and a product slab fixed index/rule. It emits the explicit
left/right PQS source boxes, product/doside slab source box, source-mode
dimensions, metadata, provenance, and diagnostics consumed by the raw-box
route producer.

This is private fixture infrastructure only. It is not a general diatomic
route geometry policy, public builder, packet-adoption seam, or operator
authority. The shifted cubic `q5/L5` and rectangular `q5/L7` samples match the
explicit-fixture producer and safe-term consumer path to roundoff. The
negative boundary remains unchanged: no shell projection, Lowdin cleanup,
support-local fallback as an algorithm, support coefficient matrices,
retained PQS weights, IDA division, packet or fixed-block adoption,
QW/Hamiltonian routing, public/default behavior, local/ECP/Gaussian/MWG/
interaction terms, IDA/MWG change, or CR2 science claim.

Commit `17dd86d` validates the same private geometry facts helper for `:z` and
`:x` fixture bond axes. The `:x` sample emits the expected left/right PQS
source boxes, product/doside slab source box, product slab fixed-axis
metadata, source-mode dimensions, retained dimension, pair count, pair policy,
and safe-term consumer output matching the explicit route-producer path to
roundoff. This is only mechanical axis-label coverage for explicit fixture
facts; it is not a general atom-centered, shell-realized current-route, or CR2
geometry builder.

The `17dd86d` checkpoint also records focused guard coverage for invalid
`bond_axis`, missing `q` when `source_mode_dims` is absent, malformed
source-mode dimensions, and source-mode axis lengths below two. The intended
flow remains:
`geometry/recipe facts -> explicit source boxes/source-mode dimensions ->
RawProductBoxPlan -> RetainedRule -> route descriptor -> source-box safe-term
consumer`.

Commit `28c3dbc` records the corresponding private input-gate checkpoint. The
gates protect private fixture construction from accidental misuse: source
boxes must be nonempty and inside `parent_dims`; `parent_dims` must be a
positive 3D integer tuple; PQS source-mode dimensions are total source-mode
lengths with at least two modes per axis for boundary selection; the identity
product/doside slab requires exactly one fixed axis; and unsupported safe
terms reject before source-box construction. These gates are not a public
route contract and do not add shell projection, Lowdin cleanup, support-local
PQS oracle behavior, support coefficient matrices, retained PQS weights, IDA
division, packet/fixed-block/QW/Hamiltonian adoption, public/default routing,
local/ECP/Gaussian/MWG/interaction terms, IDA/MWG changes, or CR2 claims.

The private route-fact adapter checkpoint adds
`_pqs_pqs_product_route_descriptor_diagnostic(route_like, metrics = nothing; ...)`.
This is diagnostic/read-path infrastructure only. It emits
`status = :descriptor_available` only when the input already supplies a left
raw-box PQS plan, a right raw-box PQS plan, and an explicit
`_CartesianNestedProductStagedByCenterUnit3D(kind = :product_doside)`. The
hand-built descriptor-available route remains consumable by the existing
private route-shaped safe-term consumer. The current high-order PQS source
construction returns `status = :descriptor_unavailable`; the focused fixture
records `pqs_descriptor_count = 1`, `pqs_raw_plan_convertible_count = 1`,
`product_doside_unit_count = 0`, and
`direct_or_support_body_piece_count = 4`. Missing facts include the second PQS
raw plan and the middle product/doside unit. Direct/support pieces, including
the contact cap, are reported as mismatches and are not reinterpreted as
product/doside. Diagnostics preserve no packet adoption, no fixed-block
construction change, no QW/Hamiltonian change, no shell projection/Lowdin, no
support-local PQS oracle, no retained-weight IDA division, and no
local/ECP/Gaussian/MWG/interaction changes.

The private current-route retained-unit inventory checkpoint adds
`_pqs_current_route_retained_unit_inventory(...)` as diagnostic/read-path
infrastructure only. It does not emit a route descriptor or whole-route
safe-term matrix consumer. For the focused q4 route, it represents all current
columns `1:487` exactly once as six ordered retained units: outer mismatch
low/high product/doside bridge slabs at `1:49` and `50:98`; left/right
atom-box support-dense direct-support units at `99:223` and `224:348`; the
contact-cap product/doside bridge at `349:373`; and the regular shared
molecular shell as a shell-realized PQS fixture at `374:487`. Coverage checks
record ordering by column range, contiguity, non-overlap, first column `1`,
last column `487`, and complete current-route representation.

This inventory keeps the shared PQS boundary honest. The shared shell is
labeled `:shell_realized_pqs_fixture`: support-local/shell-realized
coefficients describe the current shell-supported fixture representation, not
the intended operator algorithm. Raw product-box and source-box metadata remain
auxiliary compatibility metadata. A future source-box operator block must start
from an explicit object contract, not from deriving an algorithmic transform
out of this fixture. Product/product pairs use the product/doside source-box
path. Support/support and support/product pairs use support-local fallback
unless both sides are product/doside. Shell-realized PQS/product, PQS/support,
and PQS/PQS support-local contractions are labeled
`:support_local_oracle_for_shell_realization`: they are shell-row oracle/debug
validation paths, not active algorithmic pair policies.

The shell-realization fact checkpoint records compatibility evidence for each
`:shell_realized_pqs_fixture`. It records the source-box axis intervals and
total source-mode dimensions, boundary COMX-product mode selection, shell-row
projection stage, full-rank symmetric Lowdin cleanup stage,
cleanup/support-local coefficient shapes, and support/retained counts. It
checks that descriptor projection plus Lowdin reproduces the stored
support-local coefficients. It deliberately reports
`compact_source_space_transform.available = false` and
`source_box_operator_application_ready = false`: the fixture is validation and
compatibility metadata, not the object that defines the source-box-first PQS
operator algorithm.

The checkpoint is metadata/diagnostic only: no construction mutation, sidecar
installation, packet or fixed-block adoption, QW/Hamiltonian change, IDA/MWG
semantic change, local/ECP/Gaussian/interaction work, public/default route
change, CR2 status change, retained PQS positive-weight semantics, or retained
PQS IDA division is implied.

The all-pairs route checkpoint expands that q4 retained-unit inventory into 21
upper-triangular retained pairs through
`_pqs_current_route_retained_pair_inventory(...)`. Pair-policy counts are:
product/product `6`, support/support `3`, support/product `6`,
shell-realized PQS/product `3`, shell-realized PQS/support `2`,
shell-realized PQS/PQS `1`, and raw-box PQS active pair policies `0`.
Product/product pairs use the product/doside source-box path. Support-dense
pairs use the support-local fallback path. Shell-realized PQS pairs are
counted separately as support-local oracles for shell realization. In the q4
checkpoint, 15 pairs remain active algorithmic/read-path pairs, 6
product/product pairs have a source-box algorithm available, and 6
shell-realized PQS pairs are oracle-only shell-row contractions.

`_pqs_current_route_safe_term_matrices(...)` consumes that pair inventory and
builds full private diagnostic `487 x 487` matrices for overlap,
`position_x/y/z`, `x2_x/y/z`, and kinetic. Every pair is compared against a
support-local oracle. For shell-realized PQS pairs, that oracle is debug
validation only; it is not the intended route algorithm. Focused q4 evidence
was `542 / 542`, max matrix error `0.0`, helper elapsed time
`28.601491167 s`, and allocation `137,976,768` bytes.

The private current-route authority checkpoint adds
`_pqs_current_route_safe_term_authority_comparison(...)`. It compares the
private route-wide safe-term matrices against existing current-route authority
fields only. Candidate authorities are the fixed block first and the sequence
packet second. A term is compared only when the authority object carries a
finite retained-space matrix with shape `(retained_dimension,
retained_dimension)`; missing or wrong-shaped terms are recorded explicitly.
For this q4 checkpoint, overlap and kinetic are required authorities.

Focused q4 evidence: the helper compared all safe terms against fixed-block
authority: overlap, `position_x/y/z`, `x2_x/y/z`, and kinetic. The focused
probe passed `591 / 591`; the support-local oracle max error remained `0.0`;
the fixed-block authority max error was `2.6290081223123707e-13`; authority
comparison time/allocation were `0.007607708 s` and `15,227,616` bytes. The
underlying safe-term matrix build in that run took `28.357196959 s` and
allocated `137,977,696` bytes.

The active shared PQS shell remains `:shell_realized_pqs_fixture`; raw-box PQS
helpers remain reference/shadow-only for the current route. This is still
private diagnostic/shadow infrastructure only: no construction behavior
change, sidecar installation, packet/fixed-block adoption,
`_nested_shell_packet(...)` authority change, QW/Hamiltonian behavior change,
IDA/MWG behavior change, retained PQS weight division or positive
quadrature-weight claim, local/ECP/Gaussian/MWG/interaction implementation,
public/default route change, or Be2/Cr2 science claim is implied.

The focused homonuclear-style fixture uses parent/bundle shape `(5,5,7)`, left
PQS box `(1:5,1:5,1:5)`, right PQS box `(1:5,1:5,3:7)`, and a middle product
slab at `z = 4`. Both PQS source boxes use source-mode dimensions `(5,5,5)`
with retained count `98`; the product slab retained count is `25`; the total
retained dimension is `221`.

The ignored probe `tmp/work/validate_route_shaped_safe_term_consumer.jl` now
runs a focused scaling ladder and writes an ignored TSV under `tmp/work/`.
Each route has three retained units, six upper-triangular pairs, eight safe
terms, `dense_raw_source_box_pair_matrix_materialized=false`,
`dense_raw_pair_storage_avoided=true`, unsupported `:weights` rejected, and
max full/component error `0.0` against the direct shadow helper.

| route | retained dim | pairs | terms | elapsed | allocated bytes | max error |
|---|---:|---:|---:|---:|---:|---:|
| `q5_L5_slab5` | 221 | 6 | 8 | about `3.13 s` | about `60 MB` | `0.0` |
| `q5_L7_slab5` | 285 | 6 | 8 | about `0.002 s` | about `15 MB` | `0.0` |
| `q5_L9_slab5` | 349 | 6 | 8 | about `0.003 s` | about `22 MB` | `0.0` |
| `q7_L7_slab7` | 485 | 6 | 8 | about `0.005 s` | about `42 MB` | `0.0` |

The first q5/L5 timing includes compilation/warmup; the later rows are the
more useful small-fixture steady-state signal.

This checkpoint is still private shadow/reference infrastructure. It does not
adopt packet or fixed-block construction, QW/Hamiltonian assembly,
public/default behavior, CR2 status, local/ECP/Gaussian/MWG/interaction terms,
shell projection, Lowdin in raw-box operators, support-local PQS oracles,
retained PQS positive-weight semantics, IDA division, or dense raw source-box
pair matrix storage. It is not yet a Be2/Cr2 route benchmark.

A local ignored diagnostic probe records the private performance shape for a
rectangular `(5,5,7)` PQS source box and a non-identity product retained
transform. In one run with 100 repetitions, repeated single-term calls took
about `0.0067 s` and allocated about `32.3 MB`, while the multi-term pair-plan
reuse path took about `0.0014 s` and allocated about `7.9 MB`. This is local
diagnostic evidence only, not a public benchmark or readiness claim.

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

The private retained-unit fact audit adds one more diagnostic vocabulary layer
without changing construction:
`_pqs_route_retained_unit_fact_audit(construction; include_support_indices=false)`
classifies existing high-order PQS route regions before any route descriptor
is emitted. In the focused q4-style construction, contact cap is
`:product_box_constructible` with an identity-selector `q x q x 1` slab rule;
outer mismatch is `:product_box_constructible` with boundary slab construction
rules; left/right atom boxes are `:needs_direct_support_retained_unit_kind`;
and the regular shared molecular shell is `:out_of_scope` for body vocabulary
because it is still the current single PQS descriptor.

This audit keeps the contract wording strict.
`raw_product_box_operator_contract = true` means the active retained unit
already carries that raw product-box operator contract, currently expected for
already-built product/doside units.
Current direct/support pieces that are merely product-box constructible keep
`raw_product_box_operator_contract = false` and instead report
`product_box_construction_rule_available = true` when an explicit rule is
present. Product-box-constructible direct/support pieces are not
reinterpreted as product/doside, no unit is created, and the current high-order
PQS route descriptor diagnostic remains `:descriptor_unavailable`.

The contact-cap product/doside checkpoint adds one private read-path helper:
`_pqs_contact_cap_product_doside_unit(construction; ...)`. It constructs an
explicit contact-cap-only
`_CartesianNestedProductStagedByCenterUnit3D(kind = :product_doside)` from the
audited `q x q x 1` identity-selector slab rule. The helper proves
equivalence to the current direct/support contact-cap selector by matching
support indices and states, retained count, column range, identity
support-local coefficients, and parent-expanded coefficient error `0.0`.
Its diagnostics distinguish the old input fact from the helper-created unit:
`input_fact_raw_product_box_operator_contract = false` and
`created_unit_raw_product_box_operator_contract = true`.

This contact-cap helper is private diagnostic/read-path infrastructure only.
It does not install the unit into construction or sidecars, does not emit a
route descriptor, and does not convert outer mismatch, atom boxes, or PQS
descriptors. Current route descriptor diagnostics may therefore remain
`:descriptor_unavailable`. It changes no packet or fixed-block construction,
QW/Hamiltonian path, public/default behavior, retained-weight IDA division, or
local/ECP/Gaussian/MWG/interaction behavior.

The contact-cap safe-term operator checkpoint adds
`_pqs_contact_cap_safe_term_operator_comparison(construction, metrics; ...)`,
also as private diagnostic/read-path infrastructure only. It compares
contact-cap product/doside self blocks from the helper-created unit against
the current direct/support contact-cap selector oracle. The product path uses
`_product_doside_source_box_reference_block(...)`; the direct/support oracle
uses the current contact-cap build, support indices/states, and support-local
direct selector entries.

The covered safe terms are overlap, `position_x/y/z`, `x2_x/y/z`, and kinetic.
Kinetic is the signed separable sum `(K,S,S) + (S,K,S) + (S,S,K)`. The focused
q4 fixture checks all eight terms as `25 x 25` blocks, finite outputs, and max
product-vs-direct/support error `<= 1.0e-12`; unsupported terms such as
`:weights` reject. The helper consumes caller-supplied explicit axis
metric/operator data and does not claim new PGDG analytic provenance inside
the helper. It emits no route descriptor, mutates no construction, installs no
sidecar, and changes no packet/fixed-block/QW/Hamiltonian/IDA/MWG/local/
Gaussian/public behavior.

The outer-mismatch bridge now has the analogous private diagnostic checkpoint.
`_pqs_outer_mismatch_product_doside_units(construction; ...)` maps the current
boundary slab-set direct/support piece to one product/doside unit per slab
piece, not one unit for the whole slab set. In the q4 fixture the two z slabs
have column ranges `1:49` and `50:98`; each uses identity support-local
coefficients and the structural parent-expanded selector error is `0.0`.
Descriptor-piece order defines the coefficient columns, while the audited
region support is checked as an order-independent coverage set.

`_pqs_outer_mismatch_safe_term_operator_comparison(construction, metrics; ...)`
then assembles the complete `98 x 98` retained outer-mismatch block from all
four slab-pair product/doside blocks: low/low, low/high, high/low, and
high/high. The comparison target is the current direct/support outer-mismatch
selector oracle. The covered safe terms are the same eight terms as contact
cap: overlap, `position_x/y/z`, `x2_x/y/z`, and kinetic. The q4 focused
fixture checks finite product/direct outputs and max product-vs-direct/support
error `<= 1.0e-12`; unsupported terms such as `:weights` reject.

Both contact-cap and outer-mismatch bridges are private diagnostic bridges from
current direct/support route pieces into product-box retained units. They are
not installed into construction, do not emit route descriptors, and do not
change sidecars, packets, fixed blocks, QW/Hamiltonian, IDA/MWG,
local/ECP/Gaussian/interaction, public/default routes, or CR2 artifacts. Their
operator inputs are caller-supplied explicit axis data; the helpers do not
independently prove PGDG analytic provenance.

Atom boxes are now recorded through a different private diagnostic bridge:
support-dense/direct-support retained units, not product/doside units.
`_pqs_atom_box_support_dense_units(construction; ...)` emits exactly two q4
units, `:left_atom_box` and `:right_atom_box`. Both are `:support_dense`,
with column ranges `99:223` and `224:348`, retained/support count `125` each,
and local support coefficient shape `(125, 125)`. Parent-expanded coefficient
error against the current direct/support atom-box builds is `0.0`. The local
identity/direct-row error is recorded near roundoff, but it is not treated as
a product-box construction rule or raw product-box operator claim.

`_pqs_atom_box_safe_term_operator_comparison(construction, metrics; ...)`
validates those support-dense atom-box units as operator inputs using
support-local fallback. It does not use product/doside algebra or raw
product-box factorization for atom boxes. The q4 comparison assembles the full
`250 x 250` retained atom-box block from left/left, left/right, right/left,
and right/right pair blocks, including cross-atom blocks. The covered safe
terms are overlap, `position_x/y/z`, `x2_x/y/z`, and kinetic. The focused q4
fixture checks max error against the current direct/support oracle
`<= 1.0e-12`; unsupported terms such as `:weights` reject.

Atom boxes therefore remain support-local/direct-support retained units in the
private route vocabulary. No product-box construction rule is claimed for
them, and no route descriptor emission, construction mutation, sidecar
installation, packet/fixed-block/QW/Hamiltonian/IDA/MWG/local/ECP/Gaussian/
interaction, public/default route, or CR2 artifact change is implied. As with
the product-box bridge helpers, metric/operator inputs are caller-supplied
explicit axis data, not a new independent PGDG analytic provenance proof.

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
transpose for symmetric real terms, and product/product blocks. Its PQS/PQS
component uses `_pqs_pqs_source_box_reference_blocks(...)` with helper-internal
explicit boundary-column validation for supported compatible raw-box fixtures;
its PQS/product component uses the multi-term PQS/product pair-plan reuse path;
and its product/product component uses
`_product_doside_source_box_reference_block(...)`, which still compares to the
existing product-staged retained helpers as authority. The
helper records a tiny private all-pairs inventory over units `:pqs` and
`:product` with upper-triangular entries `(:pqs, :pqs)`,
`(:pqs, :product)`, and `(:product, :product)`. The supported terms are
`:overlap`, `:position_x/y/z`, `:x2_x/y/z`, and `:kinetic`. Focused tests
cover a rectangular PQS source box and a non-identity product/doside
transform.

The cross-PQS checkpoint extends the PQS/PQS source-box seam from self-only to
distinct compatible raw product-box plans. Compatibility currently requires
matching source-mode dimensions, source-mode ordering, and boundary selector
structure. Cross-axis factors are built as
`C_left' * M[left_interval, right_interval] * C_right`, and retained blocks are
assembled from boundary COMX-product selectors. The focused fixture shifts the
right cubic PQS box: left `(1:5,1:5,1:5)`, right `(3:7,1:5,1:5)`, source dims
`(5,5,5)`, retained count `98`. The explicit source-box oracle agrees to about
`5.7e-14` inside the helper, and reverse-orientation transpose consistency is
about `7.6e-15`. Dense raw product-box pair matrices are validation-only; the
retained block path still streams 1D factors and boundary selectors.

The three-unit route-like shadow checkpoint then adds
`_pqs_pqs_product_source_box_all_pairs_inventory(...)` and
`_pqs_pqs_product_source_box_shadow_blocks(...)`. It has retained units
`(:pqs_left, :pqs_right, :product)` and six upper-triangular pair entries:
`(:pqs_left, :pqs_left)`, `(:pqs_left, :pqs_right)`,
`(:pqs_left, :product)`, `(:pqs_right, :pqs_right)`,
`(:pqs_right, :product)`, and `(:product, :product)`. The shifted fixture has
full retained dimension `200`, pair count `6`, and component max error about
`5.7e-14`. Product/product, PQS/product, and PQS/PQS now share the private
source-box vocabulary for overlap, position, `x2`, and kinetic; the
product/product leg is source-box-labeled while retaining the existing
product-staged helper comparison as authority.

These checkpoints are still private/shadow-only layout evidence, not a packet
builder and not a generic route inventory framework: they do not use shell-row
projection, Lowdin, `support_coefficient_matrix` as a PQS oracle, retained PQS
weights, or IDA division, and they change no packet construction,
QW/Hamiltonian, public/default, CR2, local/ECP/Gaussian/MWG/interaction, or
IDA/MWG behavior.

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
