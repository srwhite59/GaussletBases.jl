# Cartesian Coordinate Product Box Contract

This note defines the **Coordinate Product Box** vocabulary for Cartesian
nesting work. It is a developer contract, not a public API. Its purpose is to
keep the current unification work from confusing two different jobs:

1. shellification / atom-growth, which partitions parent support into disjoint
   owned regions; and
2. lowering, which chooses coordinate-product source domains and retained
   transforms for a construction recipe.

The distinction is especially important because the low-order White--Lindsey
route and the projected q-shell (PQS) route use the same shellification regions
but very different lowering geometry.

Read this note together with:

- [`pqs_source_box_operator_framework.md`](pqs_source_box_operator_framework.md)
- [`projected_q_shell_policy.md`](projected_q_shell_policy.md) for historical inbound links

If these notes disagree, stop and update the documentation explicitly before
continuing implementation. Do not resolve the conflict silently in code.

## Core definition

A **Coordinate Product Box** (**CPB**) is an axis-aligned product of coordinate
intervals, with singleton intervals allowed:

```text
CPB = I_x x I_y x I_z
```

The word "box" here means coordinate-product structure, not necessarily a
three-dimensional filled volume. A CPB can have codimension zero, one, two, or
three depending on how many intervals are singletons.

Examples:

```text
filled rectangular box:  I_x x I_y x I_z
facet / face:           {i_x} x I_y x I_z
edge:                   I_x x {i_y} x {i_z}
corner:                 {i_x} x {i_y} x {i_z}
midpoint slab:          I_x x I_y x {i_z}
```

A shell such as

```text
B_outer \ B_inner
```

is **not** a CPB. It is an owned support region or a shell support described by
relations between CPBs. This exclusion is intentional: PQS often uses a filled
source CPB even though the final owned support is a shell.

## CPB roles and codimension

Use codimension for geometry and role for policy. Codimension alone is not
enough, because a codimension-one CPB can be a boundary facet, a midpoint slab,
a contact cap, or an outer mismatch slab.

Suggested vocabulary:

```text
VolumeCPB / FilledCPB
    Codimension 0. A filled rectangular product box.

FacetCPB / FaceCPB
    Codimension 1 boundary stratum of a shell. "Facet" is preferred in
    conceptual docs; "face" remains fine in code where existing names use it.

SlabCPB
    Codimension 1 CPB that is not necessarily a boundary facet. Midpoint,
    contact, and mismatch slabs are slab CPBs.

EdgeCPB
    Codimension 2 boundary stratum.

CornerCPB
    Codimension 3 boundary stratum.

SourceCPB
    A CPB used as a source domain for constructing retained functions or
    operator blocks.

SupportCPB
    A CPB used as an owned or diagnostic support piece.
```

A CPB record should be allowed to carry both geometry and role:

```text
intervals = (I_x, I_y, I_z)
codimension
role = :boundary_facet | :midpoint_slab | :contact_cap | ...
```

## Shellification regions are not CPBs

A `ShellificationRegion` owns a disjoint part of the parent lattice. It is the
output of atom-growth or another shellification policy. Its owned support may be
one CPB, a union of CPBs, or a shell-like set such as `B_outer \ B_inner`.

A shellification region should distinguish:

```text
owned_support Ω
    The disjoint parent rows owned by this region.

source_cpb or source_cpbs
    Coordinate-product source domains used by a lowering recipe.

lowering_recipe
    The recipe-specific rule that turns the source domain into retained
    functions and operator blocks.
```

This separation lets atom-local, diatomic, midpoint, contact, shared molecular,
and outer mismatch regions use one shellification language without forcing the
same lowering algorithm.

## Low-order White--Lindsey lowering

The low-order White--Lindsey route lowers a shellification region by decomposing
its owned shell support into disjoint boundary strata.

For an outer CPB

```text
B = I_x x I_y x I_z
```

and an inner exclusion CPB

```text
B^- = I_x^- x I_y^- x I_z^-
```

the shell support is

```text
Ω = B \ B^-
```

The low-order route decomposes `Ω` into disjoint CPB boundary strata:

```text
facets:  F_a^σ
edges:   E_a^{σ_b σ_c}
corners: C^{σ_x σ_y σ_z}
```

where `a` is an axis and `σ` is a low/high boundary side.

A conceptual stratum definition is:

```text
Facet / face:
F_x^lo = {first(I_x)} x I_y^- x I_z^-
F_x^hi = {last(I_x)}  x I_y^- x I_z^-

Edge:
E_x^{lo,hi} = I_x^- x {first(I_y)} x {last(I_z)}

Corner:
C^{lo,hi,lo} = {first(I_x)} x {last(I_y)} x {first(I_z)}
```

The lowering is then:

```text
facet CPB  -> fixed axis x doside x doside
edge CPB   -> doside x fixed axis x fixed axis
corner CPB -> direct parent site
```

Thus low-order LW is a **boundary-stratum CPB lowering**. It works on small,
disjoint CPBs and naturally produces face/facet, edge, and corner retained
pieces.

### White-Lindsey complete-shell CPB stratum enumeration

Before code enumerates `:facet_cpb`, `:edge_cpb`, and `:corner_cpb` records for
a complete shell, the complete-shell geometry must satisfy the normal
single-layer boundary contract:

```text
B_outer = I_x x I_y x I_z
B_inner = I_x^- x I_y^- x I_z^-
Ω = B_outer \ B_inner
```

For each axis `a`, `I_a^-` is the strict interior interval obtained by removing
exactly the low and high outer boundary points from `I_a`. The normal complete
shell therefore has one low and one high boundary point on every axis, and a
nonempty inner interval on every axis. If a shell does not meet this condition,
do not silently force it into the complete-shell enumeration described here.

Let

```text
L_a = first(I_a)
H_a = last(I_a)
J_a = I_a^-
P_a^low  = {L_a}
P_a^high = {H_a}
```

The complete-shell enumeration is the disjoint CPB decomposition whose union is
`Ω`:

```text
Facet on axis a, side σ:
P_a^σ x J_b x J_c

Edge on axes a,b, sides σ,τ, remaining axis c:
P_a^σ x P_b^τ x J_c

Corner on sides σ_x,σ_y,σ_z:
P_x^σ_x x P_y^σ_y x P_z^σ_z
```

Here `(a,b,c)` ranges over the Cartesian axes. Facets use the outer boundary
point on one axis and the inner intervals on the other two axes. Edges use the
outer boundary points on two axes and the inner interval on the remaining axis.
Corners use outer boundary points on all three axes.

Using inner intervals on the non-boundary axes is what makes the strata
disjoint: a parent point belongs to a facet, edge, or corner according to how
many axes sit on an outer boundary point. The union is `Ω` because every point
of `B_outer \ B_inner` has at least one axis on an outer boundary point under
the complete-shell condition.

Use explicit low/high side labels in canonical axis order. Role names should be
predictable symbols such as:

```text
:x_low_facet
:y_high_z_low_edge
:x_low_y_high_z_high_corner
```

For edges and corners, list only the axes that sit on outer boundary points,
in `x`, `y`, `z` order. This avoids role names that depend on construction
order.

Partial shells, outer mismatch slabs, contact caps, no-endcap shells, and
extra-endcap shells are not complete shells for this enumeration. Represent
them first as explicit slab or piece support decompositions, then lower those
pieces with their own roles. Do not hide such cases behind
`:facet_cpb`/`:edge_cpb`/`:corner_cpb` records that imply a complete-shell
union.

This enumeration is lowering geometry only. Shellification still owns the
disjoint shell support `Ω`. The boundary-stratum CPBs are lowering sources or
support pieces for the White-Lindsey construction, and product/doside
coefficient maps are later construction data. They are not shellification
authority and should not be recorded as if shellification had chosen
product/doside transforms.

## PQS lowering

PQS uses the same shellification regions but a different lowering source.

For a shellification region with owned support

```text
Ω = B_outer \ B_inner
```

PQS first chooses a filled source CPB:

```text
B_source = I_x x I_y x I_z
```

This source CPB is not the shell. It is a filled coordinate-product box used to
build an intermediate contracted mode space.

The PQS lowering sequence is:

```text
filled SourceCPB
-> 1D COMX/source transforms on each axis
-> product-box modes
-> boundary COMX-product mode selection
-> intermediate retained source space
-> source-space operator blocks from 1D factors
-> optional shell projection to owned_support Ω
-> Lowdin cleanup
-> final shell-realized retained unit
```

The source-space operators must be built before shell realization. Shell-row
projection and Lowdin cleanup are a final realization step, not the raw-box
operator rule.

This is why a shell can be called a shellification region while the PQS lowering
immediately goes to a filled CPB: the filled CPB is the source domain, not the
owned shell support.

## Spaces that must remain distinct

Future code and docs should keep these spaces separate:

```text
Parent space
    The full parent lattice. Do not use dense parent-space operator construction
    as the PQS algorithm.

Owned support
    Disjoint rows assigned by shellification. May be a shell, a CPB, or a union
    of CPBs.

Source CPB
    Coordinate-product source domain used by a lowering recipe. PQS source CPBs
    may overlap or nest even though owned supports are disjoint.

Intermediate retained space
    Contracted source-domain space where source operators are first built.
    PQS boundary product modes live here.

Shell realization
    Optional map from intermediate space to final owned shell support. For PQS,
    this is shell projection plus Lowdin cleanup. For many LW pieces, it is
    direct or trivial.

Final retained unit
    Column-owning unit consumed by pair planning and Hamiltonian assembly.
```

The common operator pattern should be:

```text
source CPB_i, source CPB_j
-> 1D cross-axis factors
-> source-retained operator block
-> left/right realization or retained transforms
-> final retained operator block
```

The current private mixed one-body consumer lives at the pair-block
materialization layer, after pair planning has produced a
`PairBlockMaterializationPlan`. It consumes one safe one-body term plus
caller-supplied factor/provider facts and dispatches to existing local
direct/direct, PQS/PQS raw source-space, and White-Lindsey boundary-stratum
adapter selectors. Its output is a local batch result plus compact summary.
It is not route-driver wiring, global operator assembly, Hamiltonian assembly,
Coulomb, IDA/MWG, artifact/export, PQS shell projection/Lowdin, or full
White-Lindsey route assembly.

LW often collapses the intermediate and final stages because its CPB boundary
strata already live on disjoint shell support. PQS must not collapse them unless
a reviewed contract explicitly permits it.

## Layer responsibility addendum

Each layer should have one main output type.

```text
shellification -> ShellificationRegion / owned shell support
lowering       -> CPBs plus a lowering recipe
construction   -> intermediate and final retained spaces
pair planning  -> pairs of final retained units
```

In this vocabulary, **atom-growth** is a shellification policy, not the name of
the whole low-order or PQS route. It answers "who owns which parent rows?" and
produces shells or other owned support regions.

Lowering is the first recipe-specific step. It answers "which CPBs are used to
construct retained functions from this owned support?"

For low-order White--Lindsey, lowering breaks a shell into disjoint
boundary-stratum CPBs:

```text
ShellificationRegion Ω -> facet CPBs + edge CPBs + corner CPBs + direct CPBs
```

For PQS, lowering chooses a filled source CPB associated with the shell:

```text
ShellificationRegion Ω -> filled source CPB B_source
```

This is sometimes informally described as making the box out of the shell, but
the shell itself is still not a CPB. The filled source CPB is a coordinate
product source domain; selected source modes are later realized back onto the
owned shell support by projection and cleanup.

Thus the layers should not leak into each other: shellification should not own
face/edge/corner or COMX-mode decisions, and pair planning should not start from
shells or source CPBs directly. Pair planning starts after construction, from
final retained units that keep links back to their owned support, source CPBs,
intermediate retained space, shell realization, and lowering recipe.

## Suggested object vocabulary

Use this vocabulary in future developer notes and private contracts:

```text
ShellificationRegion
    Disjoint owned support and route role.

CoordinateProductBox / CPB
    Axis-aligned product of intervals, singleton intervals allowed.

LoweringSource
    A CPB plus recipe-specific interpretation.

IntermediateRetainedSpace
    Contracted source-space object where operators can be built.

ShellRealization
    Optional projection/cleanup map from intermediate space to final support.

FinalRetainedUnit
    Column-owning retained object used by pair planning.

UnitPair
    Upper-triangular pair of final retained units.

PairOperatorBlock
    Operator block between units. It should record whether it came from a
    source-CPB algorithm, a materialized adapter, or an oracle.
```

For LW, lowering sources are commonly:

```text
LWFacetCPB
LWEdgeCPB
LWCornerCPB
LWDirectSlabCPB
```

For PQS, lowering sources are commonly:

```text
PQSFilledSourceCPB
PQSBoundaryModeRetainedSpace
PQSShellRealization
```

## Implementation guidance for Codex

Future implementation work should follow these rules:

1. Do not call `B_outer \ B_inner` a CPB. It is owned support or shell support.
2. Do not treat LW face/edge/corner strata as the universal model for pair
   blocks. They are one lowering recipe.
3. Do not treat PQS shell-row projection plus Lowdin as the raw product-box
   operator construction. It is shell realization.
4. Do not introduce dense parent-space operator construction as the PQS
   algorithm. PQS source operators should be built from 1D product factors on
   source CPBs.
5. Do distinguish `owned_support`, `source_cpb`, `intermediate_retained_space`,
   `shell_realization`, and `final_retained_unit` in new contracts.
6. When validating current LW materialization, label packet slicing or
   fixed-block slicing as an adapter or oracle unless it is the intended
   algorithmic pair-block path.
7. When adding PQS to the driver, introduce it as a lowering recipe for a
   shellification region through a filled source CPB, not as a parallel broad
   route that bypasses the shared geometry/source/realization vocabulary.

## Minimal flow diagrams

Common early route:

```text
parent lattice
-> shellification regions with disjoint owned support Ω
-> source CPBs and lowering recipes
-> intermediate retained spaces
-> shell realization, if needed
-> final retained units
-> retained-unit inventory
-> retained-unit transform contracts
-> unit pairs
-> pair operator plans
-> pair-block materialization preflight/direct-direct and PQS source pilots
-> final pair-block assembly
-> assembly
```

Low-order White--Lindsey:

```text
ShellificationRegion Ω
-> boundary-stratum CPBs: facets, edges, corners
-> product/doside contractions or direct selectors
-> final retained units
-> retained-unit transform contracts
-> unit pairs
-> pair operator plans
-> pair-block materialization direct-direct and PQS source pilots
-> LW boundary-stratum adapter preflight
-> local LW unit coefficient maps and one-body adapter pilots
-> final pair-block assembly
-> assembly
```

PQS:

```text
ShellificationRegion Ω
-> filled source CPB B
-> product-box modes
-> boundary COMX-product retained rule
-> intermediate source-retained operators
-> shell projection + Lowdin realization
-> final retained units
-> retained-unit transform contracts
-> unit pairs
-> pair operator plans
-> pair-block materialization preflight/direct-direct and PQS source pilots
-> final pair-block assembly
-> assembly
```

The layering rule remains:

```text
shellification owns support
lowering chooses CPBs and recipes
construction plans retained units and realization
pair planning starts from final retained units and transform contracts
```

Pair planning must use final retained units plus retained-unit transform
contracts. It must not rediscover realization paths by inspecting retained-unit
kinds. Current pair-block materialization is a preflight layer plus
direct/direct final local one-body pilots, PQS/PQS raw source-space safe-term
pilots, and local White--Lindsey boundary-stratum one-body adapter pilots. The
PQS source pilot does not apply shell projection, Lowdin, or final
retained-block assembly. Its current bridge summaries are metadata-only
records of the later shell-realization handoff, not the realization itself.
The final PQS pair-block readiness summary consumes those single or batch
bridge summaries and currently blocks on `:shell_realization_not_materialized`;
that blocker is now scoped to the bridge-level source-space summary. The
explicit PQS one-electron final-basis seam exists separately for retained
boundary overlap/kinetic and separated by-center nuclear matrices. The bridge
readiness summary still does not build shell projection, Lowdin, final retained
blocks, Hamiltonians, exports, artifacts, IDA/MWG data, or Coulomb. Broader
route-driver adoption, Coulomb/IDA, and many-electron assembly remain future
work.

For low-order White--Lindsey boundary-stratum retained-unit pairs, pair-block
materialization now recognizes `:white_lindsey_boundary_stratum_adapter_path`
as `:white_lindsey_boundary_stratum_adapter_preflight`. This remains the
adapter-boundary checkpoint. Behind that boundary, local old-kernel-backed
unit coefficient maps now exist for facet/face, edge, and corner strata, local
pair-level coefficient gathering exists, and local one-body adapter blocks now
exist for overlap, position_x/y/z, x2_x/y/z, and kinetic. These are local
adapter pilot blocks, not full route/operator assembly; they do not build
Hamiltonians, exports, artifacts, IDA/MWG data, or Coulomb.

`white_lindsey_boundary_stratum_adapter_summary(record)` is the internal
metadata helper that records the old-kernel reuse map for that adapter
boundary. It consumes
`:white_lindsey_boundary_stratum_adapter_preflight` records. Facet/face strata
point to `_nested_face_product`, edge strata to `_nested_edge_product`, corner
strata to `_nested_corner_piece`, and facet/edge side helpers to
`_nested_doside_1d`. For batch or plan-level inputs, the helper reports
`reuse_metadata_available_count` and `reuse_metadata_blocked_count`; these
counts are about old-kernel reuse metadata availability, not full route
assembly readiness.

`white_lindsey_boundary_stratum_unit_adapter_descriptor(unit)` is the compact
unit-level adapter descriptor. It records source-CPB and kernel-input facts
only: unit identity, stratum kind, source CPB role/codimension/count,
active/fixed axis metadata, planned old kernel symbol, and the planned
`_nested_doside_1d` helper for facet/edge strata. The materialized unit
coefficient helper is `white_lindsey_boundary_stratum_unit_coefficients(...)`;
it builds local adapter input maps, not route-global state.

`white_lindsey_boundary_stratum_pair_adapter_descriptor(record[, unit_pair])`
is the compact pair-level adapter descriptor. With unit-pair context it uses
the unit descriptors; record-only use reports only record-derived facts. It
classifies upper-triangular LW boundary-stratum pairs as facet/facet,
facet/edge, facet/corner, edge/edge, edge/corner, or corner/corner. It does
not build one-body adapter blocks. Pair-level coefficient gathering is handled
by `white_lindsey_boundary_stratum_pair_unit_coefficients(...)`.

`white_lindsey_materialized_seed_oracle_summary(...)` is a compact validation
oracle over the old materialized seed. It reports counts, ranges, roles,
packet/operator inventory availability, one-body operator availability flags,
and fixed-block matrix dimension summaries. It is validation-only: not route
authority and not adapter authority. It must not be used to make the old route
the new construction path.

`white_lindsey_boundary_stratum_one_body_block(...)` and
`white_lindsey_boundary_stratum_one_body_blocks(...)` are the local
one-body adapter pilots for overlap, position_x/y/z, x2_x/y/z, and kinetic.
`white_lindsey_boundary_stratum_one_body_adapter_summary(...)` reports compact
supported-term/stratum readiness and batch materialized/skipped summaries.
These helpers reuse old kernels as adapter inputs, not as route authority, and
do not build Coulomb, IDA/MWG data, Hamiltonians, exports, artifacts, or
production dense-parent fallback.

Focused old-seed oracle validation currently covers selected representative
local boundary-stratum pairs. The selected facet is the old seed x-low `yz`
face (face index 5, retained range `162:170`), the selected edge is the old
seed x-high/y-low `z` edge (edge index 11, retained range `210:212`), and the
selected adjacent corner is the old seed x-high/y-low/z-high corner (corner
index 6, retained range `221:221`). Representative facet/edge, facet/facet,
edge/edge, edge/corner, and corner/corner adapter blocks for overlap,
position_x/y/z, x2_x/y/z, and kinetic compare against those old fixed-block
slices within the focused test tolerance. Corner unit coefficients now expose
parent support indices when `parent_dims` and fixed coordinates are available;
if `parent_dims` is absent, the support-local metadata path remains explicit
and no parent support index is guessed. This is a local adapter-slice
checkpoint only; it is not exhaustive over all faces/edges/corners and does
not assemble a full White--Lindsey route.

The focused oracle-comparison test includes an aggregate coverage summary for
this selected checkpoint: 5 pair families x 8 one-body terms = 40 value
comparisons, with zero blocked comparisons in that focused coverage. This test
is a focused validation gate for the LW oracle comparison surface, not a
casual tiny test to run on unrelated passes.

### Next low-order White--Lindsey adapter target

The first local White--Lindsey one-body adapter surface now exists behind
`:white_lindsey_boundary_stratum_adapter_preflight`. It covers local unit
coefficient maps, local pair-level coefficient gathering, local one-body blocks
for overlap/position/x2/kinetic, and compact summaries.

Next useful targets are additional representative face/edge orientations if
needed, a compact validation summary helper, and structured plan consumption
only after review. Keep oracle tests split from routine runners if validation
runtime grows toward the ceiling. Coulomb, IDA, Hamiltonian/export work, and
route assembly remain later work requiring explicit design review. Do not make
support-row or dense-parent fallback the production algorithm, and do not treat
the old White--Lindsey route as the new route authority.
