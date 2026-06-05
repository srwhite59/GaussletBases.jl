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
- [`raw_product_source_retained_transform_policy.md`](raw_product_source_retained_transform_policy.md)
- [`projected_q_shell_policy.md`](projected_q_shell_policy.md)

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

LW often collapses the intermediate and final stages because its CPB boundary
strata already live on disjoint shell support. PQS must not collapse them unless
a reviewed contract explicitly permits it.

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
-> unit-pair inventory
-> pair operator blocks
-> assembly
```

Low-order White--Lindsey:

```text
ShellificationRegion Ω
-> boundary-stratum CPBs: facets, edges, corners
-> product/doside contractions or direct selectors
-> final retained units
-> unit-pair blocks
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
-> unit-pair blocks
-> assembly
```
