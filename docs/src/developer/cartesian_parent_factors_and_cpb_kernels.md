# Cartesian Parent Factors and CPB Block Providers

Date: 2026-06-09

This note records a proposed architectural direction for Cartesian parent
operator data, coordinate-product-box local blocks, and the staged PQS /
White-Lindsey route. It is a design note, not an implementation contract yet.

The immediate motivation is the private global-overlap input-facts work. Recent
tests show two different real dry-run report states:

- a manual-count dry report can carry parent axis counts, but has no parent
  axis-bundle factors and correctly blocks on
  `:missing_parent_axis_bundle_overlap_factors`;
- a probe-enabled dry report shows that the structured axis-bundle object
  already exists under `route_materializer_payload.parent_axis_bundle_object`.
  The private overlap facts helper can read those 1D overlap factors through a
  narrow structured source path, but the route materializer payload should not
  become the long-term owner of universal parent-axis data.

That is useful evidence, but it should not turn route-driver report payloads
into the long-term owner of parent operator data. The better direction is to
promote the underlying idea into a parent-owned axis factor packet and a CPB
block-provider layer.

## Core Correction

`CartesianParentGaussletBasis3D` should remain the route-neutral identity of the
full Cartesian parent lattice:

- the three 1D parent axes;
- the full parent lattice box;
- Cartesian flattening and center lookup;
- axis sharing;
- parent metadata.

It should not directly own all numerical operator arrays. Eagerly attaching
large matrices or Gaussian factor tensors to every parent basis object would
make the identity object too heavy and would blur construction policy with
lattice identity.

The parent layer, however, should own companion objects derived from that
parent. The important object is not a route report field and not a PQS object.
It is universal parent-axis factor data.

The guiding rule:

```text
If a factor depends only on:
  parent axes,
  coordinate mappings or distortions,
  Gaussian expansion,
  nuclei, charges, and locations,
and does not depend on:
  nesting choices,
  source boxes,
  retained units,
  shell or PQS projection,
  Lowdin cleanup,
  pair plans,
  or global placement,
then it belongs to the Cartesian parent factor layer.
```

Downstream route layers should consume this parent-owned data. They should not
rediscover or reconstruct universal axis factors as private report payloads.

## Module Boundary Direction

This design should preserve the current conceptual split, but "parent layer"
should map to a Julia module boundary. Parent-owned factors should live in the
`CartesianParentGaussletBases` module namespace, even if their implementation
is split into separate files for readability.

```text
cartesian_cpb/
  pure coordinate-product-box geometry

CartesianParentGaussletBases.jl
  route-neutral Cartesian parent identity
  parent-owned one-body and Coulomb axis factor packets
  parent factor provenance, convention, and index-domain summaries

future CartesianCPBBlockProviders
  depends on CPB geometry plus CartesianParentGaussletBases
  returns CPB-local axis blocks, slices, factored blocks, or optional dense
  local blocks

route helpers and materializers
  consume parent factors and CPB block providers
```

Parent-dependent code should not move into `CartesianCPB`. CPB geometry remains
the coordinate-window language. Parent factors are part of the parent module;
CPB block providers are the layer above both parent factors and CPB geometry.

The intended namespace direction is:

```julia
GaussletBases.CartesianParentGaussletBases.parent_overlap_axis_factor_packet(
    parent,
    axis_bundle,
)

GaussletBases.CartesianCPBBlockProviders.cpb_interval_pair(
    parent,
    left_cpb,
    right_cpb,
)
```

The old spelling with a sibling `CartesianParentAxisFactors` namespace should
be treated as a temporary implementation shape, not the long-term module
contract.

## Proposed Ownership Layers

The intended split is:

```text
CartesianParentGaussletBasis3D
  parent lattice identity:
    axes
    parent_box
    axis_sharing
    metadata

CartesianParentAxisFactorPacket3D
  universal axis-level factors tied to that parent:
    one-body factors
    electron-electron Gaussian factors
    electron-nuclear Gaussian factors
    factor index-domain and convention labels
    provenance and availability summaries

CartesianCPB
  coordinate-window geometry:
    integer product boxes
    codimension
    boundary-stratum decomposition
    support counts

CartesianCPBBlockProvider
  bound operations:
    parent factors + CPB windows -> local axis blocks
    parent factors + CPB windows -> optional local dense blocks

Shellification / PQS / retained-unit / pair-block layers
  choose owned geometry, source boxes, retained rules, transforms, pair plans,
  and placement
```

`CartesianCPB` should remain pure geometry. It should not acquire parent
factors, operator terms, Hamiltonian concepts, route reports, or retained-space
state. The CPB block-provider layer is the function group that sits above both
CPB geometry and parent-owned factors.

## Parent Axis Factor Packet

A possible top-level shape is:

```julia
struct CartesianParentAxisFactorPacket3D{P,O,C,M}
    parent::P
    one_body::O
    coulomb::C
    metadata::M
end
```

This object should satisfy:

```julia
packet.parent === parent_basis
```

when it is built directly for a parent identity object. If future construction
requires copying or rehydrating parent objects, the packet should still carry a
compact parent fingerprint and explicit provenance so consumers can check that
the factors match the parent axes and mappings they are using.

The packet is not a retained basis, not a source-box route, not a Hamiltonian
object, and not a dense 3D operator store. It owns universal 1D factor data and
summaries about how those factors were produced.

Factor packets may be partial. A parent packet can have overlap factors
available while kinetic, coordinate moments, Coulomb expansion factors, or
nuclear factors are unavailable or not requested. Availability should be
reported by category rather than inferred from a single packet-level status.

### Dependency Classes

Not every factor in the parent layer has the same dependency footprint. The
packet metadata should distinguish at least:

```text
parent_only
  depends on parent axes and mappings/distortions only
  examples: overlap, kinetic, position, x2

parent_plus_expansion
  depends on parent axes plus a Gaussian expansion policy
  example: electron-electron Gaussian axis factors

parent_plus_nuclear_configuration
  depends on parent axes, Gaussian expansion policy, nuclei, charges, and
  nuclear locations
  example: electron-nuclear Gaussian axis factors
```

All of these can be parent-layer data because none depends on shellification,
source boxes, retained units, pair plans, or placement. They should not be
treated as equally reusable, though: changing nuclei or Gaussian expansion
policy invalidates some packet parts but not parent-only one-body factors.

### One-Body Subobject

The one-body part can be represented conceptually as:

```julia
struct CartesianParentOneBodyAxisFactors3D{S,T,X,X2,M}
    overlap_1d::S
    kinetic_1d::T
    position_1d::X
    x2_1d::X2
    metadata::M
end
```

The axis-factor values are named by axis:

```julia
overlap_1d  = (x = Sx,  y = Sy,  z = Sz)
kinetic_1d  = (x = Tx,  y = Ty,  z = Tz)
position_1d = (x = Xx,  y = Xy,  z = Xz)
x2_1d       = (x = X2x, y = X2y, z = X2z)
```

The full parent 3D one-body operators remain formal factorized expressions:

```text
S3 = kron(Sx, Sy, Sz)

T3 = kron(Tx, Sy, Sz)
   + kron(Sx, Ty, Sz)
   + kron(Sx, Sy, Tz)

X3 = kron(Xx, Sy, Sz)
Y3 = kron(Sx, Xy, Sz)
Z3 = kron(Sx, Sy, Xz)
```

The parent factor packet should not materialize these full 3D parent matrices
as routine working data. Later CPB block providers and retained transforms can
slice and contract the axis factors.

### Coulomb Subobject

The Coulomb part should cover electron-electron and electron-nuclear Gaussian
factor data without implying IDA, MWG, or retained-basis semantics:

```julia
struct CartesianParentCoulombAxisFactors3D{E,N,G,M}
    gaussian_expansion::G
    electron_electron_1d::E
    electron_nuclear_1d::N
    metadata::M
end
```

For electron-electron Coulomb, if:

```text
1 / |r1 - r2| ~= sum_g c_g exp(-alpha_g |r1 - r2|^2)
```

then each Gaussian term factorizes by axis:

```text
exp(-alpha_g |r1-r2|^2)
 =
 exp(-alpha_g (x1-x2)^2)
 exp(-alpha_g (y1-y2)^2)
 exp(-alpha_g (z1-z2)^2)
```

The parent layer can own:

```julia
electron_electron_1d[g].x
electron_electron_1d[g].y
electron_electron_1d[g].z
```

where each axis factor is a parent-axis pair factor for Gaussian term `g`.
This is not a full 3D four-index Coulomb tensor. It is universal parent-axis
factor data that later source boxes and retained transforms can consume.

For electron-nuclear Coulomb:

```text
-Z_A / |r - R_A|
  ~= -Z_A sum_g c_g
      exp(-alpha_g (x-R_Ax)^2)
      exp(-alpha_g (y-R_Ay)^2)
      exp(-alpha_g (z-R_Az)^2)
```

The natural parent-layer data has both nucleus and Gaussian-term indices:

```julia
electron_nuclear_1d[nucleus, g].x
electron_nuclear_1d[nucleus, g].y
electron_nuclear_1d[nucleus, g].z
```

These factors depend on parent axes, nuclear positions/charges, and Gaussian
expansion policy. They do not depend on shellification, PQS retained rules, or
pair-block placement.

Electron-nuclear factors belong naturally in the parent Coulomb or potential
factor packet because their source is a Gaussian expansion of a Coulomb
potential. Downstream, however, electron-nuclear attraction is consumed by a
one-body CPB block provider: one bra CPB and one ket CPB. Electron-electron
factors are consumed by a two-body CPB block provider: two bra/ket CPB pairs.
Keeping that arity distinction explicit is required for later API names and
tests.

### Provenance and Metadata

A parent factor packet should carry compact metadata such as:

- `object_kind`;
- parent fingerprint or parent object identity;
- axis counts;
- axis basis family and mapping summary;
- factor space;
- factor convention;
- normalization convention;
- axis order;
- bra/ket ordering;
- factor backend, for example PGDG analytic integrals;
- Gaussian expansion identity and parameters for Coulomb factors;
- nucleus labels, positions, and charges for electron-nuclear factors;
- factor availability by category;
- materialization status, including explicit `full_3d_parent_matrices = false`
  when appropriate;
- nonclaim flags for route adoption, retained-basis IDA semantics, Hamiltonian
  assembly, exports, and artifacts.

The metadata should be summary-like. It should not become a flat cloud of
scalar report aliases copied through every route stage.

The convention labels are not optional. The same-looking 1D matrix can live in
raw parent axis space, localized PGDG intermediate space, dictionary-basis
space, density-normalized convention, raw-weighted convention, or another
future convention. Immediate overlap facts already use labels of this form:

```text
factor_space = :parent_axis_bundle_pgdg_intermediate
factor_convention = :axis_bundle_one_body_overlap
```

Every parent factor packet and CPB block result should preserve comparable
labels so downstream code does not infer semantics from field names alone.

### Index-Space Contract

Every factor must also state the index domain its rows, columns, or tensor axes
use. CPB block providers can slice a factor by CPB intervals only when the
factor is in a CPB-addressable parent-axis index domain.

This is separate from mathematical factor space. For example:

```text
factor_space = :parent_axis_bundle_pgdg_intermediate
factor_convention = :axis_bundle_one_body_overlap
index_domain = :parent_axis_indices
axis_order = (:x, :y, :z)
bra_ket_order = (:bra, :ket)
```

The important rule:

```text
CPB interval slicing is valid only for factors whose index domain is compatible
with parent axis indices.
```

A transformed/localized axis matrix may have the right dimensions but still be
unsafe to slice by raw CPB intervals if its rows and columns no longer
correspond to parent-axis indices. Such factors need either an explicit
index-map object or a different consumer. Do not infer sliceability from matrix
size.

CPB block providers must treat this as a mandatory pre-slice gate. A provider
must not slice parent factors unless the factor packet summary explicitly says
the factors are CPB-sliceable and parent-axis-indexed. For the current
overlap-only packet, the intended labels are:

```text
sliceable_by_cpb = true
index_domain = :parent_axis_indices
index_domain_source = :axis_bundle_contract
index_domain_status =
  :assumed_parent_axis_indexed_by_current_axis_bundle_contract
sliceability_source = :index_domain_contract
sliceability_status = :sliceable_by_cpb_parent_axis_index_contract
```

If any required label is missing or incompatible, the CPB provider should block
before slicing. This preserves the distinction between "has the right shape"
and "is contractually indexed by parent axis coordinates."

For two-electron factors, the index-domain contract must name all particle and
bra/ket axes. A future tensor-like electron-electron factor should carry a
domain summary comparable to:

```text
index_domain = :parent_axis_pair_indices
particle_order = (:electron1, :electron2)
bra_ket_order = (:bra1, :ket1, :bra2, :ket2)
axis_order = (:x, :y, :z)
```

The exact labels can change, but the information cannot be omitted.

## CPB Block-Provider Layer

A CPB is only coordinate-window geometry. A CPB block provider binds a parent
factor packet to operations over CPB pairs:

```julia
struct CartesianCPBBlockProvider{P,F,M}
    parent::P
    parent_factors::F
    metadata::M
end
```

The preferred module name is `CartesianCPBBlockProviders`. "Provider" is more
accurate than "matrix kernel" for this layer because the functions may return
views, copied slices, factored axis blocks, or dense local blocks. Dense matrix
materialization is only one possible service, not the defining behavior.

The first implementation should prefer a stateless, deterministic block
provider. Cache policy can be added as a wrapper:

```julia
struct CachedCartesianCPBBlockProvider{K,C,M}
    provider::K
    cache::C
    metadata::M
end
```

Keeping cache policy out of the first scientific contract makes early tests
smaller and avoids treating performance cache behavior as route semantics.

### Basic Operation

A CPB pair supplies integer windows:

```text
left_cpb  -> IxL, IyL, IzL
right_cpb -> IxR, IyR, IzR
```

The provider is the first layer that sees both parent identity and CPB geometry.
It must validate that every CPB interval lies inside the parent box. This
validation belongs here, not in `CartesianCPB`, because a bare CPB has no
parent.

For one-body terms and electron-nuclear attraction, the provider knows how to
slice parent-owned factors:

```text
Sx[IxL, IxR],  Sy[IyL, IyR],  Sz[IzL, IzR]
Tx[IxL, IxR],  Ty[IyL, IyR],  Tz[IzL, IzR]
Xx[IxL, IxR],  Xy[IyL, IyR],  Xz[IzL, IzR]
X2x[IxL, IxR], X2y[IyL, IyR], X2z[IzL, IzR]

Vx[A,g][IxL, IxR], Vy[A,g][IyL, IyR], Vz[A,g][IzL, IzR]
```

Electron-electron factors have different arity. A two-electron kernel needs two
bra/ket CPB pairs:

```text
electron 1: bra_cpb_1 -> IaB, ket_cpb_1 -> IaK
electron 2: bra_cpb_2 -> IbB, ket_cpb_2 -> IbK
```

For Gaussian term `g`, the axis factors are sliced with four one-particle
windows or an equivalent structured pair-pair object:

```text
Gx[g][Ix1B, Ix1K, Ix2B, Ix2K]
Gy[g][Iy1B, Iy1K, Iy2B, Iy2K]
Gz[g][Iz1B, Iz1K, Iz2B, Iz2K]
```

The exact storage may be matrix-like, tensor-like, or lazy, depending on the
Gaussian factor backend. The API must not squeeze electron-electron Coulomb
into the one-body `(left_cpb, right_cpb)` shape.

For one-body overlap, the axis-block result is simply:

```julia
(x = Sx_slice, y = Sy_slice, z = Sz_slice)
```

For kinetic, the result should preserve the sum-of-axis structure:

```julia
(
    x = (kinetic = Tx_slice, overlap_y = Sy_slice, overlap_z = Sz_slice),
    y = (overlap_x = Sx_slice, kinetic = Ty_slice, overlap_z = Sz_slice),
    z = (overlap_x = Sx_slice, overlap_y = Sy_slice, kinetic = Tz_slice),
)
```

The exact representation can be refined, but the key contract is that kinetic
is not a single arbitrary dense 3D object at this layer. It is an axis-factor
sum that can be materialized locally only when requested.

### Local Dense Ordering

If a CPB provider materializes a dense local block, the local indexing order
must be explicit. The default should be parent-compatible local flattening:

```text
local x index is slowest
local z index is fastest
```

This matches `parent_flat_index` ordering restricted to the CPB intervals:

```text
(ix - first(Ix)) * length(Iy) * length(Iz)
+ (iy - first(Iy)) * length(Iz)
+ (iz - first(Iz))
```

using 1-based Julia indexing after adding one as needed. Retained-unit
transforms and column ranges should not have to guess the local dense ordering.

Dense CPB blocks are still local CPB product-space operator blocks. They are
not retained blocks, not global matrices, not Hamiltonian objects, and not
evidence of route adoption. Retained-unit transforms, PQS shell realization,
Lowdin cleanup, and global placement remain downstream responsibilities.

CPB provider summaries should remain compact. If a provider materializes a
dense local block, the numerical matrix belongs on the returned object field,
for example:

```julia
dense_block.dense_block
```

The summary metadata should not duplicate that matrix. It should carry compact
facts such as:

```text
dense_block_available
dense_block_shape
dense_block_eltype
dense_local_block_materialized
```

plus provenance, factor-space, convention, index-domain, local-ordering, and
nonclaim labels. This keeps summaries usable as fingerprints without turning
metadata into a second numerical payload.

### Possible API

The provider API should start small and explicit:

```julia
cpb_axis_interval_pair(left::CoordinateProductBox, right::CoordinateProductBox)

cpb_one_body_axis_blocks(provider, left, right, term)
cpb_one_body_block(provider, left, right, term)

cpb_electron_nuclear_axis_blocks(
    provider,
    bra_cpb,
    ket_cpb;
    nucleus,
    gaussian_index,
)

cpb_electron_electron_axis_blocks(
    provider,
    bra_cpb_1,
    ket_cpb_1,
    bra_cpb_2,
    ket_cpb_2;
    gaussian_index,
)
```

Later generic dispatch can exist, but the first implementation should avoid a
large catch-all API that hides term-specific contracts.

For two-electron work, a structured helper such as `CPBPairPair3D` or
`CPBTwoParticleWindow3D` may be better than a long argument list. The important
contract is the same: electron-electron axis blocks are indexed by two
bra/ket CPB pairs, not by one CPB pair.

### Returned Objects

Returned values should be structured objects or compact named tuples, not
anonymous field clouds. A one-body axis-block object might look like:

```julia
struct CPBOneBodyAxisBlockSet{B,M}
    term::Symbol
    left_cpb::CoordinateProductBox
    right_cpb::CoordinateProductBox
    axis_blocks::B
    metadata::M
end
```

Metadata should include:

- source parent factor packet fingerprint;
- left/right CPB intervals and shapes;
- local dense ordering, when materialized;
- term;
- factor index domain and any required index maps;
- whether blocks are views or copied matrices;
- whether a dense local CPB block has been materialized;
- factor backend/provenance;
- factor space, convention, normalization, axis order, and bra/ket ordering;
- no-go flags for Hamiltonian assembly, retained transforms, route adoption,
  Coulomb/IDA/MWG where irrelevant, exports, and artifacts.

Metadata should not duplicate large numerical payloads. Dense local matrices
should live only on the owning dense-block object, while summaries record
availability, shape, element type, provenance, and convention labels.

### Local Block Records and Collections

A local block record is the bridge between one CPB provider output and later
placement or assembly layers. For overlap, the current names are:

```julia
CPBLocalOverlapBlockRecord
CPBLocalOverlapBlockCollection
cpb_local_overlap_block_record
cpb_local_overlap_block_collection
```

A record may reference either a `CPBOverlapDenseBlock` or a
`CPBOverlapAxisBlockSet` as `source_block`. The source block is kept by
reference. Record and collection summaries must remain compact fingerprints and
must not duplicate dense matrix payloads. Dense numerical data stays on the
owning block object, for example `dense_block.dense_block`.

Record metadata should include:

- `term = :overlap`;
- `block_key`;
- `source_kind`;
- compact left/right CPB summaries;
- interval-pair summary;
- axis-block summary;
- dense-block summary when materialized;
- dense availability and shape;
- factor-space, convention, normalization, and index-domain labels;
- `placement_status = :unassigned`;
- `retained_transform_status = :unassigned`;
- `route_driver_wiring = false`;
- `global_matrix_materialized = false`.

Collection metadata should include:

- `record_count`;
- `terms`;
- `block_keys`;
- `dense_block_count`;
- `blocked_record_count`;
- `blocked_record_blockers`;
- `placement_status = :unassigned`;
- `retained_transform_status = :unassigned`;
- route/global nonclaim flags, including `route_driver_wiring = false` and
  `global_matrix_materialized = false`.

This collection layer is still local CPB product-space overlap metadata. It is
not retained placement, not route-global overlap, and not Hamiltonian assembly.
It only groups local CPB provider outputs so a later, separately reviewed layer
can decide whether and how to assign placement or retained transforms.

### Placement Contract Boundary

The next boundary after a local CPB overlap block collection is
placement/retained-transform assignment. Local collection availability is not
global overlap availability. A collection says "local CPB product-space overlap
blocks exist"; it does not say where those blocks belong in a retained or
global matrix, how they should be transformed, or how multiple contributions
should be accumulated.

A future placement adapter must require explicit facts such as:

- local block collection;
- left/right retained transforms or retained-unit transforms;
- left/right retained column ranges or placement ranges;
- local-to-retained ordering convention;
- global or retained dimension;
- placement plan;
- accumulation rule for multiple local blocks.

Expected blockers at this boundary include:

- `:missing_local_overlap_collection`;
- `:missing_retained_transform`;
- `:missing_left_column_range`;
- `:missing_right_column_range`;
- `:missing_global_dimension`;
- `:missing_placement_plan`;
- `:missing_accumulation_rule`.

Until those facts are present, route-global overlap remains unavailable even if
the CPB provider has produced local overlap records or a local overlap
collection. This boundary is separate from Hamiltonian assembly and should not
claim Hamiltonian readiness, IDA/MWG semantics, PQS Lowdin/projection, exports,
or artifacts.

The private overlap local collection adapter is a boundary/fingerprint object
for this state. It recognizes `CPBLocalOverlapBlockCollection` as structured
local overlap source data, but it is not global overlap input-facts
availability and not route-global overlap stage adoption. The expected current
state is:

```text
local collection available
global overlap blocked
blocker = :missing_placement_or_retained_transform
private_global_overlap_input_facts_available = false
route_global_overlap_stage_source = false
global_matrix_materialized = false
route_driver_wiring = false
```

The adapter should remain a fingerprint until a later placement layer supplies
the required retained-transform, column-range, dimension, placement-plan, and
accumulation-rule facts. It must not become a placement engine by accumulating
local CPB blocks into a route/global matrix.

A placement candidate is the next metadata-only status carrier. It may combine
the local collection adapter or placement-requirements fingerprint with
optional placement facts, but it is still not a placement plan and not global
assembly. Its purpose is to record which requirements are available and which
are still missing, without inventing ad hoc scalar placement fields in later
adapters. With the current local collection and no placement facts, the
candidate should report only `:local_cpb_overlap_collection` available and keep
the retained transform, column ranges, global dimension, placement plan, and
accumulation rule as missing. It must keep route/global nonclaim flags false.

This lets retained-unit and pair-block code choose whether to use axis blocks
directly, materialize a dense local CPB block, apply left/right transforms, or
place into a global matrix.

## Relation to PQS Source Boxes

PQS should remain source-box-first. The source CPB is the compact product-type
parent block where operator information is gathered before shell projection and
Lowdin cleanup. The CPB block-provider layer is a natural way to obtain
source-box operator blocks from parent factors:

```text
parent factor packet
-> CPB block provider over source CPB pairs
-> source-box local operator blocks
-> retained-rule transforms
-> later shell realization / Lowdin where required
```

Support-row or shell-row contractions remain oracle/debug surfaces unless a
reviewed framework decision promotes them. The CPB block-provider layer should
not make support-row contraction the PQS algorithm.

## Relation to Current Private Overlap Facts

Current private overlap input facts are a bridge:

```text
driver report fields
-> private overlap input-facts helper
-> route-global overlap adapter
```

The real-report fingerprints show why this bridge is temporary. A probe-enabled
dry report shows that a structured axis-bundle object may exist under
`route_materializer_payload.parent_axis_bundle_object`, and the current private
helper can read that object through a narrow compatibility source path. That is
useful as a fingerprint. It should not make the route materializer payload the
natural authority for universal parent-axis factors.

The intended replacement direction is:

```text
parent basis
-> parent axis factor packet
-> CPB block provider
-> local CPB source/operator blocks
-> pair-block materialization
-> local block collection
-> optional route/global placement
```

The private overlap facts helper should eventually become a compatibility
adapter or disappear. It should not expand into a broad scalar report-field
interface for every operator term.

The current overlap-only implementation already has a local CPB product-space
path:

```text
CartesianParentGaussletBases.parent_overlap_axis_factor_packet
-> CartesianCPBBlockProviders.cpb_interval_pair
-> CartesianCPBBlockProviders.cpb_overlap_axis_blocks
-> CartesianCPBBlockProviders.cpb_overlap_dense_block
-> CartesianCPBBlockProviders.cpb_local_overlap_block_record
-> CartesianCPBBlockProviders.cpb_local_overlap_block_collection
```

That path can validate CPB interval pairs, slice parent-owned overlap axis
factors, optionally materialize a local dense CPB overlap block, and group local
provider outputs into compact collection metadata. It remains local CPB overlap
only. It is not route/global overlap adoption, does not place a global matrix,
does not assign retained transforms, and does not change private route-driver
overlap behavior.

## What This Does Not Claim

This note does not claim:

- a public API is ready;
- a production Hamiltonian route exists;
- full White-Lindsey route assembly is complete;
- PQS shell projection or Lowdin cleanup is implemented;
- Coulomb blocks are ready for production use;
- electron-electron factors have IDA or MWG semantics;
- GTO supplements are ordinary gausslet quadrature points;
- dense 3D parent matrices should be materialized routinely;
- exports or artifacts should be written from the parent factor layer.

The parent factor packet owns universal parent-axis factors. It does not own
retained-basis semantics, IDA division, shell/PQS realization, or global matrix
assembly.

## Suggested Migration Path

The migration should be incremental:

0. Keep one-body CPB-pair providers distinct from electron-electron
   CPB-pair-pair providers in names, tests, and returned objects.
1. Define a parent-layer contract for one-body axis factors, starting with
   overlap, provenance, factor convention, and index-domain labels. Map the
   existing `parent_axis_bundle_object` seed to that contract.
2. Add a small private parent factor packet object or summary for overlap-only
   facts. Keep it tied to `CartesianParentGaussletBasis3D`.
3. Add CPB interval-pair helpers that return parent-axis slice ranges for CPB
   pairs without operator semantics.
4. Add a private overlap CPB axis-block provider over the parent factor packet.
5. Repoint the private overlap input-facts path toward the parent factor packet
   or CPB provider, while preserving current real-report fingerprint tests.
6. Extend one-body factors to kinetic, position, and x2 only after overlap has
   a clean parent-factor and CPB-provider contract.
7. Treat electron-nuclear and electron-electron factor packets as a later
   design step with explicit Gaussian-expansion provenance and no IDA/MWG
   overclaims.

At each step, prefer compact objects and summaries over report-field aliases.
When a route still needs compatibility fields, derive them at the report
boundary rather than carrying them through every stage.

## Testing Strategy

Early tests should be narrow and factual:

- parent factor packet availability and provenance;
- factor convention and index-domain labels;
- parent identity/fingerprint match;
- partial packet availability, for example overlap available while Coulomb is
  not requested;
- CPB interval slicing for small boxes;
- CPB parent-box validation failures for out-of-parent intervals;
- local dense ordering for one small CPB block;
- overlap axis-block slices against known parent 1D factors;
- local dense overlap block materialization for one small CPB pair, if needed;
- two-electron CPB pair-pair arity, initially metadata-only;
- private report compatibility fingerprints that state what is still missing.

Avoid broad nested/integration tests for each small contract pass. Performance
validation becomes necessary when a provider starts serving production-like
routes or public APIs.

## Open Questions

- Should the first object be named `CartesianParentAxisFactorPacket3D`,
  `CartesianParentAxisOperatorBundle3D`, or something closer to the existing
  `parent_axis_bundle_object` language?
- Should one-body and Coulomb factor packets live in one module or separate
  submodules under the parent layer?
- Should electron-nuclear factors live under a Coulomb factor packet, a
  potential factor packet, or a one-body factor packet while still preserving
  their one-body CPB-provider consumption path?
- Which factor values should be views, copied matrices, or lazy objects?
- What parent fingerprint is stable enough for packet/parent compatibility
  checks without comparing large objects?
- How should axis-factor caches be invalidated or shared across repeated CPB
  provider calls?
- Should the CPB provider own a cache, or should caching be a wrapper around a
  stateless block provider?
- Where should electron-nuclear Gaussian expansion policy live when different
  nuclei or ECP/reference paths need different approximations?
- How should future CR2 exports refer to parent factors without defining repo
  API prematurely?
- What should the two-electron window object be called:
  `CPBPairPair3D`, `CPBTwoParticleWindow3D`, or something more explicit?

These questions should be resolved by small contracts and focused tests, not by
expanding route-driver report fields.
