# Cartesian Parent Factors and CPB Matrix Kernels

Date: 2026-06-09

This note records a proposed architectural direction for Cartesian parent
operator data, coordinate-product-box local blocks, and the staged PQS /
White-Lindsey route. It is a design note, not an implementation contract yet.

The immediate motivation is the private global-overlap input-facts work. Recent
tests show two different real dry-run report states:

- a manual-count dry report can carry parent axis counts, but has no parent
  axis-bundle factors and correctly blocks on
  `:missing_parent_axis_bundle_overlap_factors`;
- a probe-enabled dry report carries
  `route_materializer_payload.parent_axis_bundle_object`, and the private
  overlap facts helper can read its 1D overlap factors.

That is useful evidence, but it should not turn route-driver report payloads
into the long-term owner of parent operator data. The better direction is to
promote the underlying idea into a parent-owned axis factor packet and a CPB
matrix-kernel layer.

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
    provenance and availability summaries

CartesianCPB
  coordinate-window geometry:
    integer product boxes
    codimension
    boundary-stratum decomposition
    support counts

CartesianCPBBlockProvider / CartesianCPBMatrixKernel
  bound operations:
    parent factors + CPB windows -> local axis blocks
    parent factors + CPB windows -> optional local dense blocks

Shellification / PQS / retained-unit / pair-block layers
  choose owned geometry, source boxes, retained rules, transforms, pair plans,
  and placement
```

`CartesianCPB` should remain pure geometry. It should not acquire parent
factors, operator terms, Hamiltonian concepts, route reports, or retained-space
state. The CPB block-provider layer is the function group that sits above CPB
geometry.

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
as routine working data. Later CPB kernels and retained transforms can slice
and contract the axis factors.

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

### Provenance and Metadata

A parent factor packet should carry compact metadata such as:

- `object_kind`;
- parent fingerprint or parent object identity;
- axis counts;
- axis basis family and mapping summary;
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

## CPB Matrix-Kernel Layer

A CPB is only coordinate-window geometry. A CPB block provider or matrix kernel
binds a parent factor packet to operations over CPB pairs:

```julia
struct CartesianCPBBlockProvider{P,F,C,M}
    parent::P
    parent_factors::F
    cache::C
    metadata::M
end
```

or:

```julia
struct CartesianCPBMatrixKernel{P,F,M}
    parent::P
    parent_factors::F
    metadata::M
end
```

The name `provider` or `kernel` is intentional. This layer does not necessarily
allocate a matrix every time. It may return views, copied axis slices,
factorized axis-block objects, or dense local CPB blocks depending on the API.

### Basic Operation

A CPB pair supplies integer windows:

```text
left_cpb  -> IxL, IyL, IzL
right_cpb -> IxR, IyR, IzR
```

The provider knows how to slice parent-owned factors:

```text
Sx[IxL, IxR],  Sy[IyL, IyR],  Sz[IzL, IzR]
Tx[IxL, IxR],  Ty[IyL, IyR],  Tz[IzL, IzR]
Xx[IxL, IxR],  Xy[IyL, IyR],  Xz[IzL, IzR]
X2x[IxL, IxR], X2y[IyL, IyR], X2z[IzL, IzR]

Gx[g][IxL, IxR], Gy[g][IyL, IyR], Gz[g][IzL, IzR]
Vx[A,g][IxL, IxR], Vy[A,g][IyL, IyR], Vz[A,g][IzL, IzR]
```

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

### Possible API

The provider API should start small and explicit:

```julia
cpb_axis_interval_pair(left::CoordinateProductBox, right::CoordinateProductBox)

cpb_one_body_axis_blocks(provider, left, right, term)
cpb_one_body_block(provider, left, right, term)

cpb_electron_nuclear_axis_blocks(
    provider,
    left,
    right;
    nucleus,
    gaussian_index,
)

cpb_electron_electron_axis_blocks(
    provider,
    left,
    right;
    gaussian_index,
)
```

Later generic dispatch can exist, but the first implementation should avoid a
large catch-all API that hides term-specific contracts.

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
- term;
- whether blocks are views or copied matrices;
- whether a dense local CPB block has been materialized;
- factor backend/provenance;
- no-go flags for Hamiltonian assembly, retained transforms, route adoption,
  Coulomb/IDA/MWG where irrelevant, exports, and artifacts.

This lets retained-unit and pair-block code choose whether to use axis blocks
directly, materialize a dense local CPB block, apply left/right transforms, or
place into a global matrix.

## Relation to PQS Source Boxes

PQS should remain source-box-first. The source CPB is the compact product-type
parent block where operator information is gathered before shell projection and
Lowdin cleanup. The CPB kernel layer is a natural way to obtain source-box
operator blocks from parent factors:

```text
parent factor packet
-> CPB kernel over source CPB pairs
-> source-box local operator blocks
-> retained-rule transforms
-> later shell realization / Lowdin where required
```

Support-row or shell-row contractions remain oracle/debug surfaces unless a
reviewed framework decision promotes them. The CPB kernel layer should not make
support-row contraction the PQS algorithm.

## Relation to Current Private Overlap Facts

Current private overlap input facts are a bridge:

```text
driver report fields
-> private overlap input-facts helper
-> route-global overlap adapter
```

The real-report fingerprints show why this bridge is temporary. The helper can
read overlap factors when a route materializer payload happens to carry a
parent axis-bundle object, but that payload is not the natural authority for
universal parent-axis factors.

The intended replacement direction is:

```text
parent basis
-> parent axis factor packet
-> CPB block provider
-> pair-block materialization
-> optional route/global placement
```

The private overlap facts helper should eventually become a compatibility
adapter or disappear. It should not expand into a broad scalar report-field
interface for every operator term.

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

1. Define a parent-layer contract for one-body axis factors, starting with
   overlap and provenance. Map the existing `parent_axis_bundle_object` seed to
   that contract.
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
- parent identity/fingerprint match;
- CPB interval slicing for small boxes;
- overlap axis-block slices against known parent 1D factors;
- local dense overlap block materialization for one small CPB pair, if needed;
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
- Which factor values should be views, copied matrices, or lazy objects?
- What parent fingerprint is stable enough for packet/parent compatibility
  checks without comparing large objects?
- How should axis-factor caches be invalidated or shared across repeated CPB
  provider calls?
- Should the CPB provider own a cache, or should caching be a wrapper around a
  stateless matrix kernel?
- Where should electron-nuclear Gaussian expansion policy live when different
  nuclei or ECP/reference paths need different approximations?
- How should future CR2 exports refer to parent factors without defining repo
  API prematurely?

These questions should be resolved by small contracts and focused tests, not by
expanding route-driver report fields.
