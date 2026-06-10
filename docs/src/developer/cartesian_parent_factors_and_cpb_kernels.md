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

## Architecture Pivot: CPB-Local Operator Blocks

The overlap placement pilot should now be read as evidence for a CPB-local
operator-block layer, not as route-global placement machinery. The useful
organizing split is:

```text
parent / axis-factor layer
  owns universal 1D factors tied to the Cartesian parent

CPB operator layer
  builds local rectangular CPB operator blocks for source CPB pairs
  may materialize dense provider-level local matrices
  remains realization-neutral and route-neutral

realization layer
  White-Lindsey consumes CPB blocks nearly directly
  PQS applies local retained transforms:
    O_retained = T_left' * O_cpb * T_right

route/global layer
  assigns retained/global ranges
  accumulates route/global matrices later
```

The CPB layer is the right place to fill local operator blocks because its
geometry is rectangular, local, and parent-axis-addressable. It sees both the
parent-axis factor contract and the left/right CPB intervals, but it does not
own retained-unit choices, shell realization, Lowdin cleanup, pair-block
placement, or global matrix ranges.

Overlap is the current pilot term. Kinetic, position, and x2 should later
become analogous one-body CPB operator blocks over the same geometry and
axis-factor contracts. Coulomb and other two-body terms should also remain
CPB-local at this layer, but they may require pair-pair windows or factorized
records rather than the one-body `(left_cpb, right_cpb)` shape.

The CPB operator layer should keep separate kernel families instead of forcing
every term through one universal product kernel:

1. Axis-product and sum-of-axis-products kernels cover simple separable
   one-body terms such as overlap, kinetic pieces, position, and x2.
2. Electron-nuclear Gaussian-sum kernels cover CPB-local one-body potential
   blocks with the Gaussian expansion index as an inner loop.
3. Electron-electron Gaussian-sum kernels cover CPB-local two-body or
   pair-pair/factorized records, also with the Gaussian expansion index as an
   inner loop.

Coulomb-family blocks should not be built by an outer loop that repeatedly
calls `cpb_axis_product_operator_block` for each Gaussian alpha. Their common
surface with one-body blocks should eventually be compact result records and
metadata, not necessarily one universal kernel implementation. For separable
one-body terms, inactive directions still use explicit overlap factors, not
ambiguous identity labels.

The synthetic numerical pilot demonstrates local retained-transform mechanics:

```text
O_retained = T_left' * O_cpb * T_right
```

It does not demonstrate route-global overlap adoption. Do not continue adding
placement metadata layers before deciding how real White-Lindsey and PQS
realization should consume CPB operator blocks and how real retained transforms
and ranges should be produced.

## Mixed Gausslet/GTO Supplement Port Map

The mixed gausslet/GTO supplement path should be ported into the CPB operator
layer by adapting the existing polynomial-Gaussian integral code. Do not
rederive shell formulas in CPB providers.

### Existing Function Surfaces

The exact supplement representation is
`CartesianGaussianShellSupplementRepresentation3D` in
`src/cartesian_basis_representation.jl`. Each
`CartesianGaussianShellOrbitalRepresentation3D` carries:

- `label`
- `angular_powers`
- `center`
- `exponents`
- `coefficients`
- `primitive_normalization`

The active normalization convention is
`:axiswise_normalized_cartesian_gaussian`. This convention is part of the
mixed-kernel metadata contract; CPB adapters should not silently accept another
normalization.

The current mixed overlap/GTO projection surfaces are:

- `src/cartesian_cross_overlap.jl`
  - `_cartesian_basis_supplement_axis_primitive_cross` builds one-axis
    primitive overlap tables between a Cartesian basis axis and one supplement
    orbital using `_qwrg_atomic_basic_integral`.
  - `_cartesian_basis_supplement_cross` contracts those axis tables into a 3D
    Cartesian-by-supplement overlap block.
  - `_cartesian_supplement_orbital_axis_overlap_matrix` and
    `_cartesian_supplement_cross_overlap` build supplement-by-supplement
    overlap blocks.
- `src/cartesian_gto_probes.jl`
  - `_pqs_source_box_gto_axis_projection` slices/projects one-axis primitive
    supplement overlaps for a raw product-box interval.
  - `_pqs_source_box_gto_cross_overlap_shadow` forms a source-box/GTO overlap
    diagnostic from those projected axis tables. It is a shadow/fingerprint
    route, not CPB operator adoption.
  - `gto_overlap_matrix` and `_cartesian_final_gto_cross_overlap_handoff`
    remain final-basis/GTO diagnostic and handoff surfaces.

The current Qiu-White/raw mixed one-body and Gaussian-factor surfaces are:

- `src/ordinary_qw_raw_blocks.jl`
  - `_qwrg_atomic_shell_prefactor`, `_qwrg_atomic_basic_integral`, and
    `_qwrg_atomic_kinetic_integral` are wrappers around
    `GaussianAnalyticIntegrals` polynomial-Gaussian formulas.
  - `_qwrg_atomic_axis_cross_data` builds one-axis gausslet-proxy by
    supplement-primitive tables for overlap, kinetic, position, x2, and
    centered Gaussian factor terms.
  - `_qwrg_atomic_axis_aa_data` builds the corresponding
    supplement-primitive by supplement-primitive tables.
  - `_qwrg_cartesian_shell_cross_moment_blocks_3d` contracts mixed axis tables
    into 3D `overlap_ga`, `kinetic_ga`, `position_*_ga`, `x2_*_ga`, and
    `factor_ga` blocks.
  - `_qwrg_cartesian_shell_self_moment_blocks_3d` builds the analogous
    supplement self blocks.
  - `_qwrg_atomic_axis_factor_cross_data` and
    `_qwrg_atomic_axis_factor_aa_data` build shifted per-center Gaussian factor
    tables used by molecular nuclear terms.
  - `_qwrg_diatomic_cartesian_shell_overlap_blocks_3d` and
    `_qwrg_diatomic_cartesian_shell_blocks_3d` are the existing mixed
    Cartesian-shell molecular block surfaces. The latter returns by-center
    nuclear mixed pieces as `nuclear_ga_by_center` and `nuclear_aa_by_center`.
  - `_qwrg_cross_1d_blocks_proxy` and `_qwrg_cross_1d_blocks` provide the
    current proxy-based analytic Gaussian primitive cross-block path for
    simple 1D supplement data. `_qwrg_cross_1d_blocks_midpoint` is an older
    midpoint route and should be treated as reference/debug context, not the
    preferred CPB port.
  - `_qwrg_raw_overlap_blocks`, `_qwrg_raw_kinetic_cross_block`,
    `_qwrg_raw_axis_blocks`, `_qwrg_raw_factor_cross_blocks`, and
    `_qwrg_raw_one_body_blocks` show the current product contractions for
    overlap, kinetic, position/x2, Gaussian factors, and combined one-body
    terms.

The shared formula source is `src/GaussianAnalyticIntegrals.jl`:

- `polynomial_gaussian_shell_prefactor`
- `polynomial_gaussian_basic_integral`
- `polynomial_gaussian_kinetic_integral`
- `polynomial_gaussian_pair_factor_integral`
- `centered_polynomial_gaussian_pair_factor_integral`

`polynomial_gaussian_basic_integral` already accepts `xpower`, so higher
coordinate moments are formula-supported. The existing mixed QW surfaces expose
position and x2 explicitly; x3 and higher should be added later as named
wrappers over this existing formula source only after a term-level convention is
reviewed.

For pure GTO/supplement Coulomb references, `src/gaussian_coulomb_reference.jl`
contains the small-system analytic reference path:

- `gaussian_coulomb_pair_matrix`
- `_gaussian_coulomb_pair_terms`
- `_gaussian_coulomb_pair_integral`
- `_gaussian_coulomb_pair_matrix_compressed_checked`
- `_gaussian_centered_pair_terms`

These are oracle/reference surfaces for supplement pair Coulomb behavior. They
are not a large-system CPB kernel, but they define useful polynomial-Gaussian
pair conventions.

`src/radial_ylm_gto_bridge.jl` is an input-record bridge, not an operator
kernel. `radial_ylm_fit_cartesian_gto_adapter` converts centered radial/Ylm
fits into `CartesianGaussianShellSupplementRepresentation3D` plus a coefficient
map using the same axiswise normalized Cartesian primitive convention.

### CPB-Local Reuse Boundary

The CPB-local mixed supplement layer should consume prepared source records and
existing axis kernels:

1. Build or reference mixed one-axis tables with the existing
   polynomial-Gaussian helpers.
2. Slice those tables by CPB intervals or source-box intervals.
3. Use the existing CPB operator primitives for simple one-body separable
   terms:
   - overlap as one axis-product term;
   - kinetic as `Kx Sy Sz + Sx Ky Sz + Sx Sy Kz`;
   - position as one explicit coordinate axis factor times overlap factors;
   - x2 and later coordinate moments as explicit moment factors times overlap
     factors.
4. Keep nuclear attraction as a by-center Gaussian-sum one-body kernel, using
   `_qwrg_atomic_axis_factor_cross_data` /
   `_qwrg_atomic_axis_factor_aa_data` conventions for shifted centers and
   expansion terms.
5. Keep electron-electron or Gaussian pair-factor work in the Coulomb-family
   kernel family, not in the simple axis-product one-body primitive.

Reusable CPB-local axis kernel candidates are therefore:

- mixed overlap axis tables from `_cartesian_basis_supplement_axis_primitive_cross`
  for product-box/GTO overlap projections;
- mixed overlap/kinetic/position/x2/factor axis tables from
  `_qwrg_atomic_axis_cross_data`;
- supplement self overlap/kinetic/position/x2/factor axis tables from
  `_qwrg_atomic_axis_aa_data`;
- shifted center factor tables from `_qwrg_atomic_axis_factor_cross_data` and
  `_qwrg_atomic_axis_factor_aa_data`;
- pure supplement Coulomb pair reference terms from
  `gaussian_coulomb_pair_matrix` and its private term builders for tiny oracle
  comparisons only.

The CPB layer should not consume loose dense final matrices or final-basis GTO
handoff results as authority. Those are oracle surfaces or diagnostics.

### Required Input Records

A future mixed CPB supplement adapter needs compact structured inputs:

- parent identity and parent axis counts;
- left/right CPB or source-box interval metadata, in x/y/z order;
- local ordering `:parent_compatible_x_slowest_z_fastest`;
- `CartesianGaussianShellSupplementRepresentation3D` or a compact per-orbital
  record with label, angular powers, center, exponents, coefficients, and
  primitive normalization;
- a contraction convention for primitive-to-orbital coefficients;
- a source for gausslet/proxy axis rows, for example
  `_MappedOrdinaryGausslet1DBundle` plus its auxiliary proxy layer or an
  already-carried axis table object;
- `CoulombGaussianExpansion` when Gaussian factors, nuclear attraction, or
  pair factors are requested;
- center/nucleus identity, charge, and center coordinates for electron-nuclear
  by-center terms;
- explicit term labels such as `:overlap`, `:kinetic`, `:position_x`,
  `:x2_z`, `:electron_nuclear_by_center`, or
  `:electron_electron_pair_factor`.

### Metadata And Convention Labels

Mixed CPB operator records should carry compact labels rather than scalar field
clouds:

- `source_kind = :mixed_gausslet_gto_supplement_operator_block`
- `axis_kernel_source = :existing_qw_polynomial_gaussian_axis_tables`
- `formula_source = :GaussianAnalyticIntegrals_polynomial_gaussian`
- `supplement_representation = :CartesianGaussianShellSupplementRepresentation3D`
- `primitive_normalization = :axiswise_normalized_cartesian_gaussian`
- `contraction_convention = :orbital_coefficients_contract_primitive_axis_tables`
- `center_shift_convention = :explicit_axis_center_coordinates`
- `shell_power_order = (:x, :y, :z)`
- `axis_order = (:x, :y, :z)`
- `local_ordering = :parent_compatible_x_slowest_z_fastest`
- `galerkin_operator = true` for one-body operator blocks
- `by_center = true` for electron-nuclear center records
- `centers_summed = false` unless a later summation wrapper is explicitly added
- `ida_mwg_semantics = false`
- `route_driver_wiring = false`
- `route_global_matrix_materialized = false`
- `hamiltonian_assembly = false`
- `exports_or_artifacts = false`

### First Equivalence Test Later

The first implementation test should be a tiny provider-level equivalence, not
a route/Hamiltonian test:

1. Build a small Cartesian product parent plus one
   `CartesianGaussianShellSupplementRepresentation3D` orbital.
2. Select a point, edge, or face CPB interval pair.
3. Build a CPB-local mixed overlap block from the existing axis table source and
   the CPB axis-product primitive.
4. Compare only the corresponding local submatrix against the existing
   `_cartesian_basis_supplement_cross`, `_pqs_source_box_gto_cross_overlap_shadow`,
   or `_qwrg_diatomic_cartesian_shell_overlap_blocks_3d` oracle, depending on
   the fixture.

A second tiny test can do `position_x` or `x2_x` against
`_qwrg_diatomic_cartesian_shell_blocks_3d`. Nuclear attraction should remain
by-center and compare against `nuclear_ga_by_center` /
`nuclear_aa_by_center` from that existing path. Electron-electron supplement
pair behavior should use `gaussian_coulomb_pair_matrix` or the current
White-Lindsey Coulomb matrix as an oracle, not a route-global Hamiltonian.

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
The current packet implementation keeps the overlap packet-level status stable
for compatibility and also carries one-body factors when the axis bundle
already provides structured 1D matrices. The audited real-bundle source paths
are `axis.pgdg_intermediate.kinetic`, `axis.pgdg_intermediate.position`, and
`axis.pgdg_intermediate.x2`. Structured top-level fallbacks `axis.kinetic`,
`axis.position`, and `axis.x2` are also accepted. Missing or malformed optional
one-body factors are reported by category status; they do not invalidate an
otherwise available overlap packet.

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

### Coulomb-Family CPB Kernel Plan

The existing CPB-local one-body kernels are deliberately narrow:

- `cpb_axis_product_operator_block`;
- `cpb_sum_of_axis_products_operator_block`.

They are the right primitives for overlap, kinetic pieces, position, x2, and
similar simple separable one-body terms. They are not the universal Coulomb
kernel. Coulomb-family terms should share compact CPB-local result records and
metadata with the one-body operator layer, but their numerical kernels need a
different loop structure.

The relevant `PureGaussianGausslet.jl` Coulomb routines use a Gaussian
expansion and keep the Gaussian expansion index inside the optimized local
contraction. In the electron-nuclear path, `getAtomPot` builds transformed 1D
factors indexed by Gaussian term and axis, then fills the local 3D block with
an inner sum of the form:

```text
sum_g c_g * Vx[g, ix, jx, A] * Vy[g, iy, jy, A] * Vz[g, iz, jz, A]
```

with the nuclear charge and center applied for nucleus `A`. The mixed Gaussian
addition path in `addinGaussians` uses the same pattern with atom and Gaussian
indices carried on the transformed 1D factors. That structure is important:
electron-nuclear attraction should have a separate CPB-local Gaussian-sum
one-body kernel, with inputs like:

- left CPB and right CPB;
- a parent Coulomb or nuclear factor packet;
- nucleus or center identity, charge, and coordinates;
- Gaussian expansion coefficients and exponents;
- compact convention and provenance metadata.

The output should be a local rectangular CPB product-space block. The Gaussian
index `g` must be an inner loop, vectorized contraction, or otherwise optimized
inside that kernel. It should not be implemented as an outer loop that calls
`cpb_axis_product_operator_block` once per Gaussian term and allocates one dense
CPB product matrix per term.

For electron-electron Coulomb, `gethamsNoPot` builds transformed 1D Gaussian
pair factors, applies the relevant axis weights, and fills the local Coulomb
object with the same inner Gaussian sum pattern:

```text
sum_g c_g * Gx[g, ...] * Gy[g, ...] * Gz[g, ...]
```

The CPB form has different arity from the electron-nuclear one-body block. It
should consume two CPB source pairs or a structured CPB pair-pair window. Its
output may be dense, factorized, or a compact local two-body record, depending
on the reviewed storage and contraction plan. The Gaussian expansion index
again belongs inside the specialized kernel, not outside as repeated
axis-product calls.

Parent-layer Coulomb requirements follow this split. Universal 1D
Coulomb-family factors belong with the parent and Gaussian-expansion policy:

- electron-nuclear factors depend on parent axes, nuclei, nuclear charges,
  nuclear locations, and the chosen Gaussian expansion;
- electron-electron factors depend on parent axes and the Gaussian expansion
  for the pair interaction;
- neither factor family depends on WL/PQS realization, retained units, source
  pair placement, or route-global matrix assembly.

The performance contract should imitate the highly optimized
`PureGaussianGausslet.jl` structure. That code precomputes or transforms 1D
Gaussian factors, then sums over the Gaussian expansion index in the local
contraction. Reusing that organization should give excellent performance. A
future Coulomb CPB kernel should not allocate one dense CPB product matrix per
Gaussian term unless a separate performance review demonstrates that such a
materialization strategy is better for a specific backend or cache policy.

### Current Coulomb Source Audit

The current parent axis-bundle construction carries several Coulomb-family
ingredients, but not a reviewed parent Coulomb factor packet.

The current bundle shape is:

```text
parent_axis_bundle_object::_CartesianNestedAxisBundles3D
  .bundle_x::_MappedOrdinaryGausslet1DBundle
  .bundle_y::_MappedOrdinaryGausslet1DBundle
  .bundle_z::_MappedOrdinaryGausslet1DBundle
```

Each axis bundle owns a `pgdg_intermediate` with one-body and Gaussian-factor
fields:

```text
axis.pgdg_intermediate.gaussian_factors
axis.pgdg_intermediate.gaussian_factor_terms
axis.pgdg_intermediate.pair_factor_terms_raw
axis.pgdg_intermediate.pair_factors
axis.pgdg_intermediate.pair_factor_terms
axis.pgdg_intermediate.exponents
axis.pgdg_intermediate.center
axis.pgdg_intermediate.weights
axis.pgdg_intermediate.centers
```

For the route materializer reports, the axis bundle is carried through:

```text
report.route_materializer_payload.parent_axis_bundle_object
report.route_materializer_payload.axis_bundle_backend
```

The route parent report carries structured center metadata separately:

```text
report.nuclear_charges
report.atom_locations
report.center_table
```

For basis objects, the route-neutral parent metadata or QW parent object may
also carry nuclei and nuclear charges, for example:

```text
parent.metadata.nuclei
parent.metadata.nuclear_charges
report.route_materializer_payload.parent_qw_basis_object.nuclei
report.route_materializer_payload.parent_qw_basis_object.nuclear_charges
```

These are useful source facts, but they are not yet a single Coulomb factor
packet.

Audit answers:

1. The current parent axis bundle does not carry reviewed one-body
   electron-nuclear factors indexed by nucleus and Gaussian term. It carries
   centered one-body Gaussian factors for `axis.pgdg_intermediate.center`
   through `axis.pgdg_intermediate.gaussian_factor_terms` and
   `axis.pgdg_intermediate.gaussian_factors`. The existing White-Lindsey
   nuclear helpers construct per-center axis term tables on demand by rebuilding
   or caching `_mapped_ordinary_gausslet_1d_bundle(...; center = center_value)`.
   Those per-nucleus tables are not currently carried on
   `parent_axis_bundle_object`.
2. The current parent axis bundle does carry 1D Gaussian pair-factor
   ingredients at `axis.pgdg_intermediate.pair_factor_terms` and
   `axis.pgdg_intermediate.pair_factors`, plus raw weighted terms at
   `axis.pgdg_intermediate.pair_factor_terms_raw`. These are electron-electron
   Coulomb ingredients, but not yet a CPB-local pair-pair Coulomb record or a
   reviewed parent electron-electron factor packet.
3. The current axis bundle carries Gaussian exponents at
   `axis.pgdg_intermediate.exponents` and `axis.exponents`; it does not carry
   the full `CoulombGaussianExpansion` object or the expansion coefficients.
   Existing route/materializer code can receive an expansion and records
   `term_coefficients_source = :coulomb_expansion_coefficients` in materializer
   options, but the transient report payload does not make the expansion
   coefficients a durable parent-axis-bundle field.
4. Current reports carry nuclei, charges, and center data in structured route
   parent fields, and the transient materializer payload may carry a QW parent
   basis object with nuclei and charges. These facts are usable inputs for a
   future Coulomb packet, but they are separate from the Gaussian axis factors.
5. A future parent Coulomb factor packet should therefore be an explicit object
   tying together:
   - parent identity and axis counts;
   - Gaussian expansion identity, coefficients, and exponents;
   - electron-nuclear axis terms indexed by center/nucleus, Gaussian term, and
     axis;
   - electron-electron pair-axis terms indexed by Gaussian term and axis;
   - nucleus labels, charges, and coordinates for electron-nuclear terms;
   - factor-space, convention, index-domain, and backend provenance labels;
   - nonclaim flags for WL/PQS realization, route/global placement,
     Hamiltonian assembly, exports, and artifacts.

The next Coulomb implementation step should not infer this packet from loose
scalar fields or route-report aliases. Either build a metadata-only parent
Coulomb factor packet from these existing structured ingredients, or first
design the per-nucleus and pair-pair CPB kernel input records if the packet
shape is still unsettled.

The first metadata-only adapter now exists as:

```text
CartesianParentGaussletBases.parent_coulomb_axis_source_summary
```

It consumes an existing parent axis bundle, a `CoulombGaussianExpansion`, and
optional structured center metadata such as `center_table`,
`nuclear_charges`/`atom_locations`, parent metadata, or a QW parent basis
object. It reports compact availability and source-path facts only. It does
not store Coulomb arrays, build per-center numerical blocks, create CPB-local
Coulomb kernels, apply WL/PQS realization, place route/global matrices, assemble
Hamiltonians, add IDA/MWG/PQS semantics, or export artifacts.

The current expected source summary is intentionally partial:

- electron-electron pair-axis ingredients can be available through
  `axis.pgdg_intermediate.pair_factor_terms` or
  `axis.pgdg_intermediate.pair_factors`;
- expansion coefficients and exponents are available from the supplied
  `CoulombGaussianExpansion`;
- center metadata can be available from structured route or parent objects;
- electron-nuclear per-center axis term tables remain missing unless they are
  explicitly built or carried;
- a CPB pair-pair Coulomb source record remains missing.

### Coulomb Port Map From Current Cartesian/WL Code

The repo already has working Cartesian/White-Lindsey Coulomb machinery. The CPB
operator-layer work should port the local kernel structure, not rediscover it
or route it through simple axis-product one-body primitives.

Current electron-nuclear construction and consumption:

- `src/ordinary_qw_raw_blocks.jl`:
  - `_qwrg_diatomic_factor_term_cache(basis, centers_1d, expansion, backend)`
    builds a cache keyed by 1D nuclear center value. Each entry is
    `bundle.pgdg_intermediate.gaussian_factor_terms`, where the bundle is built
    by `_mapped_ordinary_gausslet_1d_bundle(...; exponents =
    expansion.exponents, center = center_value, backend)`.
  - `_qwrg_contracted_nuclear_axis_term_tables(axis_functions, basis,
    centers_1d, expansion, backend)` builds or reuses those per-center
    Gaussian term tables, then projects them through
    `_nested_factorized_axis_term_tables(...)`.
  - `_qwrg_fill_direct_contracted_nuclear_matrix!(destination, x_indices,
    y_indices, z_indices, amplitudes, term_coefficients, operator_terms_x,
    operator_terms_y, operator_terms_z)` is the important local contraction:
    it loops over local Cartesian states and sums over Gaussian terms inside
    the fill.
  - `_qwrg_bond_aligned_direct_contracted_nuclear_one_body_by_center(...)`
    consumes a factorized parent basis plus x/y/z axis bundles and returns one
    dense matrix per nucleus.
  - `_qwrg_fill_staged_nuclear_submatrix!(...)` fills a rectangular local
    parent submatrix from left/right state lists and per-axis Gaussian term
    tables.
  - `_qwrg_contract_staged_nuclear_block(left_coefficients,
    right_coefficients, parent_submatrix)` applies local retained coefficients
    to that submatrix.
  - `_qwrg_bond_aligned_staged_by_center_nuclear_one_body_by_center(...)` has
    staged and product-staged methods that iterate block pairs, fill local
    parent submatrices, contract them, and place them into final by-center
    matrices.
  - `_qwrg_diatomic_nuclear_one_body_by_center(basis, bundle_x, bundle_y,
    bundle_z, expansion)` is the full direct product parent version. It builds
    per-axis center caches and calls `_mapped_coulomb_expanded_symmetric_matrix`
    with `-expansion.coefficients`.
  - `_qwrg_fixed_block_one_body_matrix(fixed_block, expansion; Z)` consumes
    `fixed_block.gaussian_sum` as the already materialized positive Gaussian
    sum and applies the nuclear charge at one-body assembly time.

- `src/ordinary_qw_operator_assembly.jl`:
  - `_qwrg_bond_aligned_general_contracted_nuclear_one_body_by_center(...)`
    contracts parent by-center nuclear matrices with a carried-space
    contraction matrix.
  - `_qwrg_bond_aligned_nested_fixed_block_nuclear_one_body_by_center(...)`
    chooses staged-by-center, factorized direct, or general contraction.
  - `_qwrg_final_nuclear_one_body_by_center(...)` mixes carried fixed-block,
    Gaussian-residual, and residual-residual nuclear blocks into final-basis
    by-center matrices.
  - `assembled_one_body_hamiltonian(operators; nuclear_charges = ...)` can
    reassemble the one-body Hamiltonian from stored kinetic and by-center
    nuclear matrices when `nuclear_term_storage = :by_center`.

Current electron-electron construction and consumption:

- `src/ordinary_qw_raw_blocks.jl`:
  - `_qwrg_gausslet_interaction_matrix(data, expansion)` consumes
    `data.pair_gg_terms` on x/y/z and `expansion.coefficients` to build the
    direct-product density-density interaction matrix.
  - `_qwrg_diatomic_interaction_matrix(bundle_x, bundle_y, bundle_z,
    expansion)` consumes
    `bundle_*.pgdg_intermediate.pair_factor_terms` and
    `expansion.coefficients` to build the direct-product density-density
    interaction matrix with `_mapped_coulomb_expanded_symmetric_matrix`.
  - `_qwrg_fixed_block_interaction_matrix(fixed_block, expansion)` consumes
    `fixed_block.pair_sum`.

- `src/ordinary_qw_operator_assembly.jl`:
  - `_qwrg_interaction_matrix_nearest(...)` extends a carried gausslet or
    fixed-block interaction matrix to residuals by nearest-source copying.
  - `_qwrg_mwg_interaction_components(...)` builds matched-width Gaussian
    gausslet-residual and residual-residual interaction components from split
    pair matrices and analytic Gaussian blocks.
  - `_qwrg_interaction_matrix_mwg(...)` and
    `_qwrg_diatomic_interaction_matrix_mwg(...)` assemble those MWG components
    into final two-index density-density interaction matrices.
  - `_qwrg_bond_aligned_molecular_interaction_matrix(...)` chooses direct
    product versus fixed-block interaction, and `:ggt_nearest` versus `:mwg`.

These electron-electron routines currently produce two-index density-density
interaction matrices for the WL/QW operator payload. They are not CPB-local
pair-pair records yet, and they are not four-index Galerkin Coulomb tensors.

Gaussian expansion source:

- `src/ordinary_coulomb.jl` defines `CoulombGaussianExpansion` and
  `coulomb_gaussian_expansion(...)`.
- Frontend defaults such as
  `_normalized_nested_source_frontend_context(...; expansion =
  coulomb_gaussian_expansion(doacc = false))` and
  `bond_aligned_diatomic_nested_fixed_source(...; expansion = ...)` pass the
  expansion into axis-bundle construction and source assembly.
- Route materializer helpers also default or receive `white_lindsey_expansion`,
  and materializer options record
  `term_coefficients_source = :coulomb_expansion_coefficients`.

Per-axis transformed Gaussian factors and caches:

- `_mapped_ordinary_pgdg_intermediate_1d(...)` builds
  `gaussian_factor_terms`, `pair_factor_terms_raw`, `pair_factors`, and
  `pair_factor_terms` from `expansion.exponents`.
- `_qwrg_gausslet_1d_blocks(bundle)` exposes these as `factor_gg_terms` and
  `pair_gg_terms` for raw block assembly.
- `_qwrg_diatomic_factor_term_cache(...)` and
  `_qwrg_contracted_nuclear_axis_term_tables(...)` are the main per-center
  electron-nuclear caches.
- `_qwrg_mwg_interaction_components(...)` builds split gausslet-residual and
  residual-residual pair matrices for MWG electron-electron interaction.

Nuclei, charges, and center metadata:

- `BondAlignedDiatomicQWBasis3D` and related QW basis objects carry `nuclei`
  and `nuclear_charges`.
- `CartesianParentGaussletBasis3D` metadata preserves `nuclei`,
  `nuclear_charges`, and related family fields for QW parent objects.
- Route parent reports carry `report.nuclear_charges`,
  `report.atom_locations`, and `report.center_table`.
- `report.route_materializer_payload.parent_qw_basis_object` may carry the QW
  basis object, and `report.route_materializer_payload.parent_axis_bundle_object`
  carries the current axis bundles.

Existing tests and oracle surfaces:

- `test/nested/cartesian_nested_fixed_block_qw_pgdg_adapter_runtests.jl`
  checks that `fixed_block.gaussian_sum` and `fixed_block.pair_sum` match the
  nested shell packet, and that the nested operator carries a valid
  `interaction_matrix`.
- `test/nested/cartesian_basis_bundle_export_runtests.jl` checks Hamiltonian
  bundle export of `interaction_matrix`, `nuclear_term_storage = "by_center"`,
  default nuclear charges, and `nuclear_one_body_by_center` records.
- `test/nested/pqs_source_box_route_driver_report_runtests.jl` checks route
  configured WL Hamiltonian export statuses including
  `:available_low_order_density_density_interaction_matrix`.
- `test/nested/cartesian_white_lindsey_adapter_oracle_comparison_runtests.jl`
  is an overlap/one-body adapter oracle surface for selected safe terms. It is
  useful context but not a Coulomb CPB kernel test.

Smallest clean CPB adaptation boundary:

1. Parent/axis source adapter:
   - consume existing `_CartesianNestedAxisBundles3D` plus
     `CoulombGaussianExpansion` plus structured center metadata;
   - expose a compact parent Coulomb packet or source summary with expansion
     coefficients/exponents, per-axis pair-factor terms, and per-center
     one-body Gaussian term access;
   - do not copy route-report scalar fields into a flat cloud.
2. Electron-nuclear CPB-local kernel:
   - adapt the `_qwrg_fill_staged_nuclear_submatrix!` loop shape to CPB
     intervals, using left/right CPB local state lists and per-axis center
     term tables;
   - return one CPB-local rectangular block per nucleus/center or a compact
     by-center record;
   - keep alpha/Gaussian summation inside the kernel.
3. Electron-electron CPB-local kernel:
   - adapt the `_qwrg_diatomic_interaction_matrix` /
     `_mapped_coulomb_expanded_symmetric_matrix` factor product to CPB
     pair-pair windows;
   - return a local two-index density-density pair-pair record first, unless a
     reviewed design asks for another local two-body representation;
   - keep the existing global/Hamiltonian `interaction_matrix` paths as oracle
     surfaces, not the first CPB API.

Parent-owned versus WL/route-specific split:

- Parent-owned Coulomb packet data should include parent identity, axis counts,
  expansion coefficients/exponents, backend/provenance, per-axis
  electron-electron pair terms, and per-center electron-nuclear axis-term
  access keyed by center/nucleus and axis.
- WL realization should own carried-space contraction matrices, residual
  centers/widths, `:ggt_nearest` and `:mwg` residual policies, final operator
  payload assembly, and by-center Hamiltonian reassembly.
- Route/global layers should own artifact export, final Hamiltonian bundle
  writing, retained/global placement, and route status reporting.

First tiny equivalence test to add later:

- Build a tiny CPB/source-pair electron-nuclear block for one center using the
  same per-axis term tables as `_qwrg_fill_staged_nuclear_submatrix!`, then
  compare it to the corresponding submatrix from the existing WL
  `_qwrg_diatomic_nuclear_one_body_by_center` or staged-by-center path.
- For electron-electron, build a tiny CPB pair-pair density-density block from
  `bundle_*.pgdg_intermediate.pair_factor_terms` and compare it to the
  corresponding submatrix of `_qwrg_diatomic_interaction_matrix`.
- Do not make the first equivalence test build a full Hamiltonian or export an
  artifact.

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

## CPB Operator-Block Provider Layer

A CPB is only coordinate-window geometry. A CPB operator-block provider binds a
parent factor packet to local operator-block operations over CPB pairs:

```julia
struct CartesianCPBBlockProvider{P,F,M}
    parent::P
    parent_factors::F
    metadata::M
end
```

The preferred module name is `CartesianCPBBlockProviders`. "Provider" is more
accurate than "matrix kernel" for this layer because the functions may return
views, copied slices, factored axis blocks, or dense local operator blocks.
Dense matrix materialization is only one possible service, not the defining
behavior.

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

The generic local dense primitive for one product term is an axis-product CPB
operator block:

```text
axis_ops = (x = Ox, y = Oy, z = Oz)

O[(ix,iy,iz), (jx,jy,jz)] =
  Ox[ix,jx] * Oy[iy,jy] * Oz[iz,jz]
```

The implemented `cpb_axis_product_operator_block` materializes this dense local
CPB product-space block with x slowest and z fastest. It is generic over the
operator term and works for point, edge, face, and rectangular/cube CPB shapes
where one or more axis lengths may be one. Its summary is compact: the dense
matrix lives on the returned object, not in metadata. Realization and
route/global flags remain unset.

The sibling primitive for simple separable one-body sums is
`cpb_sum_of_axis_products_operator_block`. It accepts a compact list of product
terms, each with a scalar coefficient and prepared local axis operators, then
accumulates the dense CPB product-space sum. This is still an axis-product
family kernel; it is not a Coulomb Gaussian-sum kernel and should not be used as
one.

Overlap is the first thin wrapper over this primitive:

```text
parent overlap packet
-> CPB interval pair
-> local overlap axis blocks
-> cpb_axis_product_operator_block(term = :overlap)
```

It should not grow a separate overlap-only contraction kernel.

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

The current CPB-local kinetic wrapper follows this contract. It requires the
parent axis factor packet to carry both overlap and kinetic 1D factors, slices
those factors over the CPB interval pair, and delegates dense local
materialization to `cpb_sum_of_axis_products_operator_block` with the three
explicit terms `Kx Sy Sz`, `Sx Ky Sz`, and `Sx Sy Kz`. Position and x2 wrappers
are single axis-product terms: `position_x = Xx Sy Sz`,
`position_y = Sx Xy Sz`, `position_z = Sx Sy Xz`, `x2_x = X2x Sy Sz`,
`x2_y = Sx X2y Sz`, and `x2_z = Sx Sy X2z`. These wrappers are not production
Hamiltonian paths, not WL/PQS realization, and not route/global placement.

Future one-body terms should be sums or wrappers around the same axis-product
primitive, for example:

```text
kinetic    = Kx Sy Sz + Sx Ky Sz + Sx Sy Kz
position_x = Xx Sy Sz
x2_y       = Sx X2y Sz
```

Inactive directions should use explicit overlap factors such as `Sx`, `Sy`,
and `Sz`, not ambiguous "identity" labels. That preserves the tensor-product
operator meaning and prevents CPB-local blocks from silently changing the
one-body convention.

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
not retained blocks, not route/global matrices, not Hamiltonian objects, and
not evidence of route adoption. White-Lindsey may consume these local blocks
nearly directly. PQS realization applies local retained transforms such as
`O_retained = T_left' * O_cpb * T_right`. Retained/global range assignment and
route/global accumulation remain downstream responsibilities.

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
Even if placeholder values for all placement facts are supplied, the candidate
must remain blocked on `:placement_not_implemented` until a reviewed placement
engine exists.

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

The later synthetic retained-transform pilot should be treated as proof that
local CPB operator blocks can be transformed by realization-specific maps, not
as proof that route/global placement is ready. The next useful implementation,
if any, should be an overlap-only cleanup/renaming pass or a small generic
`CPBOperatorBlock` design. It should not be another placement metadata layer
unless a reviewed White-Lindsey/PQS realization design has first identified
real retained transforms and range ownership.

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
