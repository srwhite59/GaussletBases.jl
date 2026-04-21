# Cartesian Source/Build Unification Plan

## Purpose

Record the current repo-manager judgment about how the Cartesian line in
`GaussletBases` should move toward a modern equivalent of legacy:

- `atoms -> gethams`

without regressing to one giant orchestration file.

This note is intended to be durable. It should be updated as the refactor
lands, rather than replaced by ad hoc flat notes.

See also:

- [Architecture and current direction](architecture.md)
- [Cartesian nested decomposition plan](cartesian_nested_decomposition_plan.md)

## Current judgment

The main architectural duplication pressure is not in the core Cartesian
kernel. It is in the top-level orchestration surface:

- route wrappers
- backend gating
- overload/front-door expansion

The current repo split is mostly right at the geometry-policy layer:

- one-center atomic policy
- bond-aligned diatomic split policy
- chain recursive split policy
- square-lattice planar split policy

Those should remain separate.

The current repo split is too wide above that layer. The same broad operator
pipeline is expressed repeatedly through family-specific public routes and
route-specific capability checks.

## Target architecture

The target architecture is:

1. geometry/source construction
2. normalized Cartesian build contract
3. generic raw/operator pipeline
4. generic final-space mixing layer

The intended direction is:

- geometry-policy files own geometry decisions and source diagnostics
- a normalized internal Cartesian build contract owns orchestration-facing
  metadata
- a generic operator builder consumes that contract
- a generic final-space mixer handles residual closure and carried-plus-
  residual finalization

This is the clean modern equivalent of legacy `atoms -> gethams` for the
repo.

## Proposed normalized source/build contract

Introduce one internal build-context/source object for the Cartesian line.

Suggested semantic name:

- `CartesianOperatorBuildSource3D`

It should carry orchestration-facing data, not hot-kernel implementation
detail.

### Required carried data

- `basis_family`
  - examples: `:atomic`, `:bond_aligned_diatomic`,
    `:bond_aligned_homonuclear_chain`,
    `:axis_aligned_homonuclear_square_lattice`
- `geometry_policy`
  - examples: `:direct_product_only`, `:one_center_complete_shell`,
    `:bond_aligned_split`, `:chain_recursive_split`,
    `:square_planar_split`
- `carried_space_kind`
  - `:direct_product` or `:nested_fixed_block`
- `parent_basis`
  - canonical parent basis object
- `carried_object`
  - either the basis itself or the `_NestedFixedBlock3D`
- `representation`
  - canonical carried-space representation from `basis_representation(...)`
- `nuclei`
  - normalized nuclear coordinates
- `default_nuclear_charges`
  - normalized default charges or `nothing`
- `expansion_contract`
  - exponents/coefficients contract used for bundle construction and backend
    capability checks
- `contraction`
  - `nothing` for direct-product carried space, or fixed-block contraction map
    for nested routes
- `carried_one_body_data`
  - enough carried-space data for generic one-body assembly
- `capabilities`
  - supported `gausslet_backend`
  - supported `interaction_treatment`
  - pure-route support
  - supplement-bearing support
  - localized-PGDG support
  - by-center nuclear-term storage availability
- `route_metadata`
  - stable human-readable route metadata for diagnostics and errors
- `source_metadata`
  - optional geometry/source sidecars kept for glass-box diagnostics

### Design constraints

- This contract should be internal.
- It should not own geometry-specific split logic.
- It should not own raw operator kernels.
- It should not become a second version of `basis_representation(...)`.
- It should derive from `basis_representation(...)` rather than competing with
  it.

## Ownership boundaries

### What stays geometry-specific

These files should remain separate and continue to own real policy:

- `src/cartesian_nested_atomic.jl`
  - one-center shell policy
  - legacy/full-parent working-box policy
  - atomic structure diagnostics
- `src/cartesian_nested_diatomic.jl`
  - midpoint-slab policy
  - split/no-split logic
  - adaptive shell retention
  - child-box construction
  - diatomic geometry diagnostics
- `src/cartesian_nested_experimental_geometries.jl`
  - chain recursive split policy
  - odd-chain policy
  - square-lattice planar split policy
  - experimental geometry diagnostics

These files should become thin wrappers only in the sense that they should stop
participating in route-specific operator orchestration after they have built
their geometry/source objects.

### What should be unified

The main unification target is everything above geometry construction and below
the final public wrappers.

#### `src/ordinary_qw_nested_frontends.jl`

This should shrink toward:

- public nested source/fixed-block entry wrappers
- normalization into the internal Cartesian build contract
- minimal route labeling

It should stop owning repeated backend-policy logic beyond the thinnest public
sanity checks.

#### `src/ordinary_qw_operator_assembly.jl`

This should be reorganized around three generic internal phases:

1. normalize/build context
2. build raw carried/supplement blocks
3. finalize into `OrdinaryCartesianOperators3D`

The broad overload families should collapse toward a smaller number of generic
internal builders:

- pure direct-product builder
- pure nested fixed-block builder
- supplement-bearing molecular builder
- supplement-bearing atomic builder
- shared final-space mixer

#### `src/cartesian_basis_representation.jl`

This should remain the canonical metadata/representation extractor.

It should be the source of truth for:

- carried-space kind
- parent kind
- axis sharing
- support/fixed metadata
- route metadata used by generic builders

It should not absorb the operator-build pipeline itself.

## Short-range reorganization plan

### Phase 1: introduce the normalized internal build contract

Scope:

- create the internal Cartesian build-context/source abstraction
- use it only for the already-shared bond-aligned pure routes

Targets:

- direct-product pure bond-aligned ordinary routes
- pure nested fixed-block bond-aligned routes

Outcome:

- fewer family-specific wrappers in operator assembly
- fewer repeated backend validators
- one internal seam for later supplement-route unification

### Phase 2: collapse repeated backend gating onto capability metadata

Scope:

- move route-capability decisions onto the normalized build contract
- replace repeated route-specific gating helpers with generic capability
  validators

Outcome:

- less policy duplication
- clearer separation between geometry decisions and backend support decisions

### Phase 3: unify pure-route operator assembly

Scope:

- one generic pure direct-product builder
- one generic pure nested fixed-block builder

Outcome:

- direct-product and nested pure routes differ mainly in carried-space kind,
  not in route-family wrapper structure

## Longer-range reorganization plan

### Phase 4: unify supplement-bearing molecular routes

Scope:

- direct-product diatomic supplement route
- nested diatomic supplement route

Unification target:

- common raw-block preparation
- common final-space mixing path
- route differences expressed through the normalized carried-space contract

### Phase 5: unify supplement-bearing atomic routes where real algebra matches

Scope:

- atomic direct-product supplement route
- atomic nested fixed-block supplement route

Important constraint:

- do not force atomic and molecular supplement logic into a false generic layer
  before the residual/final-space algebra is genuinely common

### Phase 6: optional broader public API cleanup

Only after the internal build contract is real and stable:

- reduce public overload redundancy
- keep geometry-specific public helper names when they remain scientifically
  meaningful
- prefer thin family-specific wrappers over family-specific operator logic

## First bounded implementation chunk

The first bounded chunk should be:

- introduce the internal normalized Cartesian build contract
- route the already-shared bond-aligned pure routes through it
- leave supplement-bearing and atomic routes unchanged

Concretely, this means:

- normalize:
  - `BondAlignedDiatomicQWBasis3D`
  - `BondAlignedHomonuclearChainQWBasis3D`
  - `AxisAlignedHomonuclearSquareLatticeQWBasis3D`
  - `_NestedFixedBlock3D{<:AbstractBondAlignedOrdinaryQWBasis3D}`
- build through:
  - one generic pure direct-product builder
  - one generic pure nested fixed-block builder
- move backend-capability checks for those routes onto contract metadata

Why this chunk first:

- it cuts real orchestration duplication immediately
- it stays away from the most science-sensitive supplement/residual closure
- it validates the new architecture on routes that already share most of their
  algebra

## Main risks

- overloading the normalized contract with too much science-specific data
- confusing representation metadata with build policy
- collapsing atomic and molecular supplement logic too early
- turning capability metadata into aspirational claims rather than current
  truth
- hiding geometry diagnostics behind a falsely generic layer

## Non-goals

This reorganization effort does not mean:

- merge distinct geometry-policy files into one generic geometry engine
- collapse back to one giant orchestration file
- erase glass-box diagnostics
- force all Cartesian science lines into one public front door immediately

## Current status

As of 2026-04-20:

- kernel-level decomposition is in good shape
- geometry-policy separation is mostly right
- orchestration duplication above that layer remains the main structural bloat
- the first implementation target is the pure bond-aligned ordinary route
  family

## Update rule

When a chunk lands, update this note instead of opening a new disconnected
architecture note unless the topic has clearly split into a different
independent line.
