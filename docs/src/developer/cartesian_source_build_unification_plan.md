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

Since this note was first written, the operator-side unification has moved
materially forward. The main remaining duplication pressure is no longer the
pure/supplement operator builders themselves. It has shifted upward to the
source-backed nested front-door layer and to the lack of one clear cross-family
nested-source contract.

The current repo direction should therefore be understood as two related
unification lines:

- unify operator orchestration through normalized internal build contexts
- unify source-backed nested front doors and diagnostics through a shared
  nested-source contract

That consolidation line is now deemed successfully accomplished for the
current repo phase. Further work in this area is now optional follow-through:

- public/front-door simplification if it becomes worth the churn
- selective post-unification cleanup when new dead residue appears
- test-suite semantic cleanup and trimming now that the runner split has given
  the new architecture clearer ownership in the regression surface

## Target architecture

The target architecture is:

1. geometry/source construction
2. normalized nested-source contract
3. normalized Cartesian operator-build contract
4. generic raw/operator pipeline
5. generic final-space mixing layer

The intended direction is:

- geometry-policy files own geometry decisions and source diagnostics
- a normalized internal nested-source contract owns source/fixed-block/
  diagnostics-facing metadata
- a normalized internal Cartesian build contract owns operator-facing metadata
- a generic operator builder consumes that contract
- a generic final-space mixer handles residual closure and carried-plus-
  residual finalization

This is the clean modern equivalent of legacy `atoms -> gethams` for the
repo.

## Proposed shared nested-source contract

Introduce one internal nested-source contract for the source-backed Cartesian
line.

Suggested semantic name:

- `CartesianNestedSourceContract3D`

This should be a common internal interface by convention, not a giant
science-erasing tagged union.

The current bond-aligned diatomic source contract is the strongest glass-box
version of this idea in the repo. That should not remain diatomic-only. It
should become the model for the general source-backed nested contract.

### Common required data

- `basis_family`
  - examples: `:atomic`, `:bond_aligned_diatomic`,
    `:bond_aligned_homonuclear_chain`,
    `:axis_aligned_homonuclear_square_lattice`
- `geometry_policy`
  - examples: `:one_center_complete_shell`, `:bond_aligned_split`,
    `:chain_recursive_split`, `:square_planar_split`
- `parent_basis`
  - canonical parent basis object
- `sequence`
  - canonical nested shell/child sequence backing the fixed block
- `fixed_dimension`
  - final fixed-space dimension implied by the source
- `contract_audit`
  - explicit support/ownership/coverage audit surface
- `shell_provenance`
  - glass-box shell origin metadata
  - examples: `source_box`, `next_inner_box`, `source_point_count`,
    `retained_fixed_count`
- `source_metadata`
  - route-specific sidecars needed for reporting/export

### Split-geometry extension

Diatomic, chain, and square should additionally expose a common split-geometry
extension where it is real:

- `shared_shell_dimensions`
- `leaf_count` or equivalent subtree count
- `node_summaries`
- `topology diagnostics`

The common contract should support these as a genuine extension, not by
pretending all families have the same topology.

### Atomic as a smaller subset

Atomic should converge to the same nested-source contract only at the smaller
common subset:

- sequence
- fixed dimension
- contract audit
- shell provenance

It should not be forced into a fake recursive split-tree abstraction just to
match diatomic or chain.

### Design constraints

- This contract should be internal.
- It should not own geometry-specific split logic.
- It should not become the operator-build contract.
- It should support glass-box diagnostics and export, not hide them.

## Proposed normalized operator-build contract

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

Status: done

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

Status: done

Scope:

- move route-capability decisions onto the normalized build contract
- replace repeated route-specific gating helpers with generic capability
  validators

Outcome:

- less policy duplication
- clearer separation between geometry decisions and backend support decisions

### Phase 3: unify pure-route operator assembly

Status: done

Scope:

- one generic pure direct-product builder
- one generic pure nested fixed-block builder

Outcome:

- direct-product and nested pure routes differ mainly in carried-space kind,
  not in route-family wrapper structure

## Longer-range reorganization plan

### Phase 4: unify supplement-bearing molecular routes

Status: done

Scope:

- direct-product diatomic supplement route
- nested diatomic supplement route

Unification target:

- common raw-block preparation
- common final-space mixing path
- route differences expressed through the normalized carried-space contract

### Phase 5: unify supplement-bearing atomic routes where real algebra matches

Status: done

Scope:

- atomic direct-product supplement route
- atomic nested fixed-block supplement route

Important constraint:

- do not force atomic and molecular supplement logic into a false generic layer
  before the residual/final-space algebra is genuinely common

### Phase 6: unify source-backed nested front-door routes

Status: done

Scope:

- public bond-aligned diatomic nested source / fixed-block / diagnostics front
  doors
- experimental bond-aligned chain nested wrappers
- experimental axis-aligned square-lattice nested wrappers

Unification target:

- common reference-only backend validation
- common axis-bundle preparation
- common nested term-coefficient resolution
- common source -> diagnostics/fixed-block packaging
- route differences expressed through source contract metadata rather than
  wrapper-local orchestration

### Phase 7: converge nested-source contracts across split geometries

Status: done

Scope:

- diatomic source-backed contract
- chain source-backed contract
- square source-backed contract
- atomic on the smaller common subset

Unification target:

- one clear internal nested-source contract for source-level diagnostics,
  shell provenance, and source-backed fixed-block construction
- diatomic contract generalized where genuinely reusable
- chain and square brought onto the same glass-box shell/audit surface
- atomic converged only at the subset that is real

### Phase 8: optional broader public API cleanup

Status: deferred / optional

Only after the internal build contract is real and stable:

- reduce public overload redundancy
- keep geometry-specific public helper names when they remain scientifically
  meaningful
- prefer thin family-specific wrappers over family-specific operator logic

## Progress so far

As of 2026-04-22:

- the normalized operator-build contract is real
- operator-side orchestration is unified for:
  - pure bond-aligned direct-product routes
  - pure bond-aligned nested fixed-block routes
  - bond-aligned diatomic supplement direct-product routes
  - bond-aligned diatomic supplement nested fixed-block routes
  - atomic supplement direct-product routes
  - atomic supplement nested fixed-block routes
- source-backed nested front-door orchestration is unified for:
  - public bond-aligned diatomic nested source / fixed-block / diagnostics
    routes
  - experimental bond-aligned homonuclear chain nested wrappers
  - experimental axis-aligned homonuclear square-lattice nested wrappers
- the shared nested-source glass-box contract is real for:
  - diatomic
  - chain
  - square
  - atomic on the smaller honest subset
- immediate post-unification compatibility residue has been removed
- duplicated route-validator scaffolding has been consolidated onto
  capability-driven helpers
- the internal consolidation plan is therefore effectively complete for now

## Current bounded chunk

There is no active required bounded chunk on this consolidation line.

The next work should be treated as a separate follow-on line:

- test-suite reorganization and trimming
- optional public/front-door cleanup
- occasional cleanup passes when newly-dead scaffolding becomes obvious

## Main risks

- overloading the normalized contract with too much science-specific data
- confusing representation metadata with build policy
- collapsing atomic and molecular supplement logic too early
- collapsing atomic and split-geometry source contracts too early
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

As of 2026-04-22:

- kernel-level decomposition is in good shape
- geometry-policy separation is mostly right
- operator-side orchestration duplication has been reduced substantially
- source-backed nested front-door duplication has also been reduced
- the shared nested-source contract now exists on the intended honest subset
- the main remaining structural cleanup opportunity has moved to the test
  surface and to optional public-surface simplification

## Update rule

When a chunk lands, update this note instead of opening a new disconnected
architecture note unless the topic has clearly split into a different
independent line.
