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
- [Cartesian QW receipt wrapper status](cartesian_qw_receipt_wrapper_status.md)
- [Cartesian nested decomposition plan](cartesian_nested_decomposition_plan.md)
- [Projected q-Shell policy](projected_q_shell_policy.md)

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

As of 2026-05-26, the QW receipt-wrapper line has a clear internal coverage
contract. The covered routes share the same pattern:

- record the request as `CartesianOperatorBuildSource3D`
- run the existing authoritative `ordinary_cartesian_qiu_white_operators(...)`
  builder
- audit the result with `CartesianQWOperatorConstructionRecord3D`
- return a `CartesianQWOperatorConstructionReceipt3D`

This is a wrapper/audit layer, not a new Hamiltonian builder. See
[Cartesian QW receipt wrapper status](cartesian_qw_receipt_wrapper_status.md)
for the covered and intentionally uncovered route families.

As of 2026-05-27, the default q4 high-order atom-growth/endcap-panel recipe
has an internal opt-in construction path, but it is not a default route and is
not CR2-validated. The implementation trail is:

- `3470566`: added explicit opt-in recipe source construction for the ready
  q4 atom-growth/endcap-panel policy
- `faa5e60`: added readiness diagnostics and fixed-block handoff auditing
- `c15f8e8`: added the first PGDG QW construction smoke through the existing
  nested fixed-block builder

The ready default recipe can now be selected, realized, built into a shell
sequence, and converted into `_NestedFixedBlock3D` through:

- `_nested_bond_aligned_diatomic_high_order_recipe_source_construction(...)`
- `_nested_bond_aligned_diatomic_high_order_recipe_source_readiness(...)`
- `_nested_bond_aligned_diatomic_high_order_recipe_source_fixed_block(...)`

The current small q4 smoke fixture records:

- fixed block size: `(735, 469)`
- operator dimension: `469`
- `gausslet_backend = :pgdg_localized_experimental`
- `interaction_treatment = :ggt_nearest`
- residual count: `0`
- clean `CartesianQWOperatorConstructionReceipt3D` source/sidecar agreement
- staged by-center sidecar preserved on the fixed block
- finite symmetric overlap, one-body, and interaction matrices

The first dimensionally matched capture/H1 sanity check is also repo-native
only. It generated eight parent one-body generalized eigenvectors on the exact
same q4 smoke parent and then projected them into the opt-in fixed block. The
target is labeled `repo_native_parent_one_body_generalized_eigenvectors`; it is
not PySCF data, not CR2 evidence, and not an energy validation. The artifact
paths used for this note were:

- `tmp/work/q4_recipe_parent_reference_capture_report.txt`
- `tmp/work/q4_recipe_parent_reference_capture.tsv`

That parent-reference target records:

- parent dimension: `735`
- fixed dimension: `469`
- parent axis counts: `(7, 7, 15)`
- coefficient ordering:
  `cartesian_flat_index = (ix - 1) * ny * nz + (iy - 1) * nz + iz`
- parent and fixed backends: `:pgdg_localized_experimental`
- warning-level fallback logs rejected; no `:numerical_reference` fallback
- target vector count: `8`
- primary lowest-vector capture: `0.9999809353`
- primary H1 delta: `4.69e-5`
- worst capture over the eight vectors: `0.9876930450`
- max absolute H1 delta over the eight vectors: `0.0140`

This is useful evidence that the q4 fixed block can capture a compatible
repo-native parent one-body target. It does not replace a real target-space
occupied capture, PySCF comparison, CR2 handoff, same-density test, or chemistry
validation.

Pass 001 of the 2026-05-27 CR2 handoff added the first external-target
capture/H1 diagnostic for this exact q4 fixture. The target provenance was
duplicated one-center PySCF Cr ccECP occupied spaces, source-metric Lowdin in
the two-center GTO AO metric. The handoff matrices had shapes alpha
`(735, 20)`, beta `(735, 8)`, and combined `(735, 28)`. The repo diagnostic
verified the q4 fixture, grid order, and center coordinates exactly, with
parent and fixed backends both `:pgdg_localized_experimental`, warning log
count `0`, and `numerical_reference_fallback = false`.

That external target had lower capture into the q4 parent itself than into the
q4 fixed block: CR2-reported parent capture was alpha `0.9558876918` and beta
`0.9811637481`. Conditional on that q4 parent target, the opt-in fixed block
captured alpha `0.9994161870`, beta `0.9998427016`, and combined spin-sum
`0.9995403327`. The worst alpha column was `alpha_right_mo7_col18` with
capture `0.9987695368`; the worst beta column was `beta_left_mo0_col1` with
capture `0.9998415239`. The max absolute H1 delta was `7.249921e-03`.

This separates two questions. The q4 fixed block captures the supplied q4
parent-space target well, but the q4 parent target already loses some of the
original CR2 occupied-space norm. This is still only a fixed-space capture/H1
diagnostic. It is not relaxed Cr2 energy evidence, not a two-electron or
same-density validation, not production default support, and not public route
readiness. The q label remains a convergence/control parameter, not a validated
accuracy tier.

The next validation stage should therefore be a q-family/q-ladder study, not a
deeper one-off q4 conclusion. q4 is the first integration and acceptance
fixture because it exercises the opt-in recipe, fixed-block handoff, PGDG QW
smoke, receipt diagnostics, and external-target capture/H1 path end to end.
The useful q range is a scientific and numerical convergence question that
depends on molecule, geometry, ECP, basis, and target accuracy. Future q rows
must keep the same provenance gates: support ownership and grid/order audits,
PGDG/no-fallback behavior, carried-space and receipt diagnostics, and
per-column/per-spin capture plus H1 reporting. Reports should separate
parent-target loss from fixed-space loss, and q-specific evidence should be
treated as validation anchors rather than hard-coded architecture.

The current implementation is not yet a full q-family construction interface.
Its per-region q status is:

- atom-local cubic boxes: q/order metadata exists, but construction is still
  driven by box/support size and the existing complete-shell sequence path
  rather than by a q-retention knob
- shared endcap/panel exterior: q/L is consumed by the owned-unit product
  `doside` construction, so this is the first natural parameterization seam
- contact cap: q/order metadata exists, but construction currently uses direct
  full-support slab coefficients
- outer mismatch shells: q/order metadata exists, but construction currently
  uses direct full-support boundary slabs
- transverse annulus: optional experimental policy metadata only; there is no
  owned-unit producer yet

The smallest future implementation step should therefore be limited to shared
endcap/panel q/L parameterization unless a manager explicitly scopes broader
q-family source construction. Any CR2 q-ladder report should state which
regions actually consumed q, and should keep parent-target loss separate from
fixed-space loss.

Commits `187e327`, `a58f1c0`, and `34fc4b7` implement and smoke-test that
first shared endcap/panel q/L step. The internal opt-in path now accepts
variable shared-region q/L while keeping atom-local, contact-cap, and
outer-mismatch regions fixed at q4/order4 direct-support semantics. The
current small-fixture shared endcap/panel rows are:

- q4/order4: fixed dimension `469`, retained counts
  `[98, 96, 125, 125, 25]`
- q5/order5: fixed dimension `523`, retained counts
  `[98, 150, 125, 125, 25]`
- q6/order6: fixed dimension `589`, retained counts
  `[98, 216, 125, 125, 25]`

The q5 and q6 rows both have focused PGDG QW smoke coverage through the
existing receipt path: backend `:pgdg_localized_experimental`, zero residuals,
clean source/sidecar agreement, no dense parent matrix, no heavy metric
packet, no new Hamiltonian kernel, and warning-level fallback logs rejected.
q7/order7 is not supported for this fixture because the requested nested
`doside` retained count exceeds the interval size.

These rows are construction-smoke evidence for the shared endcap/panel seam
only. They are not Cr2 accuracy evidence, energy validation, public q-ladder
API support, or a change to atom-local/contact/mismatch semantics. The
transverse annulus remains experimental/missing.

Commit `71d7084` adds a private q-row route receipt wrapper for that smoke
path. `_nested_bond_aligned_diatomic_high_order_q_row_route_receipt(...)` lives
in the QW receipt/diagnostic layer and is not exported. It accepts an
already-built bond-aligned diatomic basis plus explicit shared q/order, then
composes the existing chain:

```text
basis -> q policy -> source construction -> readiness/fixed block ->
QW receipt -> route diagnostics
```

The wrapper is route/provenance infrastructure only. It keeps non-shared
regions fixed at q4/order4 direct-support semantics, requires
`gausslet_backend = :pgdg_localized_experimental`, rejects
`:numerical_reference`, and preserves the q7/order7 retained-count rejection.
Its focused tests pin q4/q5/q6 success and q7 failure through the existing
PGDG receipt path. This is still not a public q-ladder API, not a GTO
supplement path, not Be2/Cr2 science validation, not an MWG/HF route, not a
new Hamiltonian kernel, and not a quadrature-route change.

Commit `d36b81d` adds the next private fixture layer:
`_nested_bond_aligned_homonuclear_high_order_q_row_fixture_receipt(...)`.
This helper is still internal and unexported. It constructs the existing
`bond_aligned_homonuclear_qw_basis(...)` from explicit small-molecule geometry
and then delegates unchanged to the q-row route receipt above. The required
fixture inputs are `bond_length`, `core_spacing`, `xmax_parallel`,
`xmax_transverse`, and `shared_q`. Its explicit default provenance includes
`family = :G10`, `bond_axis = :z`, `nuclear_charge = 1.0`,
`reference_spacing = 1.0`, `tail_spacing = 10.0`,
`shared_order = shared_q`, `protected_atom_side_count = 5`, `q_min = 4`,
`nside = 5`, and `gausslet_backend = :pgdg_localized_experimental`.

The fixture forwards `basis.nuclear_charges` from the constructed homonuclear
basis and deliberately does not accept a separate `nuclear_charges` override.
It records the basis constructor inputs, nuclei, parent axis counts, parent
dimension, Cartesian flat-index convention, and delegated q-row diagnostics so
a later CR2 handoff can identify the parent grid and ordering. This is only a
parameterized construction/provenance bridge. It is not heteronuclear support,
not element-label provenance, not a supplement route, not Be2/Cr2 science
validation, not an energy result, and not a change to source builders,
Hamiltonian kernels, backend defaults, or quadrature policy.

Commit `de26059` adds a private supplement-aware variant:
`_nested_bond_aligned_homonuclear_high_order_q_row_fixture_supplement_receipt(...)`.
It composes the homonuclear q-row fixture wrapper with the existing nested
molecular supplement receipt, constructs a
`LegacyBondAlignedDiatomicGaussianSupplement` through
`legacy_bond_aligned_diatomic_gaussian_supplement(atom, basis_name, basis.nuclei; ...)`,
and forwards `fixed_block.parent_basis.nuclear_charges` into the delegated
operator construction. The wrapper requires
`gausslet_backend = :pgdg_localized_experimental` and rejects
`:numerical_reference`.

The focused q4 H/cc-pVTZ supplement test uses `lmax = 0` and
`max_width = 1.0`, then compares exact equality against the existing direct
builder for matrices, counts, raw-to-final, residual metadata,
backend/treatment/storage, kinetic sidecars, and by-center nuclear sidecars.
The focused path rejects warning-level logs so numerical fallback warnings
cannot pass silently there. This remains a private delegate/audit bridge only:
not public API, not heteronuclear support, not a default route, not a new
Hamiltonian path, not a new supplement representation, and not
Be2/Cr2/PySCF/HF/energy validation.

Commit `a525b80` adds the first private capture/H1 diagnostic for the q-row
plus supplement route:
`_nested_bond_aligned_homonuclear_high_order_q_row_fixture_supplement_capture_h1(...)`.
The helper accepts exact parent-grid target coefficient matrices only. The row
count must equal the fixture parent dimension; for the focused q4 fixture this
means `(7, 7, 15)` parent axes, `735` rows, and the repo Cartesian flat order
with `z` fastest. The helper rejects wrong row counts, non-finite target
entries, bad label or occupation vector lengths, negative occupations, and
`:numerical_reference`.

The diagnostic builds the existing q-row supplement receipt and reads back
existing QW operator/sidecar data. It reports fixed-only capture, final hybrid
capture, parent/fixed/final H1 expectations, final-overlap trust-gate
diagnostics, route diagnostics, and the first-pass omissions that are still
intentional. Final self-overlap remains diagnostic only; this pass did not add
generalized final-overlap transfer logic. The focused regression uses
repo-native q4 parent-grid columns as targets, not CR2 or PySCF artifacts.
This is capture/H1 diagnostic infrastructure only. It is not Be2/Cr2 science,
energy validation, HF/ED, same-density validation, a public API, or default
route adoption.

The first CR2 external-target run on 2026-05-27 exercised this plumbing for
the q4 parent-grid fixture. CR2 supplied parent-grid alpha/beta target
matrices matching `(7, 7, 15)` / `735` rows, and the repo driver accepted the
artifact without falling back to `:numerical_reference`. The run used the
private H/cc-pVTZ supplement fixture with `lmax = 0`, `max_width = 1.0`,
`shared_q = shared_order = 4`, and backend
`:pgdg_localized_experimental`. The reported dimensions were fixed dimension
`469`, final dimension `471`, and residual count `2`.

The external-target capture/H1 checkpoint values were:

- alpha final capture `0.9994163575`
- beta final capture `0.9998428679`
- alpha max final H1 delta `7.25` mHa
- beta max final H1 delta `2.95` mHa

This result means the external parent-grid target path is working: CR2 can
stage a target, the repo can verify the grid/order contract, the PGDG route is
enforced, and fixed/final capture plus H1 readback can be reported. It is not
Cr ECP physics validation, not Cr2 science validation, not energy validation,
and not evidence for two-electron or same-density accuracy. The H/cc-pVTZ
`lmax = 0` supplement fixture is deliberately a plumbing target, not a
chromium ECP supplement model.

The follow-on q-ladder CR2 run on 2026-05-27 used the same `735`-row
parent-grid target for q4, q5, and q6. PGDG was enforced for all three rows,
with no `:numerical_reference` fallback. The fixed/final/residual dimensions
and final capture values were:

| q | fixed/final/residual | alpha final capture | beta final capture |
|---:|---|---:|---:|
| 4 | `469 / 471 / 2` | `0.9994163575` | `0.9998428679` |
| 5 | `523 / 525 / 2` | `0.9997040846` | `0.9998748198` |
| 6 | `589 / 591 / 2` | `0.9998268146` | `0.9999020226` |

The maximum final H1 deltas decreased across the ladder:

| q | alpha max final H1 delta | beta max final H1 delta |
|---:|---:|---:|
| 4 | `7.250` mHa | `2.949` mHa |
| 5 | `5.271` mHa | `2.649` mHa |
| 6 | `4.782` mHa | `2.073` mHa |

The worst alpha orbital remained source MO 7 on the left/right labels. The
worst beta orbital shifted from source MO 0 at q4 to source MO 3 at q5 and
q6. q7 was not run and remains expected unsupported for this fixture. This
q-ladder result shows that the q-row plumbing works and that increasing q
improves capture/H1 for this target. It is still not Cr ECP energy or physics
validation: the supplement fixture remains H/cc-pVTZ with `lmax = 0`,
`max_width = 1.0`, so the result is a target-plumbing and fixed-space
diagnostic only. The cleaned local driver now writes timing reports under each
q output directory instead of emitting live timing spam to stdout.

A later Be2 total-basis closure probe corrected the target contract for
orbitals that live in the GTO supplement space. The earlier parent-grid-only
Be2 capture/H1 recipe is still useful as a parent-filtered target diagnostic,
but it is not the right closure test for a PySCF/GTO target whose closure
requires the supplement rows. The corrected total raw projection uses:

- fixed rows `<fixed(parent)|PySCF GTO> C`
- supplement rows `<supplement GTO|PySCF GTO> C`

With Be2/cc-pVDZ targets and a matching cc-pVDZ supplement at `lmax = 2`, the
total-basis route closes essentially exactly:

- q4 final capture `0.999999998974`, max final H1 error `0.001920` mHa
- q6 final capture `0.999999999397`, max final H1 error `0.001840` mHa

This supersedes the earlier interpretation that the roughly `32` mHa gap was
necessarily a parent/core-spacing limitation of the total route. Same-basis
GTO supplement closure works, and q is not limiting this closure. This is
still closure plumbing, not Be2 energy validation, not HF/same-density or
two-electron validation, and not a public/default route claim.

This status means the path is construction-smoke-ready only. It remains
explicit/internal, and active/default source builders still do not consume the
recipe policy. Legacy source-object wrapping is also not claimed: the
readiness audit records missing legacy split/source fields such as
`split_geometry`, child retention contracts, child sequences, child column
ranges, and midpoint-slab metadata. The safe handoff surface today is the
fixed block, not a replacement `_CartesianNestedBondAlignedDiatomicSource3D`.

The q5 transverse annulus recipe remains experimental/incomplete. It is
allowed to appear in planning metadata as a selected future policy, but it is
not consumed by the opt-in builder and should fail rather than silently
falling back to q4 or complete-rectangular behavior.

The longer-term high-order shell target is now Projected q-Shell (PQS), defined
as boundary COMX-product modes from the full local block transform, projected
to raw boundary rows, then cleaned with full-rank symmetric Lowdin. Cubic
atom-local shells should move toward `PQS(q, q)`, and rectangular
molecular/exterior shells should move toward `PQS(q, L)`. Existing q-row and
endcap/panel work is transitional infrastructure and validation scaffolding,
not the preferred regular shell abstraction. See
[`projected_q_shell_policy.md`](projected_q_shell_policy.md) for the policy
definition and the first q=5 C2 `R = 4.7` external preflight evidence.

As of commits `06483f8`, `bd66e51`, `e4fbcca`, and `42103e3`, mainline also
has a private opt-in PQS construction smoke. The local layer helper uses
boundary COMX-product mode selection, raw-boundary projection, and full-rank
symmetric Lowdin cleanup. PQS realization metadata is visible as a candidate,
and the opt-in source/fixed-block path can build the first small fixture with
full-parent dimension `735`.

That smoke is deliberately narrow. Defaults still use the existing
endcap/panel route where they did before, and the PQS path is private and
explicit opt-in. The pure nested QW receipt smoke uses PGDG with
`:ggt_nearest` and `:total_only`, has clean source/sidecar agreement, finite
symmetric overlap/one-body/interaction matrices, zero residuals, and rejects
warning-level numerical-quadrature logs for that path. This is construction
and QW-smoke evidence only; it is not a compression-quality, CR2, energy, or
production-default claim.

The next PQS design problem is product-staged sidecar and performance support.
PQS should not be treated as ready for by-center, supplement, or
performance-sensitive adoption until that sidecar/performance contract is
implemented and validated.

Numerical quadrature is forbidden on this PGDG smoke path. A valid smoke must
force `gausslet_backend = :pgdg_localized_experimental` and must not accept a
silent `:numerical_reference` fallback or warning-level numerical-quadrature
route. Existing numerical-reference/debug paths elsewhere in the test suite do
not change this contract.

Before CR2 use, the next scientific validation must be explicit and separate:
at minimum occupied-capture/H1 checks on the target parent and then route-
appropriate same-density or chemistry comparisons. The current smoke does not
validate Cr2 energies, CR2 workflows, public frontend support, or production
default behavior.

As of commit `928ee87`, ordinary Cartesian QW final packaging is unified for
the current builder routes. The implementation trail is:

- `16d9e9d`: introduced `_qwrg_finalize_ordinary_cartesian_operators(...)`
- `9479210`: migrated pure nested fixed-block packaging
- `7aa93df`: migrated atomic supplement packaging
- `a3d7929`: added molecular-owner, MWG-width, and by-center sidecar
  guardrails
- `928ee87`: migrated molecular supplement packaging

The helper is now used by:

- pure bond-aligned direct-product routes
- pure bond-aligned nested fixed-block routes
- one-center atomic direct-product routes with
  `LegacyAtomicGaussianSupplement`
- one-center atomic nested fixed-block routes with
  `LegacyAtomicGaussianSupplement`
- bond-aligned molecular direct-product routes with legacy Gaussian
  supplements
- bond-aligned molecular nested fixed-block routes with legacy Gaussian
  supplements

This is a final-packaging checkpoint, not a new Hamiltonian builder. The helper
is allowed to:

- validate cheap final dimensions, counts, and metadata consistency
- normalize final matrices and vectors into the stored
  `OrdinaryCartesianOperators3D` representation
- enforce explicit owner indices for multi-center residuals
- require finite positive residual widths for `interaction_treatment = :mwg`
- require kinetic and by-center nuclear sidecars when
  `nuclear_term_storage = :by_center`
- construct the final `OrdinaryCartesianOperators3D` payload from
  already-assembled inputs

The helper must not:

- build or alter overlap, one-body, interaction, residual, or raw blocks
- choose or reinterpret `gausslet_backend`
- change PGDG / numerical-reference routing
- change geometry policy, residual filtering, residual ownership, or MWG
  center/width extraction
- change timing labels or public return shapes
- become a source-driven Hamiltonian implementation

Assembly remains route-specific. Any next code seam before this final endpoint
requires a separate design pass; it should not be treated as an automatic
extension of the packaging helper migration.

### Pre-packaging seam audit

A follow-on read-only audit reviewed the likely pre-packaging extraction seams:

- one-body final mixing
- interaction symmetrization
- residual-space overlap and `raw_to_final` setup
- by-center nuclear sidecar preparation
- residual orbital metadata

No broad pre-packaging extraction is currently recommended. The remaining code
is mostly route-specific orchestration around backend policy, residual
ownership, by-center sidecars, MWG center/width semantics, and geometry
contracts. Some low-level helpers are already shared, but treating those shared
pieces as permission to merge route assembly would create conceptual drift.

Future work should start with a design and test harness for exact
pre-packaging equality fields before moving code. At minimum, that harness
would need to compare final overlap, one-body, interaction, `raw_to_final`,
residual centers, residual widths, residual owner indices, stored kinetic
sidecars, by-center nuclear sidecars, counts, and route diagnostics for each
route under consideration.

### Contracted-parent rule-seam checkpoint

The first small rule-driven seam now exists at the contracted-parent metadata
level. `CartesianContractionRule3D` can describe three current families:

- product-owned staged units from the existing endcap/panel sidecar path
- support-dense fallback units that still contract through support-local
  coefficient entries
- PQS descriptor/prototype rules, which describe the projected-q-shell idea
  but are not installed as fixed-block sidecars

For existing product-staged sidecars, contracted-parent unit creation now goes
through a rule-aware helper. That preserves the staged payload, role, support
indices, column range, and global coefficient map. The global
`CartesianContractedParent3D` coefficient matrix remains the source of truth;
rules are provenance and classification, not regenerated coefficients.

`CartesianContractionRuleInventory3D` summarizes the rules that make up a
contracted parent: rule-family counts, retained/source dimensions, parent
support coverage, metric capabilities, and whether each unit has derivable rule
metadata. Private metric dispatch shadow plans can also compare rule-derived
classification against the current staged-payload dispatch for the
product/support-dense metric packet. On the existing endcap/panel fixture, the
rule plan and payload plan agree on product-unit counts, support-fallback
counts, product/product blocks, fallback blocks, and unsupported/prototype
counts.

This is not yet a full "define a new nesting rule and everything works" layer.
Metric execution remains payload-driven, and current product-staged metric
contraction still reads concrete staged sidecar data such as axis intervals,
fixed/active axis state, 1D coefficient matrices, support states, and
support-dense block coefficients. PQS remains descriptor/prototype metadata
only and is intentionally rejected by the product-staged metric dispatch plan.

Current recommendation: keep the shadow dispatch plan as a private test/audit
helper. Do not attach full plans to packet diagnostics by default, do not add a
runtime assertion keyword, and do not route metric execution from rules until a
concrete consumer needs that guard and the payload-resolution contract is
designed. This checkpoint addresses the earlier framework concern by making
construction pieces classifiable and auditable through one rule vocabulary,
while preserving all coefficient maps, metric kernels, QW/Hamiltonian behavior,
public/default behavior, and PGDG/numerical-reference policy.

### Resolved-payload seam checkpoint

The next execution-facing seam is now in place for the existing
product/support-dense contracted-parent metric path. `CartesianShellRegion3D`
remains metadata-only: it records finite region facts and retention intent, but
it does not carry heavy coefficient blocks, full supports by default, fixed
blocks, QW operators, Hamiltonians, or metric packets.

The active execution seam for current product-staged contracted-parent metrics
is now:

```text
CartesianContractionRule3D -> resolved payload -> existing metric kernel
```

Resolved payloads normalize only executable payloads that already exist:

- product-owned staged units with `kind = :product_doside`
- support-dense staged units with `kind = :support_dense`
- `CartesianContractionUnit3D` values that carry an existing
  `metadata.staged_by_center_unit`

The public `:product_staged_metric_contraction` path is routed internally
through these resolved payloads. Its output shape and diagnostics contract are
preserved: it still reports `construction_path =
:product_staged_metric_contraction`, and it does not expose shadow-only fields
such as `resolved_payload_count` in the public packet diagnostics.

The numerical kernels remain the same authority:

- product/product blocks use the existing product-staged block fill
- mixed or support-dense blocks use the existing support-local fallback
- support-local default metric execution remains direct and unchanged

PQS remains outside production metric execution. Descriptor-only PQS is still
prototype-only metadata: it reports missing installed sidecar/product-staged
payload fields, is not installed into fixed-block sidecars, and is not
consumed by public contracted-parent metric execution. A later private
sidecar-fixture checkpoint added an explicit-`column_range` executable PQS
resolved payload for low-order reference checks only. That fixture validates
PQS self blocks and PQS/support-dense mixed blocks through private
support-local/reference contracts. PQS self-overlap may use identity only as
the post-cleanup self-block invariant; mixed overlaps do not use that shortcut.
PQS/product optimized metric blocks remain explicitly unsupported, and any
PQS/product support-local reference must be a separately named explicit
reference/debug helper rather than a silent fallback.

This PQS sidecar-fixture checkpoint does not imply kinetic, nuclear, Gaussian
or pair, interaction, QW/Hamiltonian, H1, energy, CR2, or production metric
readiness.

A follow-up private checkpoint installs that same sidecar fixture into a
synthetic single-PQS-layer `_NestedFixedBlock3D`, discoverable only through a
PQS-specific private accessor. The sidecar is stored in a fixed-block sidecar
slot for attachment/discovery testing, but ordinary by-center, QW, and
Hamiltonian consumers do not consume it; `_nested_staged_by_center_sidecar`
remains incompatible/loud for this PQS fixture. The validated checks remain
low-order and reference-scoped: PQS self blocks and PQS/support-dense mixed
blocks. The mixed q4 recipe fixed block is not covered by this fixture, and
PQS/product optimized metric blocks remain explicitly unsupported.

This checkpoint changes no QW/Hamiltonian construction, backend defaults, PGDG
policy, quadrature policy, public/default route behavior, source builders, or
coefficient maps. Any next step should be a separate scoped decision rather
than automatic continuation: either payload resolution for another existing
metric path, or a dedicated PQS sidecar design pass.

### Raw product source and retained transform direction

The next construction abstraction should treat most retained Cartesian pieces
as a raw product source space plus a retained-space transform. Operator blocks
can then be assembled in the smaller product source spaces and finished by
small dense transforms, rather than creating a new special case for every
mixed pair such as PQS/product.

This direction is documented in
[Raw product source and retained transform policy](raw_product_source_retained_transform_policy.md).
That note also records the quadrature-weight contract: gausslet-derived slabs,
shells, and PQS units start from raw product sources with positive
source-point weights in non-pathological cases. Retained-column weights need an
explicit role and finite-positive check before any IDA-style division, while
angular GTO supplements and other non-quadrature final functions must not be
treated as positive-weight quadrature carriers.

## Current bounded chunk

There is no active required bounded chunk on this consolidation line.

The next work should be treated as a separate follow-on line:

- test-suite reorganization and trimming
- optional public/front-door cleanup
- selective receipt-wrapper consumer migration where source/record/provenance
  diagnostics are useful
- occasional cleanup passes when newly-dead scaffolding becomes obvious
- a separately designed pre-packaging Hamiltonian assembly seam, if one is
  clearly identified and validated

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
