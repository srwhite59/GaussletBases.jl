# Common Terminal Shell Decomposition

Status: implemented common shellification, angular-z-extension geometry,
neutral face-product realization, and compact thin-slab lowering. The registry
owns exact ID lifecycle and source permission. The implemented source IDs are
`HP-COMP-SHELLGEOM-FN-01`, `HP-COMP-SHELLGEOM-DIAT-FN-01`,
`HP-COMP-ANGBOX-FN-01`, `HP-COMP-FACEPROD-FN-01`,
`HP-COMP-THINSLAB-FN-01`, and `HP-COMP-THINSLAB-META-FN-01`; their test IDs
are completed validation evidence with no continuing permission.
`HP-COMP-OUTERMM-*` is superseded with no source permission and must not be
restored as a separate path.

This page owns the current geometry and compact-slab contract. It does not own
PQS complete-shell source dimensions, White-Lindsey retained products, public
input defaults, or Hamiltonian semantics.

## Common Geometry Boundary

PQS and White-Lindsey share one route-family-free first operation:

```text
parent axes
+ nuclear centers
+ public ns
+ direct nucleus-centered core side
+ bond axis when present
-> direct core and contact regions
-> terminal shell and slab regions
-> deterministic owned support rows
```

For the same normalized system and parent facts, both families must enter the
common shellifier with the same parent axes, nuclear centers, public `ns`,
direct-core side, and bond axis. For z-axis diatomics, central-gap/contact,
shared-shell, angular-extension, and outer-mismatch ownership are common
geometry facts. Construction-family code must not recompute them or
reinterpret their owned rows.

Public `ns` is the common shell-size input. The direct nucleus-centered core
side is

```text
direct_core_side = isodd(ns) ? ns : ns + 1
```

and applies only to true direct core blocks. Route-local sizes are downstream
retained-construction facts: PQS uses `q = ns`, while White-Lindsey uses
`q = ns - 2`. Neither route-local value may govern common core or shell
ownership.

Common geometry must preserve:

- deterministic region and owned-support ordering;
- disjoint owned supports and complete intended terminal coverage;
- direct-core centering and atom-contact ownership;
- shell outer boxes and inner exclusions;
- native region roles, shell indices, and slab geometry;
- the distinction between real shells, direct/core sectors, and thin slabs.

Changing central-gap/contact policy, support ownership, or region ordering is
outside this contract.

## Family Boundary

Family-specific retained construction begins only after common regions exist.

PQS consumes a common complete shell as a full local source box, selects
boundary COMX/product modes, restricts them to shell-owned rows, and applies
the shell-local Lowdin correction. Its aspect-aware `(q,q,L)` source policy is
owned by `pqs_complete_shell_aspect_source_modes.md`.

White-Lindsey splits a common complete-shell boundary into facets, edges,
corners, or equivalent strata and realizes compact products of one-dimensional
contractions on each authoritative owned support. Its retained realization is
owned by `white_lindsey_terminal_basis_realization.md`.

These are different retained-construction geometries, not different first-step
shellifiers. This common contract does not choose `L`, any source-mode shape,
complete-shell retained counts, or a PQS/WL convergence policy.

## Projected-Q-Shell Staged Descriptor Retirement

`HP-RETIRE-PQS-STAGED-DESCRIPTOR-FN-01` records the completed source-only
retirement of the inert metadata/prototype descriptor attached to the otherwise
active projected-q-shell layer. `HP-RETIRE-PQS-STAGED-DESCRIPTOR-TEST-01`
records completed validation with unchanged existing owners and transient
before/after probes. Both execution grants are closed.

The descriptor was a proposed future staged sidecar. It is not consumed by the
metric packet, final basis, Hamiltonian, current product-staged sidecars, or a
public API. Its retained copies of support, source-axis coefficients, cleanup
data, prototype contractions, and false consumption flags duplicate facts from
the actual layer construction without participating in that construction.

The completed retirement deleted this private closure:

```text
_CartesianNestedProjectedQShellStagedUnitDescriptor3D
_nested_projected_q_shell_make_staged_unit_descriptor
_nested_projected_q_shell_staged_unit_descriptor
_nested_projected_q_shell_descriptor_seed_coefficients
_nested_projected_q_shell_descriptor_metric_prototype (both methods)
_nested_projected_q_shell_boundary_rectangular_pieces
_nested_projected_q_shell_piece_support_states
_nested_projected_q_shell_boundary_piece_coverage
_nested_projected_q_shell_axis_piece_pair_blocks
_nested_projected_q_shell_axis_piece_weight_vectors
_nested_projected_q_shell_boundary_mode_axis_indices
_nested_projected_q_shell_mode_matrix_from_axis_piece_blocks
_nested_projected_q_shell_cleaned_mode_matrix
_nested_projected_q_shell_descriptor_metric_product_contraction (both methods)
```

In `_nested_projected_q_shell_layer`, coefficient assembly now returns the
existing `_nested_projected_q_shell_parent_coefficients(...)` result directly.
Descriptor construction and only the diagnostics/provenance entries whose
subject was descriptor availability, descriptor metadata-only status,
prototype consumption, or the false `active_builder_consumes` claim were
removed. The two-line `pqs_staged_unit_descriptor` propagation was removed
from `src/cartesian_nested_diatomic.jl`. No deleted fact was replaced by a
smaller carrier, `NamedTuple`, accessor, alias, deprecation, serialization
record, status vocabulary, or copied diagnostic.

The active projected-q-shell construction remains authoritative and unchanged:

- boundary COMX/product-mode selection and raw-boundary projection;
- full-rank symmetric Lowdin cleanup and its deterministic gauge;
- the coefficient matrix and support indices/states;
- the metric packet and actual layer;
- numerical construction diagnostics and provenance that describe objects
  consumed by the current producer.

The generic by-center sidecar, old XY/shell-plus-core/hierarchical milestones,
active product-staged sidecars, factorized/support-reference kernels, and the
high-order branch are outside this transaction.

### Caller Proof And Accepted Boundary

At baseline `1e8e31377efb5778df8e48c9d06c4e2237fd8a81`, tracked live references
to the descriptor closure occur only in
`src/cartesian_nested_faces.jl` and the two-line propagation in
`src/cartesian_nested_diatomic.jl`. No tracked test, driver, tool, example, or
current nonhistorical contract calls it. A completed REQ-094 worker log records
an old frozen probe that used the accessor, but its retained executable replay
does not; this is historical evidence rather than a current downstream
consumer. Archived design text may retain historical names.

The implementation repeated the tracked and known paper-workspace caller scan
without finding a current non-ignored production, test, executable paper, or
downstream consumer. It required no replacement carrier or test edit and
preserved every projected-q-shell coefficient and metric-packet value.

### Accepted Budget And Validation

The implementation is confined to:

```text
src/cartesian_nested_faces.jl
src/cartesian_nested_diatomic.jl
```

Implementation commit `8fcda0086f73dbd1348aa6b261c7862f1cf64bb3`
deleted `593` and added `2` source lines, for a net reduction of `591`. It
added no file, test, helper, API, type, cache, dependency, metadata replacement,
workflow, version, or release change.

Projected-q-shell coefficient and metric-packet values, packet/basis
fingerprints, and default ordinary-QW H2 geometry, coefficients, overlap, and
kinetic matrices matched the baseline byte for byte. Matched H2+
parent/PQS/White-Lindsey dimensions, topology, energies, captures, residuals,
warnings, fingerprints, and accounting remained unchanged.

The existing core, public Cartesian, residual-GTO, projected-q-shell/nested,
and matched-release owners passed unchanged. Warm representative-layer cost
changed from `2.197 ms / 9.082 MB` to `2.154 ms / 8.986 MB`; warm complete
matched-comparison cost changed from `22.573 s / 2.535 GB` to
`22.570 s / 2.534 GB`. No material cold or warm regression was observed.

## Generic Staged By-Center Sidecar Retirement

`HP-RETIRE-BYCENTER-SIDECAR-FN-01` is deferred with no execution grant.
`HP-RETIRE-BYCENTER-SIDECAR-TEST-01` records the completed preflight and has no
continuing test grant. Main retains the generic staged by-center sidecar.

The initial main-only scan at baseline
`01ad223d3c822c54289d5d8843b3d71adbf4c520` found no tracked installer. The
required sibling scan then found an active consumer in the same repository's
`high-order/manager-lane` worktree at
`ed43ff241b16d1e95ea258843017c6638166a940`:

- `src/high_order_integrated_recipe_experimental.jl` installs the generic
  sidecar when a nonfactorized recipe would otherwise use the dense fallback;
- four integrated-recipe paths call that installer boundary;
- `test/diatomic/runtests.jl` constructs and validates the generic carrier and
  `:staged_factorized` route, including malformed-input rejection;
- the lane's current architecture explicitly retains generic staged units for
  genuinely non-product high-order representatives as part of its intended
  future mainline import.

The branch carries its own implementation, so a main deletion would not change
that worktree's immediate runtime. It would nevertheless create a known merge
and architecture incompatibility. The branch is therefore a real sibling
dependency, not archived evidence, and it may not be excluded merely to make
the main deletion appear caller-free.

### Deferred Decision Boundary

Retain unchanged on main the generic carrier, its range/builder/attach closure,
selector and path branches, nuclear-contraction overload, and current
`:staged_factorized` documentation. Also retain every fixed-block,
product-staged, factorized-final, dense-fallback, annulus, high-order, and
deferred route owner. No source, test, current algorithm documentation, or
sibling worktree edit is authorized by these records.

Retirement may be reconsidered only after a separate high-order owner decision
does one of the following:

1. accepts the generic staged sidecar as a maintained import dependency; or
2. migrates every high-order installer, call path, and test to an approved
   replacement while preserving its non-product numerical contract.

Any migration must receive its own authority and validation budget. It may not
be smuggled into a main cleanup pass. After that decision, a fresh caller scan
must cover main, the high-order worktree, and current paper/downstream
consumers before retirement can be reauthorized.

The Pass 574 stop produced no source or test edit and no numerical result to
accept. The prior approximately `214`-line deletion estimate is retained only
as audit evidence, not as an active budget or permission.

## Angular-Balanced Diatomic Geometry

For each shared z-axis diatomic shell step, shellification computes the target
box in physical parent-axis coordinates. In each bond-axis/transverse plane it
compares

```text
longitudinal margin = distance from the outer nucleus to the box end
transverse scale    = distance from the bond axis to the box side
```

The target keeps the longitudinal margin comparable to the selected
transverse scale. When the `x` and `y` transverse scales differ, the smaller
scale is the existing conservative guard. This is the physical
outer-nucleus 45-degree rule; raw index aspect is not its authority.

If an ordinary index-layer shell body underreaches the target along the bond
axis, shellification emits the difference as ordered native
`:angular_z_extension_slab` stacks. The ordinary body plus those stacks, not
the ordinary body alone, realizes target coverage. Planned extensions larger
than `ns` are split into ordered units with thickness `<= ns`.

Native angular-extension metadata remains:

```text
slab_kind = :angular_z_extension_slab
slab_normal_axis
slab_side
slab_thickness
slab_stack_index
slab_stack_count
bond_axis
reference_nucleus_index
angular_balance_rule = :outer_nucleus_45_degree
longitudinal_margin_physical
transverse_scale_physical
angular_extension_physical
```

This classification applies alongside midpoint slabs, planned boundary or
non-boundary angular extensions, and unexpected outer-mismatch fallback slabs.
It does not make any slab a real shell or a direct identity sector.

## Compact Thin Slabs

For both PQS and White-Lindsey, regions of kind
`:direct_midpoint_slab`, `:outer_mismatch_slab`, and
`:angular_z_extension_slab` use the same compact thin-slab lowering from the
same terminal region, public `ns`, native normal axis, thickness, side, stack
facts, and source support. They must never lower as full identity CPBs.

The compact unit-slice scale is

```text
ns x ns x 1
```

after one-dimensional COMX/product contraction, with `1` on the slab normal.
A thickness-`t <= ns` slab is an ordered face stack with retained scale about
`t * ns * ns`. Its support rows remain owned and disjoint; the retained
functions are compact products rather than those rows themselves.

Planned angular extensions are chunked before lowering. An unplanned fallback
slab thicker than `ns` is a policy failure: construction must stop rather than
drop the slab, retain it as identity, invent route-specific lowering, or
silently choose whole-block compression. Slab normal and thickness must come
from native metadata, never from role-string parsing.

Direct nucleus-centered and atom-contact core regions remain identity sectors.
Real complete shells remain family-specific after common shellification.

## Neutral Face Products

Compact slabs and White-Lindsey facets share the route-neutral face-product
primitive. Two active axes use retained one-dimensional contractions while
one or more parent indices are fixed on the normal axis. A thickness-one slab
is one face block; a thickness-`t` slab is an ordered stack of those blocks.

This coefficient assembly belongs to
`CartesianFinalBasisRealization`, not to PQS or White-Lindsey. Both consumers
must reuse the same internal helper. They must not duplicate the product
assembly, relabel slabs as WL boundary strata, or create a PQS-only slab
projection path.

## Inventory Contract

Terminal geometry and scaffold summaries must describe midpoint,
outer-mismatch, and angular-extension slabs as planned compact slab products.
They must not advertise stale direct-identity mappings for those regions.
These summaries are descriptive inventory only: they do not materialize
coefficients, carry Hamiltonian data, define source dimensions, or create a
parallel report or artifact payload.

The user-facing bounded inventory and due-diligence report are owned by
`terminal_shellification_due_diligence.md`; they consume native geometry and
realization facts without becoming numerical authority.

## Source Ownership

Current implementation ownership is limited to:

- `src/cartesian_shellification/terminal_geometry.jl` for common regions,
  angular targets, and native slab stacks;
- `src/pqs_source_box_route_driver_helpers.jl` for narrow common-input caller
  plumbing;
- `src/cartesian_terminal_lowering/selection.jl` and
  `src/cartesian_terminal_lowering/region_contracts.jl` for common slab
  selection and contracts;
- `src/cartesian_retained_units/lower_contract_units.jl` and
  `src/cartesian_retained_unit_transform_contracts/unit_contracts.jl` for the
  compact slab retained unit and transform contract;
- `src/cartesian_final_basis_realization/terminal_face_product_blocks.jl` and
  the PQS/WL terminal realizers for neutral face-stack realization;
- `src/cartesian_terminal_shellification_geometry.jl` for compact internal
  inventory metadata.

The registry remains authoritative for exact file permission.

## Guardrails

This contract does not change public inputs, central-gap/contact policy,
direct-core parity, complete-shell source dimensions or retained policies,
the established angular-resolution scale, artifacts, Hamiltonian operators,
RG/MWG/IDA, solvers, ECP, or Cr2 workflow semantics. In particular, it does
not select a Cr2 longitudinal `L` or promote any Cr2-specific geometry.

Mapped-COMX remains a separate PQS-only opt-in source-span facility, with
ordinary source spans as the default. Common shellification must not branch on
that choice.

Any change that requires route-specific first-step geometry, new report or
artifact fields, full-identity slabs, label-inferred slab geometry, or a new
source/retention policy requires separate authority.
