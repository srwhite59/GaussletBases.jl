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

Normalized `ns` is the common shell-size input. The direct nucleus-centered
core side is

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
from `src/cartesian/cartesian_nested_diatomic.jl`. No deleted fact was replaced by a
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
`src/cartesian/cartesian_nested_faces.jl` and the two-line propagation in
`src/cartesian/cartesian_nested_diatomic.jl`. No tracked test, driver, tool, example, or
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
src/cartesian/cartesian_nested_faces.jl
src/cartesian/cartesian_nested_diatomic.jl
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

`HP-RETIRE-BYCENTER-SIDECAR-FN-01` records the accepted retirement of the
unused generic staged by-center sidecar. `HP-RETIRE-BYCENTER-SIDECAR-TEST-01`
records its completed validation. Both grants are closed.

### Owner And Caller Closure

The high-order owner decision recorded in
`chatarchive/reports/software_reviews/high_order_lane_closure_and_generic_sidecar_release_2026-09-02.md`
classifies `high-order/manager-lane` at
`ed43ff241b16d1e95ea258843017c6638166a940` as archived, not as a maintained
future merge target, and explicitly releases its generic-sidecar claim. Its
integrated recipe and diatomic tests are historical callers. The worktree must
remain physically present until its separate `523 MB` ignored-scratch
preservation gate completes; that storage requirement grants no source
compatibility claim and this retirement must not edit or delete the worktree.

A fresh audit at main baseline
`b4e632f1417f15e7b03b4b9a165ebf519302ce22` found:

- no tracked main constructor, test, bin, example, tool, or driver installs the
  generic carrier;
- the only main references are its private definition/build/attach closure,
  selector and path branches, matching nuclear overload, and current
  documentation;
- active paper and external-work trees contain no reference to the carrier,
  installer, or `:staged_factorized` route;
- one current CR2 smoke inspects `_nested_by_center_sidecar_path`, but consumes
  only the separately preserved factorized-final and product-staged behavior;
- registered detached/release worktrees and old merged, temporary, archive, or
  roadmap refs contain frozen source copies or obsolete planning code, not a
  non-archived installer or current compatibility consumer.

The generic sidecar is therefore orphaned after the accepted owner release.
The product-staged carrier remains installed by current PQS/endcap-panel
construction and is outside this retirement.

### Accepted Deletion

Commit `f13e946dfffd47ca9bc06e003f76f54c25316415` removed exactly this
private closure:

```text
_CartesianNestedStagedByCenterSidecar3D
_nested_default_staged_sidecar_column_ranges
_nested_build_staged_by_center_sidecar
_nested_attach_staged_by_center_sidecar!
the generic carrier branch in _nested_staged_by_center_sidecar
the :staged_factorized branch in _nested_by_center_sidecar_path
the _CartesianNestedStagedByCenterSidecar3D overload of
  _qwrg_bond_aligned_staged_by_center_nuclear_one_body_by_center
```

It corrected only the directly stale current documentation that presented
`:staged_factorized` as an available route or said that the nuclear method had
both generic staged and product-staged implementations. Historical archived
design text, frozen releases, and the retained high-order worktree may keep the
retired vocabulary.

The retirement preserved unchanged:

- fixed-block sidecar storage and `_nested_staged_by_center_sidecar_cache`;
- `_CartesianNestedProductStagedByCenterSidecar3D`, all product-staged
  builders and dispatch, and the active product-staged nuclear method;
- `_nested_nonzero_coefficient_rows` and
  `_qwrg_contract_staged_nuclear_block`;
- the factorized-final route and general dense fallback;
- all current PQS, ordinary-QW, nested, mapped-COMX, and matched-H2+ behavior;
- the archived high-order worktree and its ignored or uncommitted evidence;
- the deferred XY, shell-plus-core, hierarchical, annulus, and other
  high-order material.

No replacement type, tuple, named tuple, adapter, alias, deprecation, cache,
metadata field, status vocabulary, helper, or compatibility path was added.
Git history is the main-source archive.

### Accepted Scope And Accounting

Implementation source was confined to:

```text
src/cartesian/cartesian_nested_faces.jl
src/ordinary/ordinary_qw_raw_blocks.jl
```

Direct stale-document correction was confined to:

```text
docs/src/algorithms/cartesian_nested_endcap_panel_shared_shell.md
docs/src/developer/cartesian_parent_factors_and_cpb_kernels.md
```

The exact source delta was `+0/-210`. The two current-document corrections
added three and removed four lines. No committed test, fixture, probe, new
source file, API, export, dependency, workflow, version, or release artifact
was added or edited.

### Accepted Validation

A repeated main and non-archived-consumer caller scan found no live installer.
A transient probe outside the repository froze the baseline and compared the
candidate for:

- ordinary-QW factorized-final by-center H1;
- product-staged endcap-panel by-center H1;
- every per-center matrix, their sum, and the complete one-body matrix;
- unchanged route selection as `:factorized_final` and
  `:product_staged_factorized`, respectively;
- product-staged agreement with the general dense oracle at the existing
  accepted tolerance.

All baseline-to-candidate matrices agreed bitwise. The factorized and
product-staged dimensions remained `433` and `397`, and the product-staged
dense-oracle error was `8.88e-16`. Representative factorized time changed from
`49.33 ms` to `49.98 ms`, with allocation falling by `16` bytes; product-staged
time changed from `56.01 ms` to `59.12 ms`, with allocation unchanged. These
measurements show no material performance or allocation regression.

Existing core, public Cartesian (`232/232`), nested PQS (`464/464`) and
supplemented facade (`64/64`), matched H2+ (`18/18`), and documentation
(`157/157` and `10/10`) owners passed unchanged. Package load, authority checks,
Documenter, CI run `33664855002`, and Docs run `33664855006` also passed. The
archived high-order worktree and its evidence were not modified.

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

- `src/cartesian/cartesian_shellification/terminal_geometry.jl` for common regions,
  angular targets, and native slab stacks;
- `src/cartesian/pqs_source_box_route_driver_helpers.jl` for narrow common-input caller
  plumbing;
- `src/cartesian/cartesian_terminal_lowering/selection.jl` and
  `src/cartesian/cartesian_terminal_lowering/region_contracts.jl` for common slab
  selection and contracts;
- `src/cartesian/cartesian_retained_units/lower_contract_units.jl` and
  `src/cartesian/cartesian_retained_unit_transform_contracts/unit_contracts.jl` for the
  compact slab retained unit and transform contract;
- `src/cartesian/cartesian_final_basis_realization/terminal_face_product_blocks.jl` and
  the PQS/WL terminal realizers for neutral face-stack realization;
- `src/cartesian/cartesian_terminal_shellification_geometry.jl` for compact internal
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

## Finite Collinear PQS

Pass 630 accepts Pass 629 implementation cf608bd0e151f95ea2f88960959f463184325151
and input normalization 3e32d7e50e9a232d8545fa313d2f06ac83565031.
HP-COLLINEAR-PQS-FN-01/TEST-01 are maintenance only. This additive expert
finite-chain producer is the explicit exception to the preceding no-new-geometry
guardrail; atom/diatomic geometry and public defaults remain unchanged.
The following frozen implementation boundaries remain maintenance constraints,
not authority for another capability or another 450-line extension.

### Public boundary and storage

Authorize exactly two new exports, documented at their definitions:

- `cartesian_collinear_working_basis(z, Z; core_spacing, transverse_spacing,
  padding_parallel, padding_transverse, core_side, angular_reference_count,
  outer_face_count, tail_spacing, angular_resolution_scale, expansion)`.
  All listed controls are required. Use G10, reference spacing 1 and existing
  PGDG mapped-axis construction; no alternative backend, mapping or preset.
- `cartesian_collinear_operators(working, z, Z; expansion)` returns exactly
  `(; one_body, electron_electron_ida, nuclear_repulsion)`: two complete dense
  Float64 matrices and the scalar nuclear repulsion. The explicit operator
  nuclei/charges define the potential on the supplied, unchanged working basis.

Validate nonempty equal-length real vectors, finite strictly increasing z,
finite positive charges, positive finite spacing/padding/tail/angular controls,
positive integer counts (not Bool), and odd core_side. Physical centers are
distinct although snapped sites may coincide. General positive charges are
allowed only where existing mapping validity checks succeed. No clipping,
charge replacement, amplitude repair or silently substituted mapping.
Use actual min/max z plus requested end padding; retain current coincident
transverse-center rules. Do not confuse requested padding with realized bounds.

Allow one private, unexported bare working type with exactly terminal basis and
parent axis bundles; no molecule-size-dependent type parameters or inventories.
Existing `gto_overlap_matrix` and `import_external_gto_orbitals` accept it by
delegating to existing terminal/parent mixed-overlap mathematics and packet
validation. Preserve block selection and finite-probe validation. Raw import is
X*C without renormalization or determinant cleanup. No residual/injected carrier,
global parent-by-final coefficient map, general representation engine or solver.

Accumulate each charge-weighted nuclear contribution into the kinetic matrix
using existing kernels and lexical buffers. Check expansion/axis exponent
agreement with the existing validator. Do not allocate a mandatory list of
per-center dense matrices, change CartesianIDAHamiltonian, place a sum in a
fake center slot, or claim its reweighting/artifact semantics. Retain existing
positive IDA integration-weight validation and finite output checks.
The result is small-chain numerical data, not a general Hamiltonian facade.

### Ownership and retained space

1. Start from equal-index-width odd cores snapped to the mapped parent. Reject
   cores outside the parent. Merge their overlap connected components
   simultaneously, preserving each initial component hull as a direct sector.
2. Expand all current boxes by one index on every axis only while all proposed
   boxes fit. Determine all proposed overlap components before updating any
   group; no greedy pair mutation or odd/even recipe.
3. For each component, I is the hull of its previous boxes and O the hull of
   proposed boxes. Emit I minus previous boxes as direct contact sectors and
   O minus I as the existing PQS complete shell. Do not contract cores/contacts.
4. At termination, assign interior gaps between surviving groups directly:
   hull(groups) minus union(groups). Then decompose parent minus hull(groups)
   into the existing axis-ordered outer slabs. Never confuse interior gaps with
   exterior faces or drop them when multiple groups remain.
5. Consolidate outer_mismatch_pieces into one callable private routine in its
   present owner, replacing the nested implementation. Preserve exact existing
   diatomic output, metadata, order and numerical behavior. Earlier axes own
   edges/corners; later axes use earlier axes' inner intervals.
6. Apply the unchanged compact thin-slab kernel: direct normal indices and
   outer_face_count retained functions on each in-plane axis. Validate the
   positive count against BOTH actual in-plane lengths before realization.
   Reject incompatibility; no clamp, direct-completion fallback or map change.
7. Use canonical ascending-z group order and x-major/y-major/z-fast support
   order. Inventories are vectors; spatial triples alone have fixed tuple shape.
   Production retains local blocks/ranges, not dense occupancy or Q/S oracles.

Reuse the all-nucleus angular-spacing dimension calculation, boundary COMX
selection and existing shell Lowdin/gauge/weight checks unchanged.
angular_reference_count calibrates that calculation; it is NOT a constant
source q or the existing H/H2 q/ns contract. Record actual source shapes in
acceptance evidence. outer_face_count is separate from angular calibration,
parent spacing and padding. Counts 9 and 7 below are fixture controls, not
chemically certified defaults. Scalar face counts cannot independently refine
a long in-plane direction beyond the shorter interval; reject requests that
need per-axis counts/subdivision rather than broadening this packet.

### Frozen evidence and acceptance

Reviewed reports (ignored evidence, not execution authority):
`tmp/reviews/collinear-pqs-design-2026-09-14/REPORT.md`, SHA-256
b27161c0cbb5ab65976a7076b97bde17b71c8f00ee4e003957a4ef4897a50880;
`tmp/reviews/collinear-pqs-completion-2026-09-14/REPORT.md`, SHA-256
5c2e4bd47fc3058deefbdd7e45e581ee1f7a45445832de432dea4dfb4f0f6a40.
Their passing logs and independent parent-action references were inspected.
No broad numerical replay is required for authorization.

Reuse uniform H3 (-2.4,0,2.4), H4 (-3.6,-1.2,1.2,3.6),
unequal H3 (-2.4,-1.2,2.4), contact H3 (-1.2,0,1.2), and H-He-H
(-2.4,0,2.4; charges 1,2,1). Report controls: core_side=3,
longitudinal spacing .6, tail_spacing=2.8, angular scale=1.4,
reference counts 3/5, compact45 expansion. Baseline transverse spacing/padding
are .6/3 with longitudinal padding 3. Preserve initial/growth components,
existing inner blocks/columns and source shapes; compact exterior replaces
the direct remainder, so do not assert unchanged total direct-baseline size.

Freeze enlarged H3 at transverse .45/padding 6, reference count 5,
outer count 9: parent 13x13x17, inner 617, final 1265.
Freeze unequal H3 at transverse .6/padding 3, reference count 5,
outer count 7: parent 9x9x15, inner 423, final 619.
Use direct-completion scratch references only as test oracles, never as an
alternate production mode. Complete operators must be exercised, not merely
selected blocks; at least the enlarged q_outer=9 full output must match the
independent parent-action projected checks.

Numerical acceptance is frozen as follows (absolute tolerances, rtol=0):
- Discrete coverage, ordering, group membership, direct core/contact coefficients
  and existing inner blocks/columns: exact. No missing or duplicate parent sites.
- Max-entry sampled-parent orthogonality error <=1e-10; finite positive IDA
  weights with the existing 1e-14 validity floor unchanged.
- H1/IDA symmetry and independent finite-expansion H1/IDA contraction agreement
  <=1e-10; nuclear repulsion agreement <=1e-12 Ha. Use actual positive charges.
- GTO overlap vs independent Q'X_parent <=1e-12; importer vs its X*C exactly.
- Origin-centered normalized s, px, py probes at exponents .1 and .4:
  compare each component separately, never average px/py. For the two frozen
  compact configurations, match logged captures within 1e-8 and projected
  H/IDA entries within 1e-8 Ha. Other fixtures use independent comparisons
  in their actual retained space, not direct-completion total capture targets.
- At exponent .1, enlarged compact capture is
  (.9999318254,.9997143222,.9997136299); unequal compact capture is
  (.9988312030,.9911419747,.9911418013).
- Added direct-to-compact max capture/H/IDA-self changes at exponent .1:
  enlarged <=1.1e-6 / 4.4e-6 Ha / 9.3e-6 Ha;
  unequal <=8.9e-6 / 1.1e-5 Ha / 6.1e-5 Ha.
  At exponent .4: enlarged <=1.7e-7 / 3.0e-6 Ha / 6.3e-8 Ha;
  unequal <=1.2e-7 / 3.4e-6 Ha / 4.7e-8 Ha.
These are fixture parity/approximation gates, not chemical accuracy promises.
Preserve parent loss, inner retention loss, outer truncation loss and raw
projected-operator changes as separate quantities. None is an energy-error
estimate; IDA is not exact four-index ERI. Preserve observed transverse
asymmetry, actual bounds and snap errors. No new persistent diagnostics schema.

Add one focused multiple-group termination case: index parent 7x7x25,
core_side=3, nuclear indices (4,4,5),(4,4,13),(4,4,21).
Two synchronous expansions leave three disjoint 7x7x7 boxes, 1029 sites;
termination leaves 98 direct interior-gap sites (z=9,17) and 98 exterior
sites (z=1,25). Total 1225 exactly once. Outer count 3 fits the two z-normal
7x7 faces. Check decomposition, direct gaps and slab realization using a
bounded existing-factor or synthetic orthogonal-factor fixture. This is not a
new geometry campaign. Also reject invalid outer counts and mapping-invalid
inputs without fallback; preserve tests for simultaneous/transitive merges.

Production uses only local basis storage. Dense global Q/S and independent
tensor-action oracles are allowed only in bounded tests/scratch.
For the enlarged acceptance fixture, record fresh/warm construction, transfer,
and complete-operator time/allocation separately. Complete operator build
ceilings: Nfinal<=2300, retained basis/parent/output <=512 MiB, cumulative
allocation <=4 GiB, first call <=180 s and warm <=120 s on the documented
reference environment. Do not turn machine-dependent times into universal
CI assertions or infer end-to-end cost by summing unrelated compiled phases.
Dense outputs remain quadratic; support-pair work still depends on parent size
and center count. No long-chain scaling or many-body claim.

### Exact surfaces and budgets

Only these six existing source files may change (added lines include docstrings
and relocated decomposition, not net delta):
- src/cartesian/cartesian_shellification/terminal_geometry.jl: preferred 160.
- src/cartesian/cartesian_base_hamiltonian.jl: preferred 125.
- src/cartesian/pqs_source_box_route_driver_helpers.jl: preferred 70.
- src/cartesian/pqs_source_box_low_order_materialization.jl: preferred 40.
- src/cartesian/cartesian_gto_probes.jl: preferred 20.
- src/GaussletBases.jl: preferred 2, only the two exports.
Total preferred 417; USER-APPROVED HARD EXCEPTION 450 added source lines,
for this feature only. Existing kernel bodies and all other source owners
remain unchanged; no source file, public type, cache or framework is added.

The independently selected existing test owner is
test/driver_public/cartesian_base_hamiltonian_runtests.jl, already in public
Cartesian/Supported-floor CI. Preferred/hard added test lines: 170/220,
including compact local oracles, malformed-input and multi-group checks.
No other test edits, new test owner, fixture file or runner change.
Tests must embed only compact reference values or compute bounded independent
oracles; CI must not read ignored reports or machine-local scratch files.
Reader guidance only in docs/src/manual/projected_q_shells.md and
docs/src/reference/export.md: preferred/hard total additions 40/55; curate both
docstrings, clearly label expert finite chains and required convergence controls.
Do not duplicate manuals or imply changed stable/released documentation.

Run the focused owner, existing residual-GTO public owner, matched H2+ release
owner once, relevant common-shell/compact-slab owners, package load,
docs_fast/full docs, authority/self-test, generated parity, Documenter and diff
checks. Existing atom/diatomic fingerprints and released interfaces must remain
unchanged. Normal source-bearing CI must classify full and pass all three
unchanged jobs plus Docs. No full angular suite, repeated paper example,
workflow change or benchmark campaign.

Delete duplicated outer decomposition when making it callable; do not create a
parallel completion implementation. The recursive-chain retirement conditions
are satisfied for the separately bounded Pass 637 deletion below; this original
construction grant itself authorizes no retirement. Preserve all
released ordinary-chain interfaces and the distinct sliced-chain capability.
No periodic/off-axis extension, solver, screening, artifact change, new
contraction mathematics, legacy revival, de-export, retirement, Standard60,
release or stable promotion is granted.

Failure rule: if exact ownership, frozen approximation gates, storage,
unchanged atom/diatomic behavior or the 450-line limit cannot be met, or broader
semantics/kernel edits are needed, make no implementation commit and report the
specific obstacle. No clamp, hidden fallback, tolerance relaxation or planner-only
delivery. Close implementation separately after full acceptance evidence.

## Private Recursive Chain Retirement

Pass 637 grants HP-RETIRE-RECURSIVE-CHAIN-FN-01/TEST-01 for pure deletion at
baseline 5a770a571d9527838aaaf540fed515cd03380b97. The existing public Cartesian
owner already covers H3/H4, unequal gaps, contact, mixed charges, matched-parent
operators, raw transfer, separate s/px/py capture and three unmerged groups with
1029 group +98 interior-gap +98 exterior sites, exactly once. These construction
conditions do not require H10 convergence, screening or correlated DMRG.

Independent tracked-main and bounded HFDMRG, paper, hchain, grid-study and codex
tree scans found no executable consumer outside the private closure. No dynamic
symbol construction occurs in its owner. This is bounded local clearance, not
proof about unpublished downstream code. The high-order checkout at ed43ff241
retains old callers and its own implementation. The September 2 high-order owner
closure classifies that lane as archived, not a maintained future import target;
its sidecar permission alone does not authorize this distinct deletion. Preserve
the archived worktree and all its evidence. Stop on any new non-archived caller.

Only src/cartesian/cartesian_nested_experimental_geometries.jl may change.
Delete baseline line ranges 1-63, 135-179, 227-230, 266-857 and 1482-1485,
including separator lines: exactly 708 source deletions, zero additions.
These are three private chain carriers with docstrings, three display methods,
the chain contract-audit overload, odd/even policy and candidate selection,
chain diagnostics/recursive nodes/source construction, and the chain fixed-block
overload. No replacement, alias, helper rename, fallback or relocated copy.

Frozen baseline SHA-256:
`26aba50390ce4ecb66516b700cc1c1097270d1fc20f95e5537d5744daf452c68`.
Exact concatenated deleted ranges SHA-256:
`1263baccc959ecef3b032652cb96ea34bdaf6d17f295a23e22d527772b11e0de`.
Required complete retained file (780 lines) SHA-256:
`43b987efc0478a026ac3de8fd8d70d579ad9e7a495e12ea7d1853592657ecb46`.
Preserve baseline lines 235-265, including the separator, for
`_nested_chain_three_child_boxes` SHA-256:
`a532bdb18b62e081c42876b969dc3c84ca2fcf61252a893f9a198a9497840844`.
The square-lattice ternary candidate calls this helper; preserve all square
carriers, methods and the mixed file byte-for-byte outside the deleted ranges.

Preserve ordinary product-chain exports, existing generic/prebuilt fixed-block
dispatch, sliced-chain capability, multisliced evidence, current PQS and all
atom/diatomic behavior. Released private definitions do not make these symbols
public exports; do not touch any release or tag. The sole developer symbol
reference in cartesian_nested_decomposition_plan.md is explicitly historical:
leave it and all archived evidence unchanged. No reader-document edit is needed.
Implementation scope is one source file; test/documentation additions are zero.
Normal lifecycle records and generated views belong to separate closeout.

Acceptance: repeat caller scans; verify exact removal/retained hashes and
unchanged root exports; load the package; run unchanged test/core/runtests.jl
through its existing runner and test/driver_public/cartesian_base_hamiltonian_runtests.jl.
These cover preserved ordinary chain/square geometry and finite-collinear
contracts. A bounded ignored direct check of the retained three-child helper
on x/y/z axes is permitted; no committed test or new test owner. No need to
rerun archived recursion, H10, supplementation, screening or angular campaigns
solely for orphanhood. Run docs_fast/full docs, authority/self-test, generated
parity, Documenter, log/diff checks, and unchanged source-bearing three-job CI
plus Docs. CI must classify full; do not change the classifier or selected jobs.

Failure rule: if any live caller, changed kept byte, required test adaptation,
added source line or broader dependency is found, make no implementation commit
and report. No source decomposition, other retirement, numerical policy,
workflow, dependency, API, export, release, stable or consumer-assignment change.
Repo-manager waits for this recorded grant and required checks before deletion.
