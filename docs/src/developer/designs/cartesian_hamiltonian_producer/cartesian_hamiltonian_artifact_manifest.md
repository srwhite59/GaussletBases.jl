# Cartesian Hamiltonian Artifact Manifest

This page is the canonical contract for optional provenance groups attached to
ordinary `CartesianIDAHamiltonian{Float64}` JLD2 artifacts. It does not define
a second Hamiltonian object or reader.

## Lifecycle

| ID | Lifecycle | Current boundary |
| --- | --- | --- |
| `HP-HAM-MANIFEST-FN-01` | Implemented | Matrix-order final-basis labels and recipe provenance on base and supplemented facade artifacts |
| `HP-HAM-MANIFEST-TEST-01` | Validation completed; tracked coverage is partial | Existing-reader checks plus accepted direct-JLD2 schema validation |
| `HP-HAM-MANIFEST-SRC-FN-01` | Partially implemented | Native terminal source shells/modes and retained boundary-seed relations only |
| `HP-HAM-MANIFEST-SRC-TEST-01` | Validation completed for the implemented subset; tracked coverage is partial | Accepted direct-JLD2 source-provenance checks |
| `HP-NEST-ART-FN-01` | Implemented | Truthful `nesting` and route provenance |
| `HP-NEST-ART-TEST-01` | Validation completed; tracked coverage is partial | Accepted PQS/WL provenance and readback checks |

The compact manifest landed in `f3ef53efc`, native source shell/mode groups in
`604e2323d`, retained boundary-seed relations in `5d98ffcc6`, and nesting
truth in `584e9d333`. Manager-log Passes 113, 115, 116, and 134 record the
accepted numerical and readback evidence.

## Artifact Layers

The low-level writer in `src/cartesian_ida_hamiltonian.jl` remains minimal. It
writes exactly:

```text
artifact_kind = :cartesian_ida_hamiltonian
format_version = 1
kinetic
nuclear_attraction_unit_by_center
electron_electron_ida
nuclear_charges
nuclear_positions
nup
ndn
```

`read_cartesian_ida_hamiltonian` validates the kind and version and reads only
those Hamiltonian datasets. It ignores every provenance group described below.
That compatibility behavior is intentional.

Facade producers append sidecars after the minimal artifact is written:

- base artifacts add `producer_provenance/`, `hamiltonian_manifest/`,
  `recipe_provenance/`, and `coulomb_expansion/`;
- supplemented residual-GTO/MWG artifacts add `supplement_provenance/`,
  `hamiltonian_manifest/`, `recipe_provenance/`, and
  `coulomb_expansion/`.

The R1 contract owns `producer_provenance/`; the residual workflow owns
`supplement_provenance/`; and the Coulomb policy owns
`coulomb_expansion/`. This page owns only the manifest and recipe groups and
their composition with the unchanged artifact.

Direct callers of `write_cartesian_ida_hamiltonian` do not receive these
sidecars automatically. No public manifest reader is implemented.

## Final-Basis Labels

Every facade-written artifact contains:

```text
hamiltonian_manifest/manifest_version = 1
hamiltonian_manifest/final_basis_labels/
```

`final_basis_labels/` has one row per Hamiltonian row/column in native matrix
order. Its live fields are:

```text
final_basis_col
sector
unit_label
unit_kind
source_region_label
source_region_label_status
source_box_label
source_box_label_status
owner_nucleus_index
owner_label_status
shell_label_status
shell_index
ray_label_status
ray_id
ray_family_label
radial_order_status
radial_order
center_x
center_y
center_z
center_definition
center_status
lowdin_correction_applied
supplement_label
angular_power_x
angular_power_y
angular_power_z
inferred_from_centers
inferred_from_nearest_grid
inferred_from_support_order
inferred_from_support_indices
inferred_from_raw_to_final_support
```

Current source emits `sector = :base` for terminal-basis rows and
`sector = :residual` for appended residual-GTO rows. It does not emit a third
supplement-derived sector. `final_basis_col` is `1:n`; unavailable integer
shell/ray/radial labels are `0`; unavailable angular powers are `-1`; and an
absent physical owner is `0`.

Direct terminal rows use parent-support centers. Support-local Lowdin rows use
an `abs2` coefficient centroid and are marked representative. Residual rows use
the exact diagonal position expectation in the augmented basis, their native
owner, and residual label. These centers describe a row; they do not define
its construction identity.

All five `inferred_from_*` fields are `false`. Missing source-box, shell, ray,
radial, or owner facts remain status-bearing `:unavailable` or `:mixed`. A
writer must never manufacture them from centers, nearest-grid snapping,
support order, support indices, or raw-to-final support.

Matrices and sector ranges remain in native order. This manifest does not
authorize a z-sorted matrix convention or any row permutation.

## Native Source Provenance

When terminal retained-rule facts are available, the writer also emits:

```text
hamiltonian_manifest/source_shells/
hamiltonian_manifest/source_modes/
```

### Source shells

The live `source_shells/` fields are:

```text
status
schema_version
row_count
source_shell_id
unit_label
unit_kind
final_basis_start
final_basis_stop
source_shell_label
construction_kind
axis_start
axis_stop
contracted_dims
source_mode_count
source_mode_ordering
center_definition
center_status
lowdin_correction_applied
shell_label_status
ray_label_status
radial_order_status
inferred_from_centers
inferred_from_nearest_grid
inferred_from_support_order
inferred_from_support_indices
inferred_from_raw_to_final_support
```

The current status is `:native_terminal_source_shells`, schema version is `1`,
and axis intervals and contracted dimensions are three-column integer tables.
Shell labels are native. Source-shell centers, ray labels, and radial order are
currently unavailable.

### Source modes

The live `source_modes/` fields are:

```text
status
schema_version
row_count
source_shell_id
mode_index
unit_label
native_source_id_label
local_axis_x
local_axis_y
local_axis_z
center_x
center_y
center_z
center_definition
center_status
lowdin_correction_applied
source_mode_status
shell_label_status
ray_label_status
radial_order_status
inferred_from_centers
inferred_from_nearest_grid
inferred_from_support_order
inferred_from_support_indices
inferred_from_raw_to_final_support
```

A source-mode identity is
`(source_shell_id, local_axis_x, local_axis_y, local_axis_z)` in shell-local
coordinates. It is a construction label, not a coefficient map, parent row,
support row, positive-weight claim, or reconstructable transform. Current
mode centers are `NaN` with unavailable status; ray/radial labels are also
unavailable. Parent-lattice coordinates are not written.

### Final-basis source relations

`hamiltonian_manifest/final_basis_source_relations/` is written only when at
least one retained boundary-seed relation exists. Its live fields are:

```text
final_basis_col
relation_index
relation_kind
source_shell_id
source_mode_label
local_axis_x
local_axis_y
local_axis_z
relation_status
shell_label_status
ray_label_status
radial_order_status
coefficient_status
weight_status
span_status
inferred_from_centers
inferred_from_nearest_grid
inferred_from_support_order
inferred_from_support_indices
inferred_from_raw_to_final_support
```

Current relations are only `:boundary_mode` rows with
`:native_retained_boundary_seed` status. Coefficients and weights are marked
`:not_serialized`; spans, ray labels, and radial order are unavailable.

This is the exact partial implementation boundary. Direct/support-dense rows
and residual rows have no source-relation rows. No source coefficient, weight,
span, ray, cone, or radial grouping has been implemented. The remaining
approved seam may carry only additional construction-native facts through the
existing compact `source_mode_provenance` context. It may not introduce center
inference, consumer locality policy, dense transforms, or route reports.

## Recipe Provenance

Every facade-written manifest contains `recipe_provenance/` with these live
keys:

```text
provenance_version = 1
producer
nesting
route
ns
q
q_rule
ns_source
core_spacing
padding
s_factor
mapping_s_factor
mapping_s_standard
mapping_s_effective
radius
xmax_parallel
xmax_transverse
extent_source
parent_axis_counts
atom_symbols
nuclear_charges
atom_locations
nup
ndn
basisname
basisfile
lmax
uncontracted
width_filtering
base_dimension
residual_dimension
augmented_dimension
```

Values come from validated producer input and constructed dimensions, not from
route reports, element tables, solver assumptions, or artifact inference.
Base artifacts use no supplement values, zero residual dimension, and equal
base/augmented dimensions. Supplemented artifacts record their validated
supplement and base-then-residual dimensions.

`nesting` records `:pqs` or `:wl`. `ns` is the normalized public source size;
`q` is its route-local derivative; `q_rule` and `ns_source` preserve that
normalization history. Mapping fields record the requested factor and resolved
standard/effective values. `padding`, `radius`, and the extents are provenance,
not permission to change box sizing.

Base artifacts repeat the same truthful `nesting` and base-route facts in
`producer_provenance/`. A writer must not use a PQS label for a WL construction
because an internal helper or historical route name is PQS-oriented.

Route labels are derived from system kind, nesting, and supplement state:

```text
:one_center_pqs_base
:one_center_wl_base
:z_axis_diatomic_pqs_base
:z_axis_diatomic_wl_base
:one_center_pqs_residual_gto_mwg
:one_center_wl_residual_gto_mwg
:z_axis_diatomic_pqs_residual_gto_mwg
:z_axis_diatomic_wl_residual_gto_mwg
```

Later composition authorities made all listed construction cells reachable.
That support does not originate from `HP-NEST-ART-*`; those IDs own only
truthful persisted nesting and route facts.

## Ownership And Failure Behavior

Implemented source owners are:

- `src/cartesian_base_hamiltonian.jl` for labels, source provenance, recipe,
  route truth, and facade composition;
- `src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl` for the
  supplemented artifact and `supplement_provenance/` write;
- `src/cartesian_ida_hamiltonian.jl` for the unchanged minimal writer/reader;
- `src/cartesian_final_basis_realization/CartesianFinalBasisRealization.jl`
  for the nesting-neutral module boundary description.

Manifest dimensions must equal the Hamiltonian dimension. Missing native
source provenance omits optional source groups; it does not trigger inference.
Empty relation tables are omitted. Existing writer/reader, filesystem, and
dimension errors propagate.

Tracked validation lives in:

- `test/driver_public/cartesian_base_hamiltonian_runtests.jl`;
- `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl`.

Those tests cover artifact/readback and selected provenance roots. Detailed
field/status validation was accepted through the ignored direct-JLD2 probes
recorded in manager-log Passes 113, 115, 116, and 134; it is not a complete
tracked regression of every sidecar key.

## Exclusions

This contract does not authorize:

- changes to Hamiltonian matrix keys or ordinary reader behavior;
- a public manifest reader, new artifact wrapper, or artifact schema dump;
- `T_G`, `T_A`, coefficients, dense transforms, moment matrices, support/raw
  inventories, or reconstructable final-basis coefficients;
- inferred source identity, consumer ray/locality policy, or z-sorted matrices;
- solver, restart, corrected-Hamiltonian, EGOI, protected-localized, ladder,
  Cr2-specific, or consumer-specific fields;
- public input/default changes, box-policy changes, or new composition support.

Protected-localized artifacts and their row-locality metadata are a distinct
kind and convention governed by
[the protected artifact contract](protected_localized_artifact.md). They are
not ordinary manifest variants.
