# Cartesian Hamiltonian Artifact Manifest

Status: approved for implementation under `HP-HAM-MANIFEST-FN-01` and
`HP-HAM-MANIFEST-TEST-01`.

The optional source-mode provenance seam is approved separately under
`HP-HAM-MANIFEST-SRC-FN-01` and `HP-HAM-MANIFEST-SRC-TEST-01`.

## Purpose

Canonical-driver Hamiltonian artifacts need enough sidecar information for
downstream consumers to identify matrix rows, locality, freezing groups, and the
public construction recipe. This is artifact annotation only. The Hamiltonian
object, matrix datasets, and `read_cartesian_ida_hamiltonian` behavior stay
unchanged.

## Approved Boundary

Approved source files:

- `src/cartesian_base_hamiltonian.jl`;
- `src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl`;
- `src/cartesian_ida_hamiltonian.jl` only for a small unexported sidecar
  writer/helper if local duplication would otherwise be worse.

The implementation may write additional JLD2 sidecar groups after writing the
existing `CartesianIDAHamiltonian{Float64}` artifact. It must not rename,
remove, or reinterpret the existing matrix keys:

- `kinetic`;
- `nuclear_attraction_unit_by_center`;
- `electron_electron_ida`;
- `nuclear_charges`;
- `nuclear_positions`;
- `nup`;
- `ndn`.

`read_cartesian_ida_hamiltonian` must keep reading only the Hamiltonian object
and must ignore sidecar groups. Manifest readback is validation-only through
direct JLD2 inspection or local test helpers, not a public reader API.

## Final Basis Manifest

Approved root group:

```text
hamiltonian_manifest/
```

The manifest reuses the fixed-column/source-mode provenance model from
`docs/src/developer/projected_q_shell_policy.md`. A basis-function identity is
a construction label with explicit status, not a center. Centers are
representative metadata only.

Approved root key:

```text
manifest_version = 1
```

### `final_basis_labels/`

Approved required subgroup:

```text
hamiltonian_manifest/final_basis_labels/
```

This subgroup generalizes the earlier private `fixed_column_labels` table from
"fixed columns" to "final basis columns." It has one row per Hamiltonian matrix
row/column.

Approved minimal fields:

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

Shape and status contract:

- `final_basis_col` is `1:n` in the exact matrix row/column order;
- every field has one value per final basis column;
- `owner_nucleus_index` uses one-based physical nucleus indices and `0` when
  no owner is meaningful;
- unavailable integer labels such as `shell_index`, `ray_id`, and
  `radial_order` use `0`;
- unavailable angular powers use `-1`;
- labels that are not native must use status values such as `:mixed` or
  `:unavailable`, not guessed values;
- all `inferred_from_*` flags must be `false` for production manifest rows;
- `center_*` fields are representative metadata and must not be used as
  identity labels.

Allowed `sector` values are:

- `:base`;
- `:residual`;
- `:supplement_derived`.

For the current residual-GTO/MWG path, residual rows are `:residual`; their
supplement origin is recorded by owner, supplement labels, and angular powers
where those are meaningful. `:supplement_derived` is reserved for future
non-residual supplement-derived functions if approved separately.

`source_region_label` or `unit_label` is the honest fallback grouping when
shell/ray/radial labels are unavailable. Do not infer shell, ray, radial, or
source-box labels from representative centers, nearest-grid snapping, support
order, support indices, or raw-to-final support.

### `final_basis_source_relations/`

Approved optional subgroup:

```text
hamiltonian_manifest/final_basis_source_relations/
```

Use this subgroup only for relation facts the construction natively defines.
Approved fields:

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

Relations may describe source-mode identity, product-axis tuples, boundary
modes, support spans, or Lowdin mixtures only when the construction producer
defines that relation honestly. Relation weights or spans may appear only when
they are construction-native scalar facts, not dense transforms.

### `source_shells/` And `source_modes/`

Approved optional subgroups:

```text
hamiltonian_manifest/source_shells/
hamiltonian_manifest/source_modes/
```

These subgroups are the provenance-compatible extension path from the prior
PQS sidecar contract. A source mode identity is
`(source_shell_id, local_axis_x, local_axis_y, local_axis_z)` in shell-local
axis coordinates. That identity labels a construction-native source function;
it does not claim one parent grid row, one final support row, positive weight,
or a reconstructable basis transform.

When present, `source_shells/` may contain stable shell IDs, unit links,
construction kind, axis intervals, contracted dimensions, source-mode ordering,
center definition/status, and shell/ray/radial label statuses. `source_modes/`
may contain source-shell-local mode indices, native source labels, local axis
coordinates, parent-lattice coordinates when native, representative center
metadata, `lowdin_correction_applied`, and status flags.

If these source layers are not natively available for a row, the manifest must
say `:unavailable` or `:mixed`; it must not silently infer the missing source
identity.

Approved sources for representative centers and labels are the existing
terminal basis blocks, parent axes, residual metadata, and augmented moment/MWG
descriptors already computed by the approved construction. If a source pass
cannot derive a center convention or construction label from those objects
without adding algorithmic metadata, it must stop and report the missing seam.
Do not add coefficients, dense transforms, raw inventories, route reports, or
status payloads to force the manifest to exist.

## Source-Mode Provenance Seam

The first compact manifest writer may write only `final_basis_labels/` and
`recipe_provenance/` when source-mode facts are not live at the artifact seam.
`HP-HAM-MANIFEST-SRC-FN-01` approves one narrow construction-native provenance
carrier so later source work can populate optional groups:

```text
hamiltonian_manifest/source_shells/
hamiltonian_manifest/source_modes/
```

and, where the relation is native, improve:

```text
hamiltonian_manifest/final_basis_source_relations/
hamiltonian_manifest/final_basis_labels/
```

Approved carrier concept:

```text
terminal lowering / retained-unit / raw-product source plans
-> compact source-mode provenance object
-> base working basis manifest context
-> artifact sidecar writer
```

The preferred live carrier is one internal `source_mode_provenance` field on
the object returned by `cartesian_base_working_basis(...)`. A source pass may
instead attach the same compact object to `CartesianTerminalBasisRealization`
only if that avoids duplication or loss of terminal construction ordering. In
either case, the object is artifact provenance only; Hamiltonian construction,
operator assembly, residual selection, and reader behavior must not consume it.

Approved object contents are row tables matching the optional manifest groups:

- source-shell rows with stable shell IDs, unit links, construction kind,
  source-box/source-region labels, source-mode dimensions/order, source
  intervals, center definition/status, Lowdin-correction status, and
  shell/ray/radial label statuses;
- source-mode rows with `(source_shell_id, local_axis_x, local_axis_y,
  local_axis_z)` identity, source-mode label/order fields, parent-lattice
  coordinates only when native, representative center metadata, and status
  flags;
- final-basis source-relation rows only for native relation facts, such as
  direct identity, boundary source-mode selection, product-axis tuple, or
  explicit mixed/unavailable relation status;
- optional final-basis label improvements only where the final basis column is
  directly and natively tied to a unit/source mode.

The object must not contain coefficients, dense transforms, `T_G`, `T_A`,
raw candidate inventories, raw pair inventories, route reports, allocation
probes, or diagnostic payloads. It may carry compact row IDs, labels, small
integer coordinates, source-mode dimensions, status symbols, and booleans.

Ray and radial labels remain unavailable unless already natively defined by
the construction. Do not add a repo-chosen ray/cone/radial grouping policy in
this lane.

The seam may reuse facts already present in:

- terminal lowering contracts and their source CPBs;
- retained-unit records;
- retained-unit transform contracts and raw product source metadata;
- `RawProductBoxPlan` / `PQSBoundaryProductModeRetainedRule`;
- terminal basis blocks and their matrix-order column ranges.

It must not add algorithmic metadata upstream merely to satisfy an artifact
consumer. If the native facts are missing, write explicit `:unavailable` or
`:mixed` statuses and report the missing producer seam.

## Recipe Provenance Group

Approved group:

```text
recipe_provenance/
```

Approved keys:

```text
provenance_version = 1
producer
nesting
route
q
core_spacing
padding
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

The group may repeat facts already present in `producer_provenance/` and
`supplement_provenance/` so downstream artifact consumers can use one uniform
recipe location. It must be filled from the validated public construction
contract and produced dimensions, not recovered from route reports or inferred
from element tables.

`nesting` records the public construction family (`:pqs` or `:wl`). `route`
records the truthful base route label derived from `(input.kind,
input.nesting)`, for example `:one_center_pqs_base`, `:one_center_wl_base`, or
`:z_axis_diatomic_pqs_base`. The recipe group must not write a PQS route label
for a WL artifact merely because the old route helpers were PQS-named.

For base artifacts without a supplement, supplement-specific keys may be
`nothing` and `residual_dimension` may be `0`; `augmented_dimension` equals the
base dimension. For supplemented artifacts, `basisname`, optional `basisfile`,
`lmax`, `uncontracted`, `width_filtering`, and base/residual/augmented
dimensions must match the validated supplement construction.

## Padding Scope

This lane records the current recipe; it does not change atom box-size policy.
One-center atom padding may be recorded as provenance, but source work under
this lane must not change the current one-center parent-axis counts or reinterpret
atom padding as an active size-control fix. Diatomic padding-derived extents
remain active through the existing facade contract.

## Forbidden

This lane does not approve:

- `T_G`, `T_A`, dense residual transforms, coefficients, dense moment matrices,
  or raw candidate inventories in the artifact;
- allocation probes, benchmark data, route reports, status/result payloads, or
  artifact schema dumps in the driver;
- public reader APIs, public exports, or driver public input changes;
- solver-specific, CR2-consumer-specific, Cr2-specific, ECP, RHF, or
  HFDMRG-specific algorithm fields;
- committed Cr2 fixtures or Cr2-specific branches;
- changes to Hamiltonian matrix keys or `read_cartesian_ida_hamiltonian`;
- new algorithmic metadata in terminal basis, Residual Gaussian, raw-block, or
  route objects.

## Validation

Required validation:

- `git diff --check`;
- package load;
- H atom or H2 base artifact write/readback through the existing reader;
- H2 supplemented artifact write/readback through the existing reader;
- direct JLD2 validation that `hamiltonian_manifest/final_basis_labels/` rows
  match the matrix dimension and use the approved status-bearing fields;
- direct JLD2 validation that missing shell/ray/radial/source labels are
  explicitly `:unavailable` or `:mixed`, not inferred;
- direct JLD2 validation that `recipe_provenance/` records the validated public
  system, basis, supplement, route, parent-axis counts, and dimensions;
- no Cr2 run.

Optional validation:

- Be2 supplemented artifact write/readback and manifest inspection if practical.
- Source-mode seam validation under `HP-HAM-MANIFEST-SRC-TEST-01` may add
  ignored/direct JLD2 checks that optional source groups are present only when
  construction-native provenance rows exist and that missing shell/ray/radial
  labels remain explicitly unavailable.
