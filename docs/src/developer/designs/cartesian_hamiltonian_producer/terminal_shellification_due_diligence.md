# Terminal Shellification Due Diligence

Status: implemented producer/driver reporting under
`HP-DRV-SHELLDD-FN-01` and `HP-DRV-SHELLDD-TEST-01`. The registry owns exact
ID lifecycle and source permission. This contract changes no numerical
construction, source policy, retained policy, or artifact schema.

## Role

The report is the standard bounded review surface for a constructed Cartesian
terminal basis. It exposes normalized system facts, parent axes and weights,
dimension accounting, and terminal shell rows before a consumer interprets an
energy, residual, injection, screened-Hartree, EGOI, Be2, or Cr2 result.

Warnings are advisory diagnostics. A warning does not reject construction or
change a basis unless a separate caller, test, or policy explicitly enforces
it. The report is not a route-stage dump and is not numerical authority.

## Inventory And Reporting Seam

`_cartesian_terminal_inventory_rows(...)` remains the compact terminal
inventory owner. It joins retained units, transform contracts, realized
terminal blocks, parent-axis centers, and native slab metadata into bounded
rows with these existing fields:

```text
region_key
region_kind
lowering_kind
retained_unit_kind
realization_kind
realization_class
shell_index
index_ranges
physical_ranges
support_rows
final_cols
compression_ratio
slab_axis
slab_side
slab_thickness
slab_stack_index
slab_stack_count
```

`_cartesian_terminal_due_diligence_report(...)` extends that inventory from
existing staged producer objects. The base working construction carries the
in-memory report as `terminal_due_diligence`, and the canonical driver prints
it through `print_terminal_due_diligence(...)`. Long center tables are retained
in memory while ordinary output uses a bounded preview.

This seam must not become a new module, route-report framework, persistent
payload, artifact sidecar, or metadata reconstruction path. It may read only
facts already present in the active construction.

## Report Sections

The report has four compact sections:

1. normalized system and geometry context;
2. parent axes, physical box, one-dimensional centers, and gausslet/IDA weight
   statistics;
3. final-basis dimension and compression accounting;
4. terminal region and shell rows in final-basis order.

### System And Geometry

The header exposes derived normalized facts, not merely raw input. The
existing structured context and joined axis rows provide:

```text
geometry_kind
nesting
source_span
atom_symbols
nuclear_charges
nup
ndn
validated_atom_locations
bond_axis
bond_length_physical
box_center_convention
parent_box_physical_bounds
parent_box_physical_lengths
padding_or_radius_resolved_extents
snapped_nuclear_indices
snap_errors_physical
core_spacing / reference_spacing / tail_spacing summary
parent_axis_counts
q
ns
```

For one-center atoms, radius/padding effects on parent counts and physical box
size must be visible. For z-axis diatomics, bond length, bond axis, and
transverse versus longitudinal extents must be visible.

### Parent Axes And Weights

Each `x`, `y`, and `z` axis exposes:

```text
axis
count
physical_min
physical_max
physical_length
center_preview_or_full_bounded_list
min_spacing
median_spacing
max_spacing
nearest_spacing_at_each_nucleus
nearest_index_for_each_nucleus
```

The in-memory axis row retains full one-dimensional centers; bounded printing
records when the preview is truncated. Per-axis nucleus rows carry the
nuclear coordinate, nearest index and center, physical snap error, and nearest
spacing. Core- and tail-region index spans are not emitted by the current
report and must not be inferred from these rows.

Gausslet/IDA weight summaries, where available, use only active construction
weights and expose:

```text
count
sum
min
max
abs_sum
negative_count
near_zero_count
threshold
warning
```

The current `warning` is true for any nonfinite weight, weight below
`-threshold`, or weight with absolute value at most `threshold`. It is not a
large-weight warning. These are gausslet/IDA diagnostics, not residual
integral weights, MWG weights, or proof of quadrature quality.

### Dimension Accounting

The bounded dimension summary exposes, where available:

```text
parent_grid_size
direct_or_core_columns
complete_shell_columns
slab_columns
compact_product_columns
identity_columns
base_final_dimension
supplement_dimension_when_present
residual_dimension_when_present
augmented_dimension_when_present
compression_by_class
large_identity_sector_count
```

Unavailable downstream-stage dimensions remain unavailable; the report must
not reconstruct them from route internals or artifacts.

### Terminal Rows

There is one row per terminal region or shell unit at the granularity needed
to review shellification. The structured row preserves these fields where
available:

```text
terminal_order
terminal_key / region_key
role
region_kind
shell_index
owner/contact/shared classification
index outer_box
index inner_box
index outer_shape
index inner_shape
physical bounds for x/y/z
physical side lengths for x/y/z
physical aspect ratios
source_mode_shape
expected_aspect_balanced_source_mode_shape
source_mode_count
support_rows
retained_count
expected_aspect_retained_count
final_column_range
lowering_kind
retained_rule
realization_rule_or_status
slab normal axis / side / thickness / stack index / stack count
warning_flags
warning_summary
```

Order/key fields connect a row to matrix order. Index and physical boxes make
support geometry explicit. Actual source shape, retained count, and final
column range connect shellification to basis dimension. Native slab fields are
required for midpoint, angular-extension, boundary, and fallback slab rows
when present.

The expected aspect-balanced shape and retained scale are diagnostic
comparisons only. They do not select a PQS longitudinal `L`, overwrite an
authoritative source shape, or define White-Lindsey retention.

## Warning Contract

The bounded warning vocabulary is:

```text
rectangular_physical_shell_cubic_source_modes
expected_source_shape_larger_than_actual
retained_count_below_aspect_balanced_scale
large_identity_sector
missing_shell_index
missing_physical_bounds
missing_source_mode_shape
slab_without_native_metadata
unavailable_expected_shape
axis_center_table_truncated
gausslet_weight_anomaly
padding_or_radius_not_reflected_in_box
```

`warning_flags` holds stable symbols and `warning_summary` is concise human
text. Missing or suspicious facts remain visible, but no flag silently
changes shellification, source dimensions, retained counts, lowering, or
realization. New warning symbols or table fields require separate authority.

## PQS And White-Lindsey Interpretation

The report reviews both construction families through their common terminal
basis boundary, but it does not erase their numerical distinction.

- Common direct/core, shell, slab, owned-support, and physical geometry facts
  come from common shellification.
- PQS complete shells report the actual authoritative source-mode shape from
  the separate aspect-aware `(q,q,L)` construction policy.
- White-Lindsey complete shells report their native boundary-stratum/product
  realization and retained counts; they are not judged by PQS source-box
  policy.
- Thin slabs report the common compact face-stack realization for both
  families and must not appear as direct identity sectors.

The report may compare physical aspect with actual source shape. It must not
choose source dimensions, retained policies, the established angular scale,
or mapped-COMX behavior. Ordinary source spans remain the default;
mapped-COMX remains a separate PQS-only opt-in contract.

## Separate Future Question

Local contracted-COMX nuclear-resolution diagnostics remain a separate future
design question. This report defines no such metric, threshold, warning flag,
or table field, and it does not treat contracted-COMX behavior near a nucleus
as current source-shape or retention policy. Any such diagnostic requires its
own numerical definition, ownership, and authority before entering reporting
or construction.

## Source Ownership And Boundaries

The implemented source surface is:

```text
src/cartesian_base_hamiltonian.jl
bin/cartesian_ham_builder.jl
```

`src/pqs_source_box_route_driver_helpers.jl` may supply compact access to
existing native facts; it is not a second report owner. The registry remains
authoritative for exact source permission.

The report must remain bounded and row-oriented. It must not expose route
skeletons, full source-mode inventories, pair inventories, raw blocks,
all-support-row listings, dense coefficients or transforms, full metadata,
recursive route stages, or artifact internals. It must not change public
inputs, artifact schema/readers, shellification, source selection, terminal
lowering or realization, RG/MWG/IDA, Hamiltonian semantics, solvers, ECP, or
Cr2 workflow policy.
