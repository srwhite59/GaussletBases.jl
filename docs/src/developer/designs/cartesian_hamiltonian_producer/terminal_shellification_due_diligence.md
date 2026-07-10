# Terminal Shellification Due Diligence

Status: implemented producer/driver reporting contract under
`HP-DRV-SHELLDD-FN-01` and `HP-DRV-SHELLDD-TEST-01`. It changes no numerical
construction or artifact schema.

## Purpose

Basis due diligence is part of the Cartesian/PQS producer contract. The repo
exposes a standard compact due-diligence report from driver/producer
workflows so repo consumers can confirm the derived system, parent axes,
weights, dimensions, and shellification before interpreting energies,
residual occupations, injection behavior, or Cr2-style failure modes.

The old legacy driver practice was useful because it printed enough
construction detail that a human could see suspicious basis geometry. The new
canonical driver should recover that property through structured bounded
report sections rather than a noisy route-stage dump.

## Why This Is Needed

Condensed region summaries can hide basis-construction problems. In the H2+
`ns = 5`/`ns = 7` audit, a due-diligence probe caught the suspicious
`complete_shell_1` immediately:

```text
physical side lengths:          3.464 x 3.464 x 6.646
source_mode_shape:              (5, 5, 5)
expected aspect-balanced shape: (5, 5, 10)
retained_count:                 98
expected retained scale:        178
```

This was not a normal code bug. It was inadequate basis construction: a
rectangular physical shell was represented by cubic `q x q x q` source modes.
That kind of problem should be visible before a consumer reasons about
energies or residual/injection behavior.

## Contract

- The repo/driver must expose a terminal basis due-diligence report for
  Cartesian/PQS terminal bases.
- Consumers are expected to inspect this report before interpreting energies,
  residual-Gaussian behavior, injection behavior, or high-cost production
  artifacts.
- Warning flags are advisory diagnostics. They are not automatic construction
  failures unless a caller, test, or later policy explicitly chooses to enforce
  them.
- The report is a user-facing basis review surface, not a route-stage
  diagnostic dump. It must remain bounded and mostly row-oriented.

## Required Report Sections

The due-diligence report should have four compact sections:

1. normalized system and geometry context;
2. parent axes, physical box, 1D centers, and gausslet/IDA weight statistics;
3. final-basis dimension and compression accounting;
4. shell-by-shell terminal region table.

The implemented report may print compact summaries with wrapped rows. It also
builds an in-memory row/section representation so callers can inspect the same
facts without scraping text.

## System And Geometry Header

The header should report derived, normalized facts, not just echo user input:

```text
geometry_kind
nesting
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
```

For one-center atoms, the report should make it clear whether public
padding/radius actually changed parent counts and physical box size. For
z-axis diatomics, the report should make the bond length, bond axis, and
transverse versus longitudinal extents visible.

## Parent Axis And Weight Tables

For each parent axis `x`, `y`, and `z`, report:

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
core_region_index_span
tail_region_index_span
```

The full 1D center locations are useful derived facts. The normal driver
printout may use a bounded preview when an axis is long, but the in-memory
report row should keep enough structured data to emit a full per-axis table in
a focused report without changing artifact schema.

Gausslet/IDA weight statistics should be reported where available:

```text
axis or total scope
weight_count
weight_sum
weight_min
weight_max
weight_abs_sum
negative_weight_count
near_zero_weight_count
near_zero_threshold
large_weight_warning
```

If meaningful, include per-axis 1D weight stats and aggregate 3D/base IDA
weight stats. These are diagnostics only. They are not residual integral
weights, not MWG weights, and not automatic proof of quadrature quality.

## Dimension And Compression Accounting

Report a compact dimension summary:

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

This section should let a consumer see whether the basis size is coming from
core/contact sectors, complete shells, slabs, residuals, or supplement-derived
columns before inspecting detailed shell rows.

## Required Shell Table Fields

The due-diligence table should include one row per terminal region or shell
unit at the granularity needed to review shellification. The implemented
report uses compact rows and wrapped printing, while the in-memory row carries
these facts where available:

```text
terminal_order
terminal_key
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
retained_count
final_column_range
lowering_kind
retained_rule
realization_rule_or_status
slab normal axis / side / thickness / stack index / stack count
warning_flags
warning_summary
```

Field meanings:

- `terminal_order` and `terminal_key` identify the row in matrix/final-basis
  order.
- `role`, `region_kind`, and `shell_index` identify shellification role
  without requiring route-internal reconstruction.
- `owner/contact/shared classification` distinguishes atom-local/contact
  sectors from shared molecular shells and slab-like regions.
- Index boxes/shapes report the discrete parent support.
- Physical bounds and side lengths report the real geometry represented by the
  index support.
- Physical aspect ratios should make rectangular shells visible. For z-axis
  diatomics, the transverse-to-longitudinal relationship must be visible.
- `source_mode_shape` is the actual source-mode shape used by the basis
  construction.
- `expected_aspect_balanced_source_mode_shape` is an advisory diagnostic shape
  for shell-like regions where physical aspect ratio makes a cubic source
  shape suspicious. This field is reporting only; it does not implement
  aspect-balanced source modes.
- `source_mode_count`, `retained_count`, and `final_column_range` connect the
  shellification row to final-basis size.
- `lowering_kind`, `retained_rule`, and `realization_rule_or_status` state how
  support became final columns.
- Slab metadata is required for midpoint, angular z-extension, boundary, and
  fallback slab rows when present.
- `warning_flags` should be short stable symbols; `warning_summary` should be
  concise human text.

## Initial Warning Flags

The first reporting pass should support advisory warning flags such as:

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

This list is not a new enforcement policy. It is a bounded diagnostic
vocabulary so humans and repo consumers can see suspicious construction facts
consistently.

## Implemented Seam

The implementation extends
`src/cartesian_base_hamiltonian.jl`'s
`_cartesian_terminal_inventory_rows(...)`.

Current shape:

- joins existing terminal inventory rows with terminal retained-rule
  plan/support records;
- gathers normalized system/geometry context and parent-axis summaries from
  existing staged producer objects;
- gathers gausslet/IDA weight statistics only from existing weights already
  present in the construction path;
- produces an in-memory/report object;
- has the canonical driver print the bounded due-diligence report through the
  existing driver summary path;
- keeps the report compact enough for normal driver output, with bounded
  previews for long axis-center lists;
- does not change artifact schema.

Implemented source surface:

```text
src/cartesian_base_hamiltonian.jl
bin/cartesian_ham_builder.jl
```

Optional only if a compact accessor is directly required:

```text
src/pqs_source_box_route_driver_helpers.jl
```

No new module, route-report framework, persistent payload object, or artifact
sidecar group is approved by this design.

## Registered IDs

### HP-DRV-SHELLDD-FN-01

Approved behavior:

- add one helper/report surface for terminal due diligence;
- expose the report from canonical driver/producer workflows;
- include the required system, axis/box/weight, dimension-accounting, and
  shell-by-shell fields listed above where available;
- compute advisory warning flags, including rectangular physical shells
  represented by cubic source modes and anomalous derived axis/weight facts;
- keep warning flags advisory by default;
- preserve existing compact terminal-region inventory behavior unless it is
  intentionally extended by this report;
- preserve all numerical construction behavior.

Forbidden:

- artifact schema/provenance/reader changes;
- public input or driver semantic changes;
- shellification policy changes;
- source-mode selection changes;
- aspect-balanced source-mode implementation;
- terminal lowering, retained-unit, terminal-realization, RG/MWG/IDA,
  Hamiltonian, raw-block, solver, or Cr2 workflow changes;
- route skeleton exposure, source-mode inventory dumps, pair inventories,
  raw-block details, all-row support listings, full metadata dumps, or
  recursive route-stage reports;
- dense coefficient, transform, pair, or raw support dumps;
- automatic failure on warning flags unless a later policy approves it.

Ongoing guardrail: maintenance must keep using
`_cartesian_terminal_inventory_rows(...)` and compact accessors. If a requested
extension would require a broad report/payload framework, artifact fields, or
shellification policy changes, stop and request new authority.

Historical implementation budget: the accepted lane targeted at most `180`
added `src`/`bin` lines for one report/table surface, not a reporting subsystem.

### HP-DRV-SHELLDD-TEST-01

Implemented validation contract:

- `git diff --check`;
- package load if source is touched;
- bounded H2 or H2+ driver/producer smoke showing due-diligence report
  sections and shell rows;
- focused row inspection showing a rectangular physical shell warning when an
  existing bounded fixture has one;
- focused inspection of normalized system/geometry, axis/box summaries, and
  gausslet/IDA weight statistics;
- confirm ordinary compact terminal inventory output remains bounded;
- confirm artifact/readback matrix deltas are unchanged if artifact writing is
  exercised;
- no Cr2 run required.

No committed fixtures or tests are approved by default.

## Open Follow-Up

Aspect-balanced complete-shell source modes are a separate source-policy fix
under `HP-PQS-ASPECTSHELL-FN-01` and `HP-PQS-ASPECTSHELL-TEST-01`. This
due-diligence lane must not implement that policy. It should make the problem
visible by reporting actual source-mode shape, expected aspect-balanced shape,
retained count, final columns, and warning flags.
