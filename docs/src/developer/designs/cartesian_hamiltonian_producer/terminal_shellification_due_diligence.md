# Terminal Shellification Due Diligence

Status: design and future source authority under `HP-DRV-SHELLDD-FN-01` and
`HP-DRV-SHELLDD-TEST-01`. This document approves no production implementation
in this pass and no artifact schema change.

## Purpose

Basis due diligence is part of the Cartesian/PQS producer contract. The repo
should expose a standard shell-by-shell table from driver/producer workflows so
repo consumers can inspect the actual terminal basis construction before
interpreting energies, residual occupations, injection behavior, or Cr2-style
failure modes.

The old legacy driver practice was useful because it printed enough
construction detail that a human could see suspicious basis geometry. The new
canonical driver should recover that property through a structured bounded
table rather than a noisy route-stage dump.

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

- The repo/driver must expose a shell-by-shell basis due-diligence table for
  Cartesian/PQS terminal bases.
- Consumers are expected to inspect this table before interpreting energies,
  residual-Gaussian behavior, injection behavior, or high-cost production
  artifacts.
- Warning flags are advisory diagnostics. They are not automatic construction
  failures unless a caller, test, or later policy explicitly chooses to enforce
  them.
- The table is a user-facing basis review surface, not a route-stage
  diagnostic dump. It must remain bounded and row-oriented.

## Required Table Fields

The due-diligence table should include one row per terminal region or shell
unit at the granularity needed to review shellification. The first
implementation may use compact rows and wrapped printing, but the in-memory
row must carry these facts where available:

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
```

This list is not a new enforcement policy. It is a bounded diagnostic
vocabulary so humans and repo consumers can see suspicious construction facts
consistently.

## Intended Implementation Seam

The first source pass should extend or wrap
`src/cartesian_base_hamiltonian.jl`'s
`_cartesian_terminal_inventory_rows(...)`.

Implementation shape:

- join existing terminal inventory rows with terminal retained-rule
  plan/support records;
- produce an in-memory/report table first;
- have the canonical driver print the bounded due-diligence table through the
  existing driver summary path;
- keep the table row-oriented and compact enough for normal driver output;
- do not change artifact schema in the first implementation.

Approved source surface for the later implementation:

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

- add one helper/table surface for terminal shellification due-diligence rows;
- expose the table from canonical driver/producer workflows;
- include the required shell-by-shell fields listed above where available;
- compute advisory warning flags, including rectangular physical shells
  represented by cubic source modes;
- keep warning flags advisory by default;
- preserve existing compact terminal-region inventory behavior unless it is
  intentionally extended by this table;
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
- automatic failure on warning flags unless a later policy approves it.

Failure rule: if the due-diligence table cannot be built by extending/wrapping
`_cartesian_terminal_inventory_rows(...)` and compact accessors without adding
a broad report/payload framework, artifact fields, or shellification policy
changes, stop and report the missing seam.

Line budget: target at most `120` added `src`/`bin` lines. This should be a
small table/report surface, not a new reporting subsystem.

### HP-DRV-SHELLDD-TEST-01

Approved validation:

- `git diff --check`;
- package load if source is touched;
- bounded H2 or H2+ driver/producer smoke showing due-diligence rows;
- focused row inspection showing a rectangular physical shell warning when an
  existing bounded fixture has one;
- confirm ordinary compact terminal inventory output remains bounded;
- confirm artifact/readback matrix deltas are unchanged if artifact writing is
  exercised;
- no Cr2 run required.

No committed fixtures or tests are approved by default.

## Open Follow-Up

Aspect-balanced complete-shell source modes are likely a separate source-policy
fix. This due-diligence lane must not implement that policy. It should make
the problem visible by reporting actual source-mode shape, expected
aspect-balanced shape, retained count, final columns, and warning flags.
