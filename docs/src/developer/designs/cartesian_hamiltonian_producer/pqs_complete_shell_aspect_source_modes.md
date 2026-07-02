# PQS Complete-Shell Aspect-Aware Source Modes

Status: design and future source-policy authority under
`HP-PQS-ASPECTSHELL-FN-01` and `HP-PQS-ASPECTSHELL-TEST-01`. This document
approves no source implementation in this pass.

## Purpose

The terminal due-diligence report may warn when a rectangular physical
complete shell is represented by cubic `(q, q, q)` source modes. Reporting the
mismatch is not a basis-policy change. Fixing it is: it changes retained
counts, final dimensions, transforms, Hamiltonian matrices, and downstream
energies.

This design records the separate source-policy lane for restoring explicit
aspect-aware PQS complete-shell source dimensions.

## Problem

Current PQS terminal lowering hard-codes complete-shell source modes as cubic:

```julia
source_mode_shape = ntuple(_ -> policy.q, 3)
```

in `src/cartesian_terminal_lowering/region_contracts.jl`'s
`_pqs_complete_shell_contract(...)`.

For a z-axis diatomic, a shared complete shell can be physically rectangular.
The H2+ due-diligence audit showed a shell with physical side lengths roughly:

```text
x/y/z = 3.464 x 3.464 x 6.646
```

but current PQS lowering used:

```text
actual source_mode_shape = (5, 5, 5)
```

The diagnostic aspect-balanced estimate was:

```text
expected source_mode_shape = (5, 5, 10)
retained count: 98 -> 178
```

That warning is only a diagnostic today. The source-policy fix is to make the
PQS complete-shell source mode shape explicitly aspect-aware.

## Old Code To Recover

The repo still contains the older angular-resolution machinery that separated
transverse `q` from the bond-axis source length `L`.

Important existing primitives:

- `src/cartesian_nested_diatomic.jl`
  - `_nested_diatomic_reference_band(...)`
  - `_nested_diatomic_shared_shell_reference_band(...)`
  - `_nested_diatomic_choose_shell_axis_retain_count(...)`
  - `_nested_diatomic_adaptive_shell_retention(...)`
  - `_nested_diatomic_source_box_dimension_plan(...)`
  - `_nested_diatomic_projected_q_shell_adaptive_source_dimensions(...)`
  - `_nested_projected_q_shell_source_mode_plan(...)`

- `src/cartesian_nested_faces.jl`
  - `_nested_projected_q_shell_layer(...)`

The important old concepts are:

```text
selected_q:        requested transverse PQS source size
raw_source_dims:   actual source-mode dimensions, e.g. (q, q, L)
raw_q:             transverse raw source dimension
raw_L:             bond-axis raw source dimension
```

The older adaptive path did not choose `L` by arbitrary shell depth. It chose
axis retained counts by comparing angular coverage statistics to a reference
angular-resolution band:

```text
reference band -> per-axis retained count -> source_mode_dims
```

Then it fed:

```text
raw_source_dims
selected_q
raw_q
raw_L
axis_selector_retained_counts
```

into `_nested_projected_q_shell_layer(...)`. That layer already accepts
non-cubic `raw_source_dims` and records `raw_q_matches_selected_q`,
`source_mode_dims`, and `raw_L`.

There is also a simpler aspect-aware rule already used for central distorted
product boxes in `src/cartesian_shellification/terminal_geometry.jl`:

```julia
aspect_ratio = physical_sizes[axis] / transverse_size
L = max(shell_side, round(Int, shell_side * aspect_ratio))
source_mode_shape = ntuple(a -> a == axis ? L : shell_side, 3)
```

That central-box rule is useful as a diagnostic cross-check, but the preferred
complete-shell repair should restore or re-express the older angular-band
selection explicitly.

## Current Blockers

The source-policy blocker is narrow but real:

- `_pqs_complete_shell_contract(...)` currently records `(q, q, q)`;
- `src/pqs_multilayer_shell_source_plan.jl` currently rejects non-cubic
  `raw_source_dims` and passes `L = q`;
- `src/pqs_multilayer_shell_region_plan.jl` can carry `source_mode_shape`
  from terminal lowering into shell layers, but must be validated for
  non-cubic source shapes;
- raw-product source-mode indexing and retained rules already support general
  three-axis source dimensions, but need focused validation in this path;
- terminal realization validates source-mode shape consistency and should
  continue doing so.

## Approved Source Policy

For z-axis diatomic PQS complete shells, source-mode shape should match the
physical angular-resolution requirement of the shell.

The intended shape is:

```text
source_mode_shape = (q, q, L)
```

for bond axis `z`, with the analogous axis permutation for `x` or `y` if a
future approved geometry uses another bond axis.

Policy:

- `q` remains the selected transverse PQS source size.
- `L` is a bond-axis source length chosen by an explicit angular-resolution
  rule, not by a hidden cubic default.
- The first implementation should restore or re-express the older
  angular-band retained-count logic:
  - build a reference angular-resolution band;
  - evaluate candidate retained counts along each axis;
  - choose the count whose angular coverage matches the reference band;
  - convert retained counts to source-mode dimensions with the established
    source-mode padding convention;
  - pass `raw_source_dims`, `selected_q`, `raw_q`, and `raw_L` explicitly into
    the projected-q shell layer.
- A simple physical-aspect estimate such as
  `L = max(q, round(Int, q * aspect_ratio))` may be used only as a diagnostic
  or fallback candidate if the source pass explicitly validates that it matches
  the angular-band rule for the approved fixtures.

## Registered IDs

### HP-PQS-ASPECTSHELL-FN-01

Approved future source surface:

```text
src/cartesian_terminal_lowering/region_contracts.jl
src/pqs_multilayer_shell_source_plan.jl
src/pqs_multilayer_shell_region_plan.jl
```

Optional only if directly needed to reuse the old angular-resolution helpers:

```text
src/cartesian_nested_diatomic.jl
src/cartesian_nested_faces.jl
```

Approved behavior:

- change PQS complete-shell source-mode shape from cubic `(q, q, q)` to an
  explicit aspect-aware `(q, q, L)` shape for z-axis diatomic complete shells;
- derive `L` from the restored angular-resolution rule, or from a documented
  validated equivalent;
- pass non-cubic `raw_source_dims` through the multilayer shell source plan;
- call `_nested_projected_q_shell_layer(...)` with explicit `q`, `L`,
  `raw_source_dims`, and `selected_q`;
- preserve shell support ownership and shell-local projection/Lowdin cleanup;
- preserve due-diligence reporting of actual and expected source-mode shapes.

Forbidden:

- artifact schema/provenance/reader changes;
- public input or driver semantic changes;
- WL source-mode or retained-basis policy changes;
- thin-slab, angular z-extension, direct/core identity, or residual/MWG/IDA
  changes;
- global residual/injection changes;
- old route-global materialization revival;
- broad source-mode framework or report/payload expansion;
- Cr2 production claim.

Expected consequences:

- retained counts and final dimensions may change;
- Hamiltonian matrices and energies may change;
- old scalar targets tied to cubic complete-shell source modes must be
  remeasured rather than preserved.

Failure rule: if non-cubic complete-shell source modes require changing
support ownership, terminal realization semantics, artifact schema, public
driver inputs, or a broad route/report framework, make no source commit and
report the blocker.

Line budget: target at most `160` added `src` lines. If restoring the angular
selection needs a broader extraction from the old diatomic high-order path,
stop and request a narrower helper-authority amendment.

### HP-PQS-ASPECTSHELL-TEST-01

Approved validation:

- `git diff --check`;
- package load;
- focused complete-shell source-shape probe showing a rectangular physical
  shell uses `(q, q, L)` rather than `(q, q, q)`;
- retained count matches `prod(source_mode_shape) -
  prod(source_mode_shape .- 2)` for the selected shell;
- due-diligence report shows actual shape, expected aspect shape, retained
  count, and no stale cubic-shape warning for the repaired shell;
- bounded H2 or H2+ artifact/readback smoke;
- finite/symmetric base Hamiltonian matrices if an artifact is written;
- no Cr2 run required.

No committed fixtures or tests are approved by default.

## Boundary With Due Diligence

`HP-DRV-SHELLDD-*` may report:

```text
actual source_mode_shape
expected aspect-balanced source_mode_shape
retained-count delta
warning flags
```

It must not change construction behavior.

`HP-PQS-ASPECTSHELL-*` is the source-policy lane that may change construction
behavior by making PQS complete shells aspect-aware.
