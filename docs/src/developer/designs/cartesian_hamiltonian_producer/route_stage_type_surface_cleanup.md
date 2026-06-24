# Route/Stage Type-Surface Cleanup

Status: approved for implementation under `HP-ROUTE-STAGE-TYPE-FN-01` and
`HP-ROUTE-STAGE-TYPE-TEST-01`.

## Purpose

The Be2 q5 p10 supplemented driver path is fast when warm, but cold
construction remains dominated by compile/type latency. Trace attribution on
`b17b9161` points to oversized route/stage type surfaces rather than package
load, artifact writing, Gaussian raw blocks, or numerical kernels. This lane
allows a narrow cleanup of stale route/stage compatibility inventories and
runtime-sized carriers in the attributed surfaces.

## Approved Boundary

Approved source files:

```text
src/pqs_source_box_route_driver_helpers.jl
src/cartesian_terminal_shellification_geometry.jl
```

Approved targets:

- `_pqs_source_box_route_driver_terminal_lowering_contract_inventory_from_plan`;
- `cartesian_units` route/stage return surfaces that carry oversized
  compatibility inventories;
- `_pqs_source_box_route_driver_transform_stage_low_order_summary`;
- `cartesian_transforms` route/stage return surfaces that carry oversized
  compatibility inventories;
- `_cartesian_terminal_shellification_region_unit_inventory`;
- related terminal-region lowering inventory summary surfaces in
  `src/cartesian_terminal_shellification_geometry.jl` only where the same
  runtime-sized type-surface pattern appears.

Approved changes:

- delete stale route/stage compatibility inventories that no active approved
  caller needs;
- replace remaining runtime-sized `NamedTuple` / `Tuple` carriers with
  vector-backed compact internal objects, stable dictionaries, accessors, or
  smaller summaries;
- shrink wide internal stage return signatures only where all live approved
  callers can be updated within the approved source files;
- preserve deterministic terminal shellification/lowering order and existing
  behavior.

This lane is about compile/type-surface cleanup. It is not authority to change
route semantics, numerical construction, artifacts, or user-visible driver
behavior.

## Required Preservation

The source pass must preserve:

- H2 base artifact/readback behavior;
- H2 supplemented artifact/readback behavior;
- deterministic terminal shellification/lowering order;
- existing public driver contract;
- existing artifact schema and manifest behavior;
- existing numerical matrices.

## Forbidden

This lane does not approve:

- source files outside the two approved files;
- driver changes;
- artifact schema or manifest changes;
- public API/export changes;
- numerical kernel changes or matrix value changes;
- raw-block changes;
- Residual Gaussian, MWG, or IDA semantic changes;
- route semantic changes or shellification behavior changes;
- route diagnostic/status/report expansion;
- broad route-stage redesign;
- new public contracts;
- PackageCompiler, PrecompileTools, sysimage, or precompile workload work;
- new committed tests, tools, benchmarks, input fixtures, Cr2 runs, or
  Cr2-specific workflow.

No compatibility adapter may preserve the old runtime-sized type surface merely
under a new name.

## Validation

Required validation:

- `git diff --check`;
- package load;
- H2 base artifact write/readback;
- H2 supplemented artifact write/readback;
- H2 R3 endpoint if the pass touches terminal realization behavior;
- focused terminal shellification/lowering order parity;
- focused scan for newly introduced `NamedTuple{...}`, variable-size
  `Tuple(...)`, `Tuple{Vararg{...}}`, and runtime-keyed inventories in the
  approved files;
- no Cr2 run.

Optional after correctness passes:

- Be2 q5 compile/timing comparison using the same explicit p10-style
  supplemented driver fixture that produced the attribution.

Existing committed tests may be adjusted only if they directly assert the old
stale compatibility inventory shape. No new committed test file is approved.

Line budget: at most `200` added `src` lines, with net simplification expected.

Failure rule: if cleanup requires source files outside the approved boundary,
broad route-stage redesign, new public contracts, artifact changes, numerical
changes, or a precompile/sysimage mechanism, make no source commit and report
the exact blocker.
