# Raw Product Source-Mode Inventory Cleanup

Status: approved for implementation under `HP-RAW-SRCMODE-FN-01` and
`HP-RAW-SRCMODE-TEST-01`.

## Purpose

`RawProductBoxPlan` still carries source-mode inventory length in concrete
tuple types. That was deliberately left out of the retained-unit route
inventory cleanup because raw product source planning crosses different
ownership boundaries. This lane removes that route-size type surface while
preserving deterministic source-mode ordering, retained-rule behavior, and
manifest source-mode provenance.

## Approved Boundary

Approved source files:

```text
src/cartesian_raw_product_sources/records.jl
src/cartesian_raw_product_sources/source_mode_indices.jl
src/cartesian_raw_product_sources/summaries.jl
```

Narrow consumer files, only as required by the storage change:

```text
src/cartesian_retained_unit_transform_contracts/unit_contracts.jl
src/cartesian_final_basis_realization/pqs_terminal_basis_realization.jl
src/cartesian_base_hamiltonian.jl
```

Approved changes:

- replace `RawProductBoxPlan.source_mode_indices::Tuple{Vararg{NTuple{3,Int}}}`
  with vector-backed storage;
- replace `RawProductBoxPlan.source_mode_column_indices::Tuple{Vararg{Int}}`
  with vector-backed storage, or remove it if it is exactly `1:count` and
  accessors can supply the same ordered column numbers;
- preserve fixed `NTuple{3,Int}` source-mode coordinates and source-mode
  dimensions;
- preserve deterministic source-mode order and column-number association;
- update summaries and the listed narrow consumers to use the stable accessor
  contract.

Accessor semantics mean ordered mode facts, length, indexing/iteration where
currently used, and retained-rule parity. They do not require preserving a
variable-length `Tuple` return type. Callers must not depend on source-mode
inventory size being encoded in the concrete type.

## Deferred

This lane does not approve:

- terminal-lowering `contracts` or `available_contracts` tuple cleanup;
- retained-unit transform-contract tuple cleanup outside the listed narrow
  consumer wiring;
- broad pair-block or source-box rewrites;
- public input `NamedTuple` changes;
- fixed `NTuple{3,Int}` coordinate/dimension changes;
- artifact schema changes or matrix-key changes.

## Forbidden

This lane does not approve:

- public API/export changes;
- canonical driver changes;
- Hamiltonian object changes or reader changes;
- numerical kernel changes, matrix value changes, or route semantic changes;
- terminal lowering, shellification, terminal basis, Residual Gaussian, raw
  Gaussian block, IDA/MWG, or Qiu-White semantic changes;
- route-stage objects, report/status/payload expansion, persistent caches, or
  compatibility layers that preserve the old tuple-backed shape;
- new committed tests, tools, benchmarks, input fixtures, Cr2 runs, or
  Cr2-specific workflow.

## Validation

Required validation:

- `git diff --check`;
- package load;
- H2 base artifact write/readback;
- H2 supplemented artifact write/readback;
- H2 R3 endpoint;
- focused raw-product source order and retained-rule parity;
- manifest source-mode and final-basis source-relation inspection;
- focused search confirming `RawProductBoxPlan` no longer stores
  source-mode inventories as `Tuple{Vararg{...}}`;
- no Cr2 run.

Existing committed tests may be adjusted only if they directly assert the old
tuple-backed source-mode inventory shape. No new committed test file is
approved.

Line budget: at most `150` added `src` lines, with net simplification expected.

Failure rule: if vectorizing the raw product plan requires source files outside
the approved surfaces, broad pair-block/source-box rewrites, public API or
artifact changes, numerical changes, or compatibility layers preserving the old
tuple-backed shape, make no source commit and report the exact caller/blocker.
