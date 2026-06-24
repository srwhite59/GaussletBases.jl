# Route Inventory Type-Surface Cleanup

Status: first lane approved for implementation under `HP-ROUTE-INV-FN-01` and
`HP-ROUTE-INV-TEST-01`.

## Purpose

The base/supplemented producer path still carries some route-size facts in
runtime-keyed concrete types. That creates unnecessary compile/type-surface
cost and keeps stale route-inventory vocabulary alive in the live producer
contract. The first cleanup lane removes the clearest case without changing
physics, artifacts, public inputs, or driver behavior.

## First Approved Lane

Approved source file:

```text
src/pqs_source_box_route_driver_helpers.jl
```

Approved targets:

- `_pqs_source_box_route_driver_named_tuple_from_units(...)`;
- runtime-keyed retained-unit inventory fields produced from unit labels,
  including `source_boxes`, `source_dimensions`, `retained_counts`, and
  `ranges`;
- runtime-keyed `pair_family_counts = NamedTuple{families}(...)`;
- internal same-file consumers that currently expect those runtime-keyed
  `NamedTuple` shapes.

Approved replacement shape:

- vector-backed records or tables with stable field names;
- stable dictionaries keyed by unit or family labels where lookup by label is
  genuinely needed;
- helper accessors that hide the storage shape from same-file callers;
- compact summaries that expose counts/order without encoding route size in the
  concrete type.

The retained-unit vector remains the primary ordered inventory. Retained-unit
labels may remain data values, but they must not become type parameters.

## Deferred Lanes

Not approved by `HP-ROUTE-INV-*`:

- `TerminalLoweringPlan.available_contracts::Tuple{Vararg{...}}`;
- `TerminalLoweringPlan.contracts::Tuple{Vararg{...}}`;
- `RetainedUnitTransformContractPlan.contracts::Tuple{Vararg{...}}`;
- public input `NamedTuple`s;
- small fixed `NTuple{3,Int}` coordinate or dimension values;
- artifact sidecar tables.

Raw product source-mode inventory cleanup is separately approved under
`HP-RAW-SRCMODE-FN-01` / `HP-RAW-SRCMODE-TEST-01` because it crosses different
ownership boundaries.
Terminal-lowering and retained-unit transform contract-plan vector cleanup is
separately approved under `HP-CONTRACT-VEC-FN-01` /
`HP-CONTRACT-VEC-TEST-01`.

## Forbidden

This lane does not approve:

- source files outside `src/pqs_source_box_route_driver_helpers.jl`;
- numerical kernel changes;
- route recipe behavior changes;
- terminal lowering, shellification, terminal basis, Residual Gaussian, raw
  product source, or raw-block changes;
- canonical driver changes;
- Hamiltonian object, matrix-key, reader, or artifact schema changes;
- public API/export changes;
- report/status/payload expansion;
- new route-stage objects or compatibility adapters;
- new committed tests or fixtures;
- Cr2 runs or Cr2-specific workflow.

## Validation

Required validation:

- `git diff --check`;
- package load;
- H2 base artifact write/readback through the existing reader;
- H2 supplemented artifact write/readback through the existing reader or
  canonical driver path;
- focused check that no `NamedTuple{unit_keys}` or `NamedTuple{families}` route
  inventory remains in `src/pqs_source_box_route_driver_helpers.jl`;
- no Cr2 run.

Existing tests may be adjusted only if they directly assert the old
runtime-keyed inventory shape. No new committed test file is approved.

Line budget: at most `120` added `src` lines, with net simplification expected.

Failure rule: if this requires source files outside the approved file, broader
route/stage rewiring, public API changes, artifact changes, or preserving the
old runtime-keyed shape through an adapter, make no source commit and report
the blocker.
