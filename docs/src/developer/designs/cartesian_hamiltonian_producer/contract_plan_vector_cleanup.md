# Contract-Plan Vector Cleanup

Status: approved for implementation under `HP-CONTRACT-VEC-FN-01` and
`HP-CONTRACT-VEC-TEST-01`.

## Purpose

After retained-unit route inventories and raw product source-mode inventories,
the next type-surface cleanup target is plan-level contract inventories that
cross stage boundaries as variable-length tuple fields. This lane replaces
those plan inventories with vector-backed storage while preserving the existing
accessor contracts, iteration order, and behavior.

## Approved Boundary

Approved source files:

```text
src/cartesian_terminal_lowering/contracts.jl
src/cartesian_terminal_lowering/selection.jl
src/cartesian_terminal_lowering/summaries.jl
src/cartesian_retained_unit_transform_contracts/records.jl
src/cartesian_retained_unit_transform_contracts/unit_contracts.jl
src/cartesian_retained_unit_transform_contracts/summaries.jl
```

Narrow consumer files, only as required by the storage change:

```text
src/pqs_source_box_route_driver_helpers.jl
src/cartesian_base_hamiltonian.jl
src/cartesian_final_basis_realization/pqs_terminal_basis_realization.jl
```

Approved changes:

- replace
  `TerminalLoweringPlan.available_contracts::Tuple{Vararg{TerminalLoweringContract}}`
  with vector-backed storage;
- replace
  `TerminalLoweringPlan.contracts::Tuple{Vararg{TerminalLoweringContract}}`
  with vector-backed storage;
- replace
  `RetainedUnitTransformContractPlan.contracts::Tuple{Vararg{RetainedUnitTransformContract}}`
  with vector-backed storage;
- preserve the accessor surface:
  `available_contracts(plan)`, `selected_contracts(plan)`, `contracts(plan)`,
  and `transform_contracts(plan)`;
- preserve iteration order, selected-contract semantics, transform-contract
  semantics, summaries, and existing behavior.

Accessor compatibility means the same ordered contract facts and iteration
semantics. It does not require preserving variable-length `Tuple` concrete
field types or accessor return types.

## Explicitly Out Of Scope

This lane does not approve changing:

- `source_cpbs::Tuple{Vararg{CoordinateProductBox}}`;
- raw product source-mode storage;
- retained-unit route inventories;
- public input `NamedTuple`s;
- fixed coordinate/product-box value objects whose tuple shape is small and
  mathematically meaningful.

`source_cpbs` is a smaller per-contract tuple and is not the primary
stage-boundary plan inventory targeted by this pass.

## Forbidden

This lane does not approve:

- source files outside the approved surfaces;
- public API/export changes;
- canonical driver changes;
- Hamiltonian object changes, matrix-key changes, reader changes, or
  artifact/manifest schema changes;
- numerical kernel changes, matrix value changes, route semantic changes, or
  shellification behavior changes;
- raw product source-mode changes, raw Gaussian block changes, Residual
  Gaussian changes, IDA/MWG changes, or Qiu-White semantic changes;
- route-stage objects, report/status/payload expansion, persistent caches, or
  compatibility adapters preserving the old tuple-backed plan field types;
- new committed tests, tools, benchmarks, input fixtures, Cr2 runs, or
  Cr2-specific workflow.

## Validation

Required validation:

- `git diff --check`;
- package load;
- H2 base artifact write/readback;
- H2 supplemented artifact write/readback;
- H2 R3 endpoint;
- focused terminal-lowering contract order parity;
- focused retained-unit transform-contract order parity;
- focused search confirming the targeted plan inventories no longer store
  contracts as `Tuple{Vararg{...}}`;
- no Cr2 run.

Existing committed tests may be adjusted only if they directly assert the old
tuple-backed contract-plan field shape. No new committed test file is approved.

Line budget: at most `150` added `src` lines, with net simplification expected.

Failure rule: if vectorizing the plan inventories requires source files outside
the approved surfaces, broad route/stage rewrites, public API or artifact
changes, numerical changes, or compatibility layers preserving the old
tuple-backed plan field types, make no source commit and report the exact
caller/blocker.
