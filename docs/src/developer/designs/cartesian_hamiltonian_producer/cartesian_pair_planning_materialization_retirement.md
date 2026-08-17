# Cartesian Pair Planning And Materialization Retirement

## Status And Authority

This page records the completed retirement of:

- `HP-RETIRE-PAIR-LADDER-FN-01`;
- `HP-RETIRE-PAIR-LADDER-TEST-01`.

Status: source retirement completed at `32bb3e2a7` and was accepted in Pass
470. The function ID is retired and the validation ID is completed; neither
grants further work.

## Decision

Retire the unused development-era ladder:

```text
retained units
-> upper-triangular unit-pair inventory
-> pair-operator plans
-> pair-block preflight and local materialization pilots
```

The current source-box-first PQS and White--Lindsey producer does not consume
this ladder. It constructs terminal blocks from the common shellification and
lowering records, then builds exact one-body and IDA operators through the
current terminal owners. Keeping a second metadata/materialization story adds
source volume and misleading architecture without protecting a live endpoint.

The qualified internal APIs are not public compatibility obligations. Do not
preserve them through aliases, stubs, deprecations, adapters, moved helpers,
or replacement tests.

## Retired Source

Pass 469 deleted these complete source families:

```text
src/cartesian_unit_pairs/
src/cartesian_pair_operator_plans/
src/cartesian_route_core/unit_pairs.jl
src/cartesian_route_core/pair_operator_plans.jl
```

In `src/cartesian_pair_block_materialization/`, it deleted every file except:

```text
pqs_source_axis_transforms.jl
```

It also deleted the orphaned final-basis prototype and standalone wrappers:

```text
src/cartesian_final_basis_realization/pqs_source_shell_final_basis.jl
modules/CartesianCPB.jl
modules/CartesianRouteCore.jl
tools/cartesian_driver_ladder_lib.jl
```

It removed only the matching includes, exports, aliases, imports, and obsolete
module descriptions from:

```text
src/GaussletBases.jl
src/cartesian_route_core/CartesianRouteCore.jl
src/cartesian_final_basis_realization/CartesianFinalBasisRealization.jl
src/cartesian_pair_block_materialization/CartesianPairBlockMaterialization.jl
```

`CartesianPairBlockMaterialization` remains loaded only as the narrow internal
owner of the retained mapped-COMX axis-transform helper. Its wrapper was
reduced from `247` to `19` lines containing only the imports, constant aliases,
export, and include required by that helper. This is retained numerical
ownership, not compatibility glue for the retired ladder.

## Mapped-COMX Exception

The sole live numerical exception is:

```text
src/cartesian_pair_block_materialization/pqs_source_axis_transforms.jl
pqs_source_axis_transform_facts_from_pgdg_axes
```

`src/pqs_source_box_route_driver_helpers.jl` calls this function when it builds
mapped-COMX source-axis facts. `HP-MCOMX-FILE-01` and `HP-MCOMX-WIRE-01`
continue to own the file and function. The retirement must not move, duplicate,
rename, reinterpret, or numerically modify them.

The reduced module must retain the existing doside/PGDG source transform,
mapped source-span selection, `AxisSourceTransformFact` output, overlap checks,
and ordinary fallback behavior exactly.

## Caller Proof

At baseline `140ac6a8f1346426d1cca4aa364ebfcf4cebbf90`:

- no committed source, test, driver, tool, example, or standalone module calls
  the unit-pair, pair-plan, pair-materialization, direct-block, safe-term,
  source-shell bridge, or readiness APIs outside the retiring cluster;
- the source-shell final-basis functions are referenced only by their defining
  module and CPBM aliases;
- the two standalone `modules/` wrappers have no committed consumer;
- the one-line ladder tool has no caller;
- the retiring code owns no artifact, serialization, registration, `__init__`,
  or load-time side effect beyond its root includes;
- ignored historical probes may still name the old APIs, but they are not
  compatibility obligations.

If implementation review finds a live committed caller or side effect, stop
without committing and report it. Do not narrow the deletion by adding glue.

## Accepted Accounting

The audited file content is:

```text
pair-ladder inventory, including current CPBM wrapper/helper    9,431
retained mapped-COMX helper                                     -173
adjacent orphaned final-basis/module/tool files                  +664
direct old file content to delete or replace                    9,922
```

The accepted source delta was `+8/-9,950`, net `-9,942` source lines. Including
the obsolete tool, the complete pass was `+8/-9,951`. The retained
`pqs_source_axis_transforms.jl` remained byte-for-byte unchanged. No source,
test, tool, module, helper, metadata, status vocabulary, artifact, or
compatibility file was added.

## Preserved Owners

The retirement must not delete, move, or reinterpret:

- `src/cartesian_cpb/` geometry;
- current shellification, terminal lowering, retained units, and retained-unit
  transform contracts;
- current terminal PQS and White--Lindsey realization;
- mapped-COMX construction and source-axis facts;
- raw product sources and Gaussian raw blocks;
- exact terminal one-body operators and IDA;
- Residual Gaussians, MWG, PRFs, represented Hartree, artifacts, drivers, or
  solvers.

Current route documentation should describe the live numerical path:

```text
parent and common shellification
-> terminal lowering and retained contracts
-> PQS or White--Lindsey terminal realization
-> exact one-body and IDA assembly
-> optional residual/supplement composition
```

It must not present the retired pair ladder as the intended future architecture.

## Accepted Validation

The source pass:

1. confirmed all retired includes, modules, types, functions, and tool paths
   absent from committed source;
2. passed package load and the unchanged `core`, `ida`, and `cartesian` groups
   at `12,720/12,720`;
3. passed the public Cartesian endpoint at `232/232`, with PQS/WL dimensions
   `487/487`, parent axes `9 x 9 x 15`, and no padding or weight warnings;
4. reproduced the mapped-COMX construction exactly: source shape `(5,5,5)`,
   `125` source modes, terminal shape `1331 x 98`, and unchanged ordinary,
   mapped, and terminal-coefficient fingerprints;
5. preserved the retained helper path and file bytes exactly;
6. passed the deleted-symbol, anti-bloat, and `git diff --check` gates;
7. passed GitHub numerical CI. The first Docs run failed only because machine
   authority still named the deleted paths; Pass 470 closes that lifecycle.

The existing tests are validation inputs, not edit surfaces. Add no committed
probe or replacement test for deleted internal APIs.

## Non-Goals

This closed authority does not permit:

- a new pair-planning, pair-materialization, provider, or route framework;
- changes to terminal coefficients, mapped-COMX output, one-body operators,
  IDA/MWG, residual selection, or public behavior;
- deletion of surviving CPB, shellification, lowering, retained-unit,
  transform-contract, raw-product, terminal, or Hamiltonian owners;
- changes to artifacts, metadata schemas, reports, drivers, public APIs,
  numerical defaults, or physics policy;
- deletion of `CartesianParentAxisFactors.jl` or any broader QW, protected
  ladder, chain/square, or large-file surface.

Reintroduction of any retired API or implementation requires a new docs-only
authority amendment tied to a current numerical consumer.
