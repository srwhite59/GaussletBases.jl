# Cartesian Pair Planning And Materialization Retirement

## Status And Authority

This page owns the approved retirement contract for:

- `HP-RETIRE-PAIR-LADDER-FN-01`;
- `HP-RETIRE-PAIR-LADDER-TEST-01`.

Status: approved for exact source and tool retirement in Pass 468. The source
pass is separate. Until it lands, the old modules remain loaded but grant no
new implementation direction.

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

## Exact Retirement Surface

Delete these complete source families:

```text
src/cartesian_unit_pairs/
src/cartesian_pair_operator_plans/
src/cartesian_route_core/unit_pairs.jl
src/cartesian_route_core/pair_operator_plans.jl
```

In `src/cartesian_pair_block_materialization/`, delete every file except:

```text
pqs_source_axis_transforms.jl
```

Delete the orphaned final-basis prototype and standalone wrappers:

```text
src/cartesian_final_basis_realization/pqs_source_shell_final_basis.jl
modules/CartesianCPB.jl
modules/CartesianRouteCore.jl
tools/cartesian_driver_ladder_lib.jl
```

Remove only the matching includes, exports, aliases, imports, and obsolete
module descriptions from:

```text
src/GaussletBases.jl
src/cartesian_route_core/CartesianRouteCore.jl
src/cartesian_final_basis_realization/CartesianFinalBasisRealization.jl
src/cartesian_pair_block_materialization/CartesianPairBlockMaterialization.jl
```

`CartesianPairBlockMaterialization` remains loaded only as the narrow internal
owner of the retained mapped-COMX axis-transform helper. Reduce its wrapper to
the imports, constant aliases, export, and include required by that helper.
This is retained numerical ownership, not compatibility glue for the retired
ladder.

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

## Accounting And Budget

The audited file content is:

```text
pair-ladder inventory, including current CPBM wrapper/helper    9,431
retained mapped-COMX helper                                     -173
adjacent orphaned final-basis/module/tool files                  +664
direct old file content to delete or replace                    9,922
```

The 247-line CPBM wrapper is replaced by a minimal owner around the retained
173-line helper. Root include/export/docstring cleanup adds further deletions.
The source pass must reduce tracked source by at least approximately 9,900
lines.

At most 30 added source lines are allowed for the reduced CPBM owner. Add no
new source, test, tool, module, helper, metadata, status vocabulary, artifact,
or compatibility file. No test edit is authorized.

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

## Validation Contract

The source pass must:

1. confirm all retired includes, modules, types, functions, and tool paths are
   absent from committed source;
2. run package load and the unchanged `core`, `ida`, and `cartesian` groups;
3. inspect the public Cartesian terminal due-diligence endpoint;
4. run one transient mapped-COMX terminal construction before and after the
   deletion and require exact dimension plus coefficient/transform fingerprint
   parity;
5. confirm the retained helper path and mapped-COMX authority paths are
   unchanged;
6. report exact added/deleted line accounting and the anti-bloat diff scan;
7. pass `git diff --check`.

The existing tests are validation inputs, not edit surfaces. Add no committed
probe or replacement test for deleted internal APIs.

## Non-Goals

This authority does not permit:

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
