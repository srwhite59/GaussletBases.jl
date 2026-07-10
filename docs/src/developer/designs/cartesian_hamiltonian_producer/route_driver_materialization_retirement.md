# Route-Driver Materialization Workflow Retirement

Status: completed by `e2e164e9b`, with dangling ladder runners removed by
`77fa2700b`. The `HP-RETIRE-DRV-MAT-*` and
`HP-RETIRE-LADDER-RUNNERS-*` IDs are closed historical deletion authority.

## Decision

The old route-driver materialization/report/save wrapper workflow was retired.
The canonical Cartesian producer is the staged human-facing driver plus the
`CartesianIDAHamiltonian` artifact path; CR2-facing artifact production does
not depend on route-stage wrapper choreography.

## Removed Surface

Commit `e2e164e9b` removed the old wrappers and their private helpers:

```text
cartesian_materialization
cartesian_print_summary
cartesian_print_details
cartesian_save
_pqs_source_box_route_driver_materialization
_pqs_source_box_route_driver_print_materialization
_pqs_source_box_route_driver_save
```

It also removed the obsolete reporting module, old harness/probe/smoke tools,
and stale docs-policy assertions. Commit `77fa2700b` then removed:

```text
tools/run_cartesian_driver_ladder.jl
tools/run_cartesian_line_ladder.jl
```

`tools/cartesian_driver_ladder_lib.jl` remains quarantined historical support,
not a canonical workflow or source of driver authority.

## Current Guardrail

Current focused scans find no live wrapper or deleted-runner reference in
`src`, `bin`, `test`, or `tools`. Do not restore wrappers, runner entrypoints,
replacement adapters, status/payload objects, or their behavior in the
canonical driver. A new docs-only amendment is required before related source,
test, tool, or workflow work.
