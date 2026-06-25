# Route-Driver Materialization Workflow Retirement

Status: approved retirement/quarantine authority under
`HP-RETIRE-DRV-MAT-FN-01`, `HP-RETIRE-DRV-MAT-TOOL-01`,
`HP-RETIRE-DRV-MAT-DOC-01`, and `HP-RETIRE-DRV-MAT-TEST-01`.

## Decision

The old route-driver materialization/report/save wrapper workflow is approved
for retirement or quarantine. The current Cartesian producer workflow is the
canonical staged human-facing driver plus the `CartesianIDAHamiltonian`
artifact path. CR2-facing artifact production does not depend on the old
route-stage wrapper choreography.

Retired wrapper names:

```text
cartesian_materialization
cartesian_print_summary
cartesian_print_details
cartesian_save
_pqs_source_box_route_driver_materialization
_pqs_source_box_route_driver_print_materialization
_pqs_source_box_route_driver_save
```

## Evidence

Focused search found no hits for the retired wrapper names in
`bin/cartesian_ham_builder.jl`. The current CR2-facing artifact workflow uses
canonical staged producer functions and `CartesianIDAHamiltonian` artifacts.

Live hits for the audited names are old wrapper definitions, old tools and
harnesses, a stale docs-policy test, and stale compact-doc references. Those
surfaces can pull agents back toward old route-stage choreography and away from
the canonical driver/artifact path.

## Approved IDs

- `HP-RETIRE-DRV-MAT-FN-01` - source wrapper retirement.
- `HP-RETIRE-DRV-MAT-TOOL-01` - old tool/harness quarantine or deletion.
- `HP-RETIRE-DRV-MAT-DOC-01` - active docs/index policy cleanup.
- `HP-RETIRE-DRV-MAT-TEST-01` - validation and stale docs-policy test cleanup.

## Approved Source Surface

Approved later source files:

```text
src/pqs_source_box_route_driver_helpers.jl
src/pqs_source_box_low_order_materialization.jl
src/pqs_source_box_route_driver_reporting.jl
src/GaussletBases.jl
```

`src/GaussletBases.jl` is allowed only if the reporting include becomes unused
after the wrapper/report/save path is removed.

Approved source behavior:

- remove `cartesian_materialization`, `cartesian_print_summary`,
  `cartesian_print_details`, and `cartesian_save`;
- remove implementation helpers used only by those wrappers:
  `_pqs_source_box_route_driver_materialization`,
  `_pqs_source_box_route_driver_print_materialization`, and
  `_pqs_source_box_route_driver_save`;
- remove a related old White-Lindsey atomic pure-gausslet materialization helper
  only if the source pass proves it becomes uncalled;
- do not add replacement wrappers, adapters, compatibility status objects, or
  tests.

## Approved Tool And Test Surface

Approved later tool files:

```text
tools/cartesian_driver_harness.jl
tools/cr2_cartesian_ida_stage_probe.jl
tools/cartesian_driver_ladder_lib.jl
tools/h2_pqs_base_hamiltonian_smoke.jl
```

Approved tool behavior:

- delete or quarantine old tools that exist only to drive the retired wrapper
  workflow;
- do not move diagnostic or ladder behavior into the canonical driver.

Approved later test file:

```text
test/docs/cartesian_ham_builder_policy_runtests.jl
```

Approved test behavior:

- remove or update only old route-stage wrapper assertions that assume the
  canonical driver should avoid calling these now-retired wrappers;
- do not add new committed tests or fixtures.

## Approved Docs Surface

Approved later docs files:

```text
docs/src/developer/algorithm_implementation_index.md
docs/src/developer/designs/cartesian_hamiltonian_producer/implementation_slices.md
docs/src/developer/designs/cartesian_hamiltonian_producer/registry.md
docs/src/developer/designs/cartesian_hamiltonian_producer/r1_public_base_producer.md
docs/src/developer/pqs_manager_running_log.md
```

Approved docs behavior:

- stop describing the old wrapper workflow as canonical or active current
  authority;
- keep historical references historical;
- keep the canonical staged driver and current producer/artifact path unchanged.

## Forbidden

This retirement lane does not approve:

- changes to `bin/cartesian_ham_builder.jl`;
- changes to current staged producer functions;
- artifact schema, provenance, reader, or manifest changes;
- route, shellification, terminal-lowering, raw-block, Residual Gaussian, MWG,
  IDA, or Hamiltonian assembly changes;
- changes to `pqs_multilayer_complete_core_shell_h1.jl`;
- changes to `pqs_complete_core_shell_final_basis.jl`;
- broad ordinary or Qiu-White donor-kernel retirement;
- replacement wrappers, adapters, status fields, payloads, or tests;
- Cr2 workflow.

## Validation

`HP-RETIRE-DRV-MAT-TEST-01` approves only:

- `git diff --check`;
- package load;
- focused `rg` over `src`, `bin`, `test`, and `tools` showing no remaining
  live references to the retired wrapper names;
- canonical small base artifact/readback smoke;
- canonical small supplemented artifact/readback smoke;
- unchanged H2 Residual Gaussian endpoint;
- docs-policy test passing after update or the stale wrapper assertion removed;
- no Cr2 run.

## Failure Rule

If any current canonical producer path or public artifact workflow depends on
these wrappers, make no source commit and report the exact dependency. Do not
preserve the wrapper workflow through compatibility adapters.

Expected result: substantial net deletion, including hundreds of source lines
and likely old tool/harness lines. No new implementation framework should
appear.

## Dangling Ladder Runner Follow-Up

Status: approved deletion authority under `HP-RETIRE-LADDER-RUNNERS-FN-01` and
`HP-RETIRE-LADDER-RUNNERS-TEST-01`.

After the route-driver materialization/report/save workflow retirement, the
ladder library is quarantined and no longer a live workflow. Two runner scripts
remain only as entrypoints into that retired workflow:

```text
tools/run_cartesian_driver_ladder.jl
tools/run_cartesian_line_ladder.jl
```

`HP-RETIRE-LADDER-RUNNERS-FN-01` approves deleting only those two files. It
does not approve replacements, canonical driver edits, source edits, test
edits, new tools, or changes to `tools/cartesian_driver_ladder_lib.jl`.

`tools/cartesian_driver_ladder_lib.jl` remains under the previous
`HP-RETIRE-DRV-MAT-TOOL-01` quarantine unless a later amendment explicitly
approves deleting it.

`HP-RETIRE-LADDER-RUNNERS-TEST-01` approves only:

- `git diff --check`;
- package load;
- focused `rg` over `src`, `bin`, `test`, and `tools` for
  `run_cartesian_driver_ladder`, `run_cartesian_line_ladder`, and
  `cartesian_driver_ladder_lib`;
- canonical small base artifact/readback smoke;
- no Cr2 run.

Failure rule: if any live source, canonical workflow, or approved tool still
depends on these runner scripts, make no commit and report the exact
dependency. Do not preserve the runners through an adapter.

After the runner deletion pass, pause this cleanup lane unless a separate
docs-only amendment names another remaining stale surface.
