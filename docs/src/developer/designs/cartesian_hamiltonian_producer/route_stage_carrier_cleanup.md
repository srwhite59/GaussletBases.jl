# Route/Stage Carrier Cleanup

Status: approved for implementation under `HP-ROUTE-STAGE-CARRIER-FN-01` and
`HP-ROUTE-STAGE-CARRIER-TEST-01`.

## Purpose

After `HP-ROUTE-STAGE-TYPE-FN-01`, warm Be2 q5 p10 supplemented construction
remains fast, but post-cleanup attribution on `118a639b` still shows cold
latency in route/stage type specialization. The prior targeted inventory
owners moved out of the dominant set. The remaining route/stage owners are now
broader stage signatures and plan carriers around `cartesian_shells`,
`cartesian_units`, `cartesian_transforms`, terminal topology support-region
planning, and terminal retained-rule planning.

This lane allows a narrow follow-up cleanup of stage-carried shellification,
route-skeleton, support-plan, retained-rule-plan, and terminal-plan carriers
where they inflate route/stage function signatures.

## Approved Boundary

Approved source files:

```text
src/pqs_source_box_route_driver_helpers.jl
src/pqs_source_box_diatomic_complete_core_shell.jl
```

Optional source file, only if directly required to slim the approved path:

```text
src/cartesian_final_basis_realization/pqs_terminal_basis_realization.jl
```

Approved targets:

- `cartesian_shells` stage carrier and return signature;
- `cartesian_units` stage carrier and return signature;
- `cartesian_transforms` stage carrier and return signature;
- terminal topology support-region planning in
  `src/pqs_source_box_diatomic_complete_core_shell.jl`;
- terminal retained-rule planning in
  `src/pqs_source_box_diatomic_complete_core_shell.jl`;
- terminal realization plan carriers in
  `src/cartesian_final_basis_realization/pqs_terminal_basis_realization.jl` only
  where directly required to avoid reintroducing a large stage-carried plan
  shape through the approved route/stage path.

Approved changes:

- stop carrying giant shellification, route-skeleton, support-plan,
  retained-rule-plan, and terminal-plan `NamedTuple` / tuple shapes across the
  approved stage function signatures;
- replace necessary carriers with compact typed/vector-backed records, stable
  dictionaries, accessors, or smaller summaries;
- recompute small derived summaries from canonical objects inside the approved
  path where that is simpler than carrying wide stage payloads;
- delete stale compatibility carriers with no active approved caller;
- preserve deterministic terminal support, shellification, and lowering order.

Route skeleton construction semantics are not changed by this lane. If changing
`src/pqs_source_box_route_driver_skeletons.jl` is required, stop and request a
separate amendment.

## Required Preservation

The source pass must preserve:

- H2 base artifact/readback behavior;
- H2 supplemented artifact/readback behavior;
- H2 R3 endpoint if terminal realization is touched;
- deterministic terminal support/shellification/lowering order;
- existing public driver contract;
- existing artifact schema and manifest behavior;
- existing route semantics and numerical matrices.

## Forbidden

This lane does not approve:

- source files outside the approved boundary;
- edits to `src/pqs_source_box_route_driver_skeletons.jl`;
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

No compatibility adapter may preserve the old runtime-sized carrier merely
under a new name.

## Validation

Required validation:

- `git diff --check`;
- package load;
- H2 base artifact write/readback;
- H2 supplemented artifact write/readback;
- H2 R3 endpoint if terminal realization is touched;
- focused terminal support/shellification/lowering order parity;
- focused scan for newly introduced runtime-sized `NamedTuple{...}`,
  `Tuple(...)`, `Tuple{Vararg{...}}`, and runtime-keyed inventories in the
  approved files;
- no Cr2 run.

Optional after correctness passes:

- Be2 q5 post-cleanup compile/timing comparison using the same explicit
  p10-style supplemented driver fixture that produced the attribution.

Existing committed tests may be adjusted only if they directly assert the old
stale carrier shape. No new committed test file is approved.

Line budget: at most `250` added `src` lines, with net simplification expected.

Failure rule: if cleanup requires source files outside the approved boundary,
broad route-stage redesign, public API changes, artifact changes, numerical
changes, or precompile/sysimage machinery, make no source commit and report the
exact blocker.
