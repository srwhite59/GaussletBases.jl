# QW And High-Order Experimental Cluster Retirement

## Status And Authority

This page owns the retirement contract for:

- `HP-RETIRE-QW-DONOR-FN-01`
- `HP-RETIRE-QW-DONOR-TEST-01`

Status: source retirement completed and validated in Pass 409. The function ID
is retired, and the validation ID is completed; neither grants further work.

This authority removes an unsuccessful experimental implementation. It does
not abandon homonuclear linear atomic chains as a future scientific target and
does not change the current `nesting=:wl` producer.

## Retirement Decision

The retirement deleted these four files together:

- `src/cartesian_qw_operator_carried_spaces.jl`
- `src/ordinary_qw_experimental_paths.jl`
- `src/cartesian_high_order_doside_experimental.jl`
- `src/cartesian_high_order_doside_ida_experimental.jl`

They contained `6,008` source lines at the approved baseline. Their exported
surfaces are explicitly experimental, have no committed source or test
consumers outside this cluster, and are not compatibility obligations.

The same retirement removed from `src/GaussletBases.jl`:

- the four includes for those files;
- the two exported experimental path types;
- the six exported chain/square nested operator, diagnostic, and report names;
- the corresponding six empty generic declarations;
- the `CartesianQWOperatorCarriedSpaces` submodule import/export surface and
  its 24 exported names.

Do not retain any of these through stubs, deprecation wrappers, aliases, moved
implementations, or compatibility modules. Ignored probes may remain stale
historical evidence; source adapters must not be added to keep them runnable.

## Caller And Ownership Audit

At baseline `78cb6f806`:

- committed `src/` and `test/` contain no caller of the eight top-level
  experimental exports outside the four-file cluster and `GaussletBases.jl`;
- `CartesianQWOperatorCarriedSpaces` is consumed only by
  `ordinary_qw_experimental_paths.jl`;
- current Cartesian base and `nesting=:wl` construction do not include or call
  this cluster;
- the active chain/square basis types and geometry diagnostics are owned by
  `ordinary_qw_types_and_bases.jl`, not by the retiring implementation.

No live committed caller was found. Restoring any retired surface requires new
authority; compatibility glue remains forbidden.

## Preserved Source

The retirement must not delete or reinterpret:

- chain/square basis types and constructors in
  `src/ordinary_qw_types_and_bases.jl`;
- shared generics and surviving diatomic methods in
  `src/ordinary_qw_nested_frontends.jl`;
- active geometry and retained-unit helpers in
  `src/cartesian_nested_faces.jl`,
  `src/cartesian_nested_owned_units.jl`, and
  `src/cartesian_nested_diatomic.jl`;
- `src/ordinary_qw_residuals.jl`, `src/ordinary_qw_raw_blocks.jl`, and
  `src/ordinary_qw_operator_assembly.jl`;
- `src/cartesian_qw_hybrid_representation.jl`;
- `src/cartesian_carried_spaces.jl`, whose API and callers require a separate
  audit.

The current WL, PQS, shellification, residual-GTO, MWG, Hamiltonian, artifact,
and driver behavior must remain unchanged.

## Historical Scientific Evidence

The retired implementation established useful but insufficient evidence:

- a doside ladder and full tensor-shell construction;
- the distinction between full shell basis (FSB) and full basis union (FBU);
- successful undistorted He+ and He control calculations;
- an unresolved benefit question on distorted mapped parents;
- bounded Cr parent-orbital capture diagnostics;
- endcap and panel concepts later represented in active adjacent geometry and
  retained-unit code.

Those results remain design history, not source authority. Atomic chains remain
a long-term target, but any future chain implementation starts from a new
design and the preserved basis/geometry contracts rather than restoring this
cluster.

## Accepted Validation

The source retirement pass verified:

1. no audited include, export, generic, module, type, or function remained in
   committed source;
2. package load succeeded;
3. existing focused core tests continued to protect chain/square basis construction
   and geometry diagnostics;
4. the public Cartesian base test protected the current `nesting=:wl` path;
5. terminal due diligence for that WL endpoint was inspected and reported;
6. `git diff --check` passed;
7. source deletion was `6,026` net lines, with no new tests.

Existing tests may be edited only to remove a stale reference discovered during
deletion. Do not add replacement tests for retired APIs.

## Non-Goals

This retirement does not authorize:

- redesigning atomic chains or square lattices;
- changing active high-order adjacent geometry helpers;
- deciding the fate of `cartesian_carried_spaces.jl`;
- modifying WL/PQS numerical behavior, source modes, shellification, residuals,
  MWG, IDA, artifacts, drivers, or solvers;
- touching `src/hamiltonian_corrections.jl` or successor handoffs.
