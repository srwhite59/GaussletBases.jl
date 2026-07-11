# Public ns Direct-Core Side Parity

Status: implemented under `HP-COMP-NSCORE-FN-01`; the bounded validation
contract `HP-COMP-NSCORE-TEST-01` is completed.

## Contract

`ns` is the public requested cube/source/nesting size. Route-local `q` is
derived by the selected nesting family:

```text
nesting = :pqs  -> q = ns
nesting = :wl   -> q = ns - 2
```

Direct nucleus-centered core identity blocks need an odd side so the nucleus is
centered on a grid point. Their side is therefore derived from public `ns`,
not route-local `q`:

```text
direct_core_side = isodd(ns) ? ns : ns + 1
```

This rule applies only to direct nucleus-centered core identity blocks.
Boundary shells, WL boundary strata, and other non-direct support regions keep
their route-local retained construction and are not oddized by this rule.

Examples:

```text
ns = 6, nesting = :pqs -> q = 6, direct_core_side = 7
ns = 6, nesting = :wl  -> q = 4, direct_core_side = 7
```

The WL boundary retained-count policy remains:

```text
ns = 4 -> 4^3 - 2^3 = 56
ns = 5 -> 5^3 - 3^3 = 98
ns = 6 -> 6^3 - 4^3 = 152
```

Thus same-`ns` PQS/WL comparisons share direct-core centering without
forcing equal boundary construction or final dimensions.

## Ownership

Source owners:

- `src/cartesian_base_hamiltonian.jl`:
  `_cartesian_base_direct_core_side(ns)` and the one-center parent minimum;
- `src/pqs_source_box_route_driver_helpers.jl`: direct-core route
  construction and truthful `:public_ns_direct_core_side` provenance.

The broader public `ns` and route-local `q` contract is owned by
[Nesting/supplement composition](nesting_supplement_composition_plan.md).

## Validation And Evidence

Implementation commit: `e41aba6eb` (`Fix public ns direct-core parity`).
Accepted bounded PQS/WL atom parity, diatomic smoke, dimension, provenance, and
readback evidence is recorded as Pass 160 in
`docs/src/developer/pqs_manager_running_log.md`.

No dedicated committed test file was added. Current public base regression
coverage lives in
`test/driver_public/cartesian_base_hamiltonian_runtests.jl`.

## Failure Behavior And Non-Goals

Invalid derived sizes fail before terminal construction. A future change must
not repair direct-core parity by changing boundary retained counts.

This contract does not authorize driver or public-input changes, route or
shellification redesign, terminal lowering, retained-unit or realizer changes,
artifact or manifest changes, old WL materialization, committed fixtures, or
Cr2-specific workflow.
