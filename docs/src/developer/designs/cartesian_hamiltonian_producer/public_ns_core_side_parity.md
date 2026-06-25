# Public ns Direct-Core Side Parity

Status: approved narrow source authority under `HP-COMP-NSCORE-FN-01` and
`HP-COMP-NSCORE-TEST-01`.

## Problem

`ns` is the public requested cube/source/nesting size. Route-local `q` is
derived from `ns`:

```text
nesting = :pqs  -> q = ns
nesting = :wl   -> q = ns - 2
```

The shared route setup still derives the nucleus-centered direct/core side from
route-local `q`. That makes even-`ns` same-size PQS/WL comparisons use
different direct core boxes. For example, `ns = 6` gives PQS a `7^3` direct
core side but WL a `5^3` direct core side, while the WL boundary shell retained
count remains `6^3 - 4^3 = 152`.

This is an odd/even convention mismatch, not a physics difference.

## Decision

The odd-side parity rule applies only to direct nucleus-centered core identity
blocks. It is needed there so the nucleus remains centered on a grid point.

The direct core side must be derived from public `ns`, not route-local `q`:

```text
direct_core_side = isodd(ns) ? ns : ns + 1
```

Boundary/shell retained sizes still use the appropriate route-local
construction. They must not inherit direct-core oddization.

Examples:

```text
ns = 6, nesting = :pqs -> q = 6, direct_core_side = 7
ns = 6, nesting = :wl  -> q = 4, direct_core_side = 7
```

The WL boundary retained count policy remains:

```text
ns = 4 -> 4^3 - 2^3 = 56
ns = 5 -> 5^3 - 3^3 = 98
ns = 6 -> 6^3 - 4^3 = 152
```

## Approved IDs

- `HP-COMP-NSCORE-FN-01` - direct core side from public `ns`.
- `HP-COMP-NSCORE-TEST-01` - parity validation gates.

## Approved Source Surface

Approved later source files:

```text
src/pqs_source_box_route_driver_helpers.jl
src/cartesian_base_hamiltonian.jl
```

`src/cartesian_base_hamiltonian.jl` is approved only if needed to keep
one-center parent minimum sizing consistent with the same direct-core rule.

## Approved Behavior

- keep public `ns` as the requested cube/source/nesting size;
- keep route-local `q` derivation unchanged:
  - PQS: `q = ns`;
  - WL: `q = ns - 2`;
- derive only direct nucleus-centered core side from public `ns`:
  `direct_core_side = isodd(ns) ? ns : ns + 1`;
- do not apply this oddization rule to boundary shells, WL boundary-stratum
  retained products, or any non-direct support region;
- preserve odd-`ns` dimensions where this cleanup should be a no-op;
- preserve WL boundary retained count policy, including `ns = 6` boundary shell
  count `152`;
- update internal provenance or summary labels such as `:odd_q_core_side` to an
  `ns`-truthful rule if they are still written.

## Forbidden

This amendment does not approve:

- driver changes;
- public input changes;
- route skeleton redesign;
- shellification geometry rewrites beyond direct-core side authority;
- terminal lowering, retained-unit, or terminal-realizer changes;
- artifact schema changes;
- manifest expansion;
- PQS/WL special cases in the driver;
- old WL materialization revival;
- committed tests or fixtures;
- Cr2-specific workflow.

## Validation

`HP-COMP-NSCORE-TEST-01` approves only:

- `git diff --check`;
- package load;
- small one-center atom base artifact/readback for `ns = 5`, `6`, and `7` with
  `nesting = :pqs` and `nesting = :wl`;
- same-`ns` PQS/WL atom dimensions match for `ns = 5`, `6`, and `7` under a
  bounded common fixture;
- `ns = 6` no longer has the reported `66`-row skew;
- small H2 or Be2 smoke to confirm the diatomic path still constructs;
- existing H2 Residual Gaussian endpoint smoke if touched code crosses the
  supplemented path;
- no Cr2 run.

## Failure Rule

If fixing the parity requires changing route skeleton semantics, terminal
lowering, retained-unit records, WL boundary coefficient construction, artifact
schema, or driver inputs, make no source commit and report the exact blocker.

Separate follow-up: WL manifest source-shell/source-mode provenance is still
weaker than PQS. That should be handled by a later provenance lane, not bundled
with this parity fix.
