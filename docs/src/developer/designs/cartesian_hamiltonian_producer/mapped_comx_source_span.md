# Mapped-COMX Source Span

Status: approved narrow mainline source-span authority under
`HP-MCOMX-FILE-01`, `HP-MCOMX-OBJ-01`, `HP-MCOMX-FN-01`,
`HP-MCOMX-WIRE-01`, and `HP-MCOMX-TEST-01`.

## Purpose

The high-order lane is an experimental proving ground. It produced the
mapped-COMX source-span idea and scratch evidence that the idea is worth
installing once in mainline so real atom gates can test it through the current
carried-space / PQS / PGDG machinery.

Mainline owns the facility and contract. High-order remains a consumer and
benchmark lane, like CR2: it should call the installed mainline option rather
than maintain a duplicate scratch implementation.

## Numerical Rule

The approved first source-span option is:

```text
protected physical P2
+ mapped Chebyshev enrichment T_k(s_lambda(u))
+ lambda = 0.5
+ no sqrtJ
+ physical-u COMX localization
```

For a one-dimensional local physical coordinate `u`:

```text
s_lambda(u) = (1 + lambda) * u / (1 + lambda * u^2)
```

Default source-span specification:

```text
protected_degree = 2
lambda = 0.5
mapped_family = :chebyshev_s
include_sqrt_jacobian = false
localization_coordinate = :physical_u
```

Construction semantics:

1. Build the protected physical polynomial block `1, u, u^2`.
2. Add mapped Chebyshev columns `T_k(s_lambda(u))` until the requested source
   mode count is reached.
3. Project mapped columns against the protected physical block in the local
   parent metric.
4. Orthonormalize the combined source span.
5. Localize by physical-coordinate COMX, not mapped-`s` COMX.

The option is general in `n_s` and `protected_degree`. The tested
`n_s = 5, 6, 7` values are evidence points, not hard-coded cases.

## Evidence Summary

High-order-manager scratch scans on 2026-06-25 found that pure mapped columns
improve angular center spacing on cube faces, while protected `P2 + T` gives a
better first occupied-mode source span than stronger `P3`/`P4` protection at
small retained counts.

The proposed default `lambda = 0.5` was a balanced setting across the reported
low angular proxy rows:

```text
n_s = 5: L1 = 1.431e-3, L2 = 7.529e-3, Eg = 2.238e-3, T2g = 9.547e-3
n_s = 6: L1 = 3.391e-4, L2 = 1.446e-3, Eg = 2.238e-3, T2g = 3.863e-4
n_s = 7: L1 = 9.702e-5, L2 = 2.818e-4, Eg = 4.105e-4, T2g = 1.417e-4
```

No Hamiltonian, IDA, HF, DMRG, ECP, EGOI, MWG, or Cr2 calculation was used as
approval evidence. Those are later consumer gates after the source-span option
is installed in mainline.

## Approved IDs

- `HP-MCOMX-FILE-01` - source files for the mainline mapped-COMX source-span
  option.
- `HP-MCOMX-OBJ-01` - compact `MappedCOMXSourceSpec` or equivalent source
  specification object.
- `HP-MCOMX-FN-01` - source-span and axis-transform construction.
- `HP-MCOMX-WIRE-01` - narrow wiring into existing raw product source / PQS
  axis-transform construction.
- `HP-MCOMX-TEST-01` - validation gates.

## Approved Source Surface

Primary owner:

```text
CartesianRawProductSources
```

Approved files:

```text
src/cartesian_raw_product_sources/CartesianRawProductSources.jl
src/cartesian_raw_product_sources/mapped_comx_source_span.jl
src/cartesian_raw_product_sources/axis_transform_facts.jl
src/cartesian_raw_product_sources/records.jl
src/cartesian_pair_block_materialization/pqs_source_axis_transforms.jl
```

`records.jl` is approved only for narrow accessors, metadata fields, or record
validation needed to carry source-span provenance. `pqs_source_axis_transforms.jl`
is approved only for narrow compatibility wiring from existing PQS raw-source
axis-transform construction to the new mainline source-span option.

No root include change is approved except the module-local include in
`CartesianRawProductSources.jl`.

## Approved Behavior

The implementation may:

- add a compact `MappedCOMXSourceSpec` or equivalent typed source-span object;
- construct the one-dimensional source columns for protected physical
  polynomials plus mapped Chebyshev enrichment;
- project mapped columns against the protected physical block in the local
  parent metric;
- build materialized `AxisSourceTransformFact`s whose coefficient matrices are
  compatible with existing `RawProductBoxPlan` and PQS boundary product-mode
  retained rules;
- localize by physical `u` COMX;
- record compact metadata:
  - `source_span_family = :mapped_comx`;
  - `protected_degree`;
  - `lambda`;
  - `mapped_family`;
  - `mapped_orders`;
  - `include_sqrt_jacobian = false`;
  - `localization_coordinate = :physical_u`;
  - requested/resolved source mode count;
  - rank and overlap/orthogonality diagnostics;
- leave ordinary polynomial source spans available and unchanged;
- expose the option only through internal construction controls needed by the
  approved validation gates.

The Hamiltonian, operator, and artifact layers should consume the same
carried-space / raw product source facts as before. They should not branch on
whether a source span was ordinary polynomial or mapped-COMX, except for
descriptive provenance already carried by source facts.

## Forbidden

This amendment does not approve:

- changing default source spans;
- public API or export changes;
- canonical driver input changes;
- artifact schema, manifest, or reader changes;
- Hamiltonian, one-body, IDA, MWG, Residual Gaussian, raw Gaussian block, or
  solver changes;
- ECP, EGOI, RHF, ED, DMRG, or Cr2 workflow;
- explicit `Y_lm` / angular injection;
- `sqrtJ` weighting;
- mapped-`s` COMX as the production localization gauge;
- importing high-order branch scaffolding, scripts, route wrappers, status
  objects, diagnostics, or reports;
- duplicate high-order-maintained implementation of the same mainline option;
- committed Cr fixtures or broad high-order benchmark fixtures.

## Validation

`HP-MCOMX-TEST-01` approves only:

- `git diff --check`;
- package load;
- local source-span validation for `n_s = 5`, `6`, and `7`:
  - full retained rank;
  - protected `P2` span preserved;
  - source-axis overlap approximately identity after construction;
  - physical-`u` COMX off-diagonal residual reported;
  - metadata records the approved source-span rule;
- a bounded cubic H and He+ one-electron gate comparing ordinary polynomial
  and mapped-COMX source spans with fixed support and retained count;
- a bounded He `1s^2` fixed-orbital IDA gate if the already-supported
  analytic path can run without new solver or artifact workflow;
- high-order-manager consumer benchmarks on the installed mainline option for
  Cr occupied capture, reported back as evidence rather than committed mainline
  fixtures;
- no Cr2 run.

No committed test file is approved by default. A later implementation blurb may
name a small standalone script or ignored probe when needed for the approved
gates. Committed fixtures, public driver tests, solver tests, and Cr/Cr2
benchmark fixtures require a separate amendment.

## Failure Rule

If the source option cannot be installed through `CartesianRawProductSources`
and existing PQS raw-source axis-transform wiring without changing Hamiltonian
assembly, artifact schemas, public driver inputs, or high-order-specific
workflow, make no source commit and report the exact missing mainline seam.

If the real atom gates fail after a faithful implementation of the approved
source-span rule, do not tune defaults or add injection in the same pass.
Report the measured failure and request a separate design amendment.
