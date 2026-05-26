# Cartesian QW Receipt Wrapper Status

## Purpose

Record the current internal status of the Cartesian QW construction receipt
line.

The receipt layer is now coherent enough to be treated as the internal unified
wrapper for the route families it explicitly covers. It is still not a new
Hamiltonian implementation.

## Bottom line

`cartesian_qw_operator_construction_receipt(...)` is a safe outer wrapper
around existing QW builders.

For covered routes it does exactly this:

1. builds a `CartesianOperatorBuildSource3D`
2. calls the existing matching `ordinary_cartesian_qiu_white_operators(...)`
   builder
3. derives a `CartesianQWOperatorConstructionRecord3D`
4. returns a `CartesianQWOperatorConstructionReceipt3D`

The existing builder remains the authority for all numerical Hamiltonian
construction.

## Covered route families

The source of truth in code is
`cartesian_qw_operator_receipt_coverage()` in
`src/cartesian_qw_operator_carried_spaces.jl`.

As of this note, the receipt layer covers:

- bond-aligned direct-product routes without a Gaussian supplement
- bond-aligned nested fixed-block routes without a Gaussian supplement
- one-center atomic direct-product routes with `LegacyAtomicGaussianSupplement`
- one-center atomic nested fixed-block routes with
  `LegacyAtomicGaussianSupplement`
- bond-aligned diatomic direct-product routes with molecular legacy Gaussian
  supplements, including homonuclear and heteronuclear supplement objects
- bond-aligned diatomic nested fixed-block routes with molecular legacy
  Gaussian supplements, including homonuclear and heteronuclear supplement
  objects

These are receipt-covered because the repo already has matching authoritative
builders for them. The receipt layer only delegates and audits.

## Intentionally uncovered route families

The receipt layer deliberately does not accept:

- alias/frontend function objects such as `ordinary_cartesian_product_operators`
  or `nested_cartesian_operators`
- already-built `OrdinaryCartesianOperators3D` payloads
- non-diatomic molecular supplement inputs, such as chain or square routes with
  molecular legacy Gaussian supplements
- mapped one-center atomic routes without a Gaussian supplement

Alias/front-door APIs remain direct-builder-only for now. If a caller wants a
receipt, it should call the receipt helper with the underlying canonical basis,
fixed block, and optional supplement inputs.

Already-built operator payloads should use the post-build audit helpers:

- `cartesian_qw_operator_construction_record(...)`
- `cartesian_qw_operator_carried_space_sidecar(...)`

Uncovered scientific route families should stay uncovered until a matching
authoritative builder exists and the route is deliberately added to the
receipt coverage table.

## Why this is safe

The receipt layer does not:

- implement Hamiltonian kernels
- choose a new backend
- alter geometry policy
- build metric packets
- materialize dense parent matrices
- change residual or supplement semantics
- change numerical-reference fallback behavior

Tests assert that receipt-built operators match direct-builder operators for
representative covered routes. Receipt diagnostics also record that the wrapper
delegated to an existing builder and did not use a new Hamiltonian kernel.

## PGDG and numerical-reference boundary

PGDG analytic routing remains the production contract where a route supports
it. Numerical quadrature remains an explicit `:numerical_reference` or
diagnostic route, not something the receipt wrapper silently selects.

The receipt layer preserves whatever backend semantics the existing builder
already enforces. It does not widen or reinterpret backend support.

## Relation to Cartesian/Hamiltonian unification

The receipt line is the current bridge between the normalized Cartesian build
metadata and the existing Hamiltonian builders.

It proves that a single internal shape can:

- record a request before construction
- run the current route-specific builder
- audit the produced operator payload against the request

This is the appropriate seam for later consumer migration or a future facade.
That later work should still call proven builders until a generic Hamiltonian
builder exists and is validated.

## Next engineering step

The next step is not to add more science routes opportunistically.

The useful next step is selective consumer migration: use receipt wrappers in
internal callers that benefit from source/record/provenance diagnostics, while
continuing to treat existing builders as authoritative.

New route coverage should be added only when:

- the underlying direct builder already exists
- the route has a clear coverage-table entry
- receipt output is tested for exact numerical equality against the direct
  builder
- unsupported aliases or ambiguous inputs still fail clearly
