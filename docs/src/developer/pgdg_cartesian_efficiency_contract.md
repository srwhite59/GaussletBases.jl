# PGDG Cartesian Efficiency Contract

## Purpose

This note records the intended efficiency contract for the PGDG Cartesian
production line in `GaussletBases`.

It is deliberately scoped.

This contract is about:

- basis construction
- exact raw matrix construction
- parent-to-final contraction structure

This contract is not about:

- final IDA closure
- residual-Gaussian nearest/MWG interaction rules
- numerical-reference validation backends
- the quarantined legacy 1D hybrid route

## Contract

On the PGDG production line:

1. distorted 1D primitives are replaced by plain Gaussian proxies before matrix
   assembly
2. primitive 1D operator tables are analytic
3. contracted 1D basis tables are formed only by contraction from those
   primitive analytic tables
4. Cartesian product-parent matrices are formed only as products or sums of
   products of those contracted 1D tables
5. product-parent bases are formed by contraction from a single Cartesian
   product basis
6. final bases are formed by direct sums of parent blocks followed by one final
   mixing contraction

On PGDG-capable pure molecular routes, per-center nuclear one-body terms are
therefore contracted directly into the exposed parent/final space; avoidable
dense parent 3D nuclear matrix materialization is outside the intended
contract.
On carried-plus-residual molecular supplement routes, the final one-body mix
should likewise exploit the carried identity block instead of treating the
`raw_to_final` map as an arbitrary dense congruence.

In the shorthand used during planning:

- `P1dF -> CP1dF -> CPB -> PPB -> FB`

## Scope

This contract applies directly to:

- the mapped PGDG proxy line in `src/ordinary_pgdg.jl`
- the contracted 1D PGDG intermediate/bundle layer in
  `src/ordinary_mapped_backends.jl`
- the compact nested Cartesian fixed-block production line built from those 1D
  tables

This contract does **not** yet claim that every current public ordinary/QW
route runs on the PGDG production lane. Many QW-facing entry points still
intentionally require `gausslet_backend = :numerical_reference`.

## Hybrid And Residual Routes

The broadened final-basis language is still useful for the current hybrid and
residual-Gaussian routes:

- the raw parent side can be read as a direct sum of product-parent sectors
- the final `raw_to_final` map is the final mixing contraction

But the downstream residual interaction closures remain outside this contract.
The bond-aligned diatomic molecular supplement route may still reuse that same
direct-contracted GG/fixed nuclear backbone, but its GA/AA supplement blocks
and residual-space closure remain outside this contract.

That separation is intentional. It lets the repo unify basis and raw-matrix
construction without pretending that all later approximate interaction models
already share one algebra.

## Validation And Reference Boundary

Numerical quadrature remains allowed in:

- validation/reference backends
- compatibility/reference paths that are explicitly marked as such

Numerical quadrature is **not** part of the PGDG production contract and should
not be reached by silent fallback on that lane.

The public mapped PGDG backends therefore use the local-linear analytic
Gaussian proxy layer:

- `:pgdg_experimental` builds `:mapped_pgdg_primitives`
- `:pgdg_localized_experimental` localizes that same primitive layer

Sampled log-fit or derivative-fit proxy helpers may remain as explicit
diagnostic/refinement tools, but they are outside the public production
backend path.

## Adoption Plan

The intended repo adoption order is:

1. document the contract clearly
2. add guardrails so the PGDG experimental/production lane refuses silent
   numerical primitive fallback or sampled proxy construction
3. keep `:numerical_reference` explicit and separate
4. later tighten metadata and naming once the production lane is fully wired
   through the public ordinary/nested surfaces

The first concrete implementation step is therefore narrow:

- enforce analytic primitive-layer expectations on the PGDG experimental lane
- leave broader QW and hybrid reformulation for later passes
