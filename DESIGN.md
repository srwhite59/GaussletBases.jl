# GaussletBases.jl — Design Snapshot

_Status: working draft based on current interface discussion. This file is meant to be the single evolving design note for the first public version._

## Goals

GaussletBases.jl should be a **small multilevel library** with a clear interface.

It should support three levels of use:

1. **function objects**
2. **basis objects**
3. **ready-to-use operator matrices**

It should also expose the exact primitive stencil behind each function and each basis, so users can work at the analytic primitive layer when needed.

The design is intentionally smaller than a full workflow package. The first release should focus on the parts that are already conceptually clear and avoid unnecessary infrastructure.

---

## Main design decisions locked for now

- Julia is the launch reference implementation.
- The library uses **callable function objects** as the main low-level interface.
- The public interface uses **one public stencil level**.
- `stencil(f)` means the exact defining linear combination of the **lowest public primitive functions**.
- `basis[i]` returns a lightweight concrete callable function object.
- Radial work keeps **basis and quadrature as separate objects**.
- The primitive analytic layer is exposed basis-wide as:
  - `primitives(basis)`
  - `stencil_matrix(basis)`
- Primitive-to-basis contraction is exposed explicitly with:
  - `contract_primitive_vector`
  - `contract_primitive_diagonal`
  - `contract_primitive_matrix`
- Distortion/mapping is represented explicitly as a wrapper:
  - `Distorted(f, mapping)`
- For radial primitives, the preferred public names are:
  - `HalfLineGaussian`
  - `XGaussian`
- The radial reference-coordinate spacing is called `reference_spacing`, not `u_spacing`.

---

## Public object model

```julia
abstract type AbstractFunction1D end
abstract type AbstractPrimitiveFunction1D <: AbstractFunction1D end
abstract type AbstractBasisFunction1D    <: AbstractFunction1D end

abstract type AbstractCoordinateMapping end
abstract type AbstractBasisSpec end
```

Common protocol:

```julia
(f::AbstractFunction1D)(x::Real) = value(f, x)

value(f::AbstractFunction1D, x::Real)
direct_value(f::AbstractFunction1D, x::Real)
derivative(f::AbstractFunction1D, x::Real; order::Int = 1)

center(f::AbstractFunction1D)
reference_center(f::AbstractFunction1D)
moment_center(f::AbstractFunction1D, grid)
integral_weight(f::AbstractFunction1D)

stencil(f::AbstractFunction1D)
```

### Meaning of the common interface

- `f(x)` is the primary user-facing evaluation style.
- `value(f, x)` is the named form.
- `direct_value(f, x)` means evaluation directly from the exact defining stencil, bypassing interpolation or caches.
- `center(f)` means the physical-space X-center.
- `reference_center(f)` means the center before coordinate mapping.
- `moment_center(f, grid)` is a real diagnostic quantity, especially near boundaries or the origin.
- `integral_weight(f)` means the integral of the basis function, not a stencil coefficient or quadrature weight.

---

## Primitive function objects

```julia
Gaussian(; center::Real = 0.0, width::Real)

HalfLineGaussian(; center::Real = 0.0, width::Real)

XGaussian(; alpha::Real)

Distorted(f::AbstractPrimitiveFunction1D,
          mapping::AbstractCoordinateMapping)
```

### Notes

- Primitive objects describe **shapes**, not weighted terms.
- Scalar coefficients live in stencil terms, not in primitive objects.
- `HalfLineGaussian` is preferred over names like `PositiveGaussian`.
- `XGaussian` is the near-origin primitive used in the radial construction.
- `Distorted(f, mapping)` is the current design choice for mapped primitives.

---

## Gausslet families and low-level gausslet objects

```julia
GaussletFamily(name::Symbol)

Gausslet(family::Union{GaussletFamily,Symbol};
         center::Real = 0.0,
         spacing::Real = 1.0)

Gausslet(coefficients::AbstractVector{<:Real};
         center::Real = 0.0,
         spacing::Real = 1.0)
```

### Notes

- Built-in family names include `:G4`, `:G6`, `:G8`, `:G10`.
- Custom coefficient-vector construction remains available for expert workflows.
- A `Gausslet` is a callable function object whose exact definition is its Gaussian stencil.
- Fast evaluation may use interpolation internally, but interpolation is not the mathematical definition.

---

## Coordinate mappings

```julia
IdentityMapping()

AsinhMapping(; c::Real,
             s::Real,
             tail_spacing::Real = 10.0)

AsinhMapping(; a::Real,
             s::Real,
             tail_spacing::Real = 10.0)
```

Mapping interface:

```julia
uofx(mapping, x)
xofu(mapping, u)
dudx(mapping, x)
du2dx2(mapping, x)

(mapping::AbstractCoordinateMapping)(x) = uofx(mapping, x)
```

### Notes

- Public call syntax means the **forward map**, not the density.
- Keep `uofx`, `xofu`, `dudx`, and `du2dx2` as explicit math-facing names.
- Mappings are first-class public objects.

---

## Basis recipes

```julia
UniformBasisSpec(family::Union{GaussletFamily,Symbol};
                 xmin::Real,
                 xmax::Real,
                 spacing::Real = 1.0)

HalfLineBasisSpec(family::Union{GaussletFamily,Symbol};
                  xmax::Real,
                  reference_spacing::Real = 1.0,
                  tails::Int = 6,
                  mapping::AbstractCoordinateMapping = IdentityMapping())

RadialBasisSpec(family::Union{GaussletFamily,Symbol};
                rmax::Real,
                mapping::AbstractCoordinateMapping,
                reference_spacing::Real = 1.0,
                tails::Int = 6,
                odd_even_kmax::Int = 6,
                xgaussians::AbstractVector{XGaussian} = XGaussian[])

RadialBasisSpec(family::Union{GaussletFamily,Symbol};
                count::Int,
                mapping::AbstractCoordinateMapping,
                reference_spacing::Real = 1.0,
                tails::Int = 6,
                odd_even_kmax::Int = 6,
                xgaussians::AbstractVector{XGaussian} = XGaussian[])
```

Construction entry point:

```julia
build_basis(spec::AbstractBasisSpec)
build_basis(spec::HalfLineBasisSpec; grid_h=nothing, refine_grid_h=true)
build_basis(spec::RadialBasisSpec; grid_h=nothing, refine_grid_h=true)
```

### Notes

- `BasisSpec` stays explicit; do not shorten it to bare `Spec` in user-facing docs.
- Most users should prefer `rmax`; `count` is mainly for reproducibility, testing, and benchmarking.
- `tails`, `odd_even_kmax`, and `xgaussians` are public because they are method-level choices, not tiny implementation details.
- `grid_h` is a build-time control keyword on `build_basis`, not part of the mathematical `BasisSpec`.

---

## Basis objects

```julia
UniformBasis
HalfLineBasis
RadialBasis
```

Shared accessors:

```julia
length(basis)
basis[i]

basis_spec(basis)
family(basis)
mapping(basis)

centers(basis)
reference_centers(basis)
integral_weights(basis)

basis_diagnostics(basis)
basis_diagnostics(basis, grid)
```

### `basis[i]` return types

```julia
UniformBasis  -> Gausslet
HalfLineBasis -> BoundaryGausslet
RadialBasis   -> RadialGausslet
```

The object returned by `basis[i]` is a lightweight concrete function object that behaves like a standalone callable function but may share cached data with the parent basis.

---

## One public stencil layer

The library exposes one public stencil level.

### Public rule

`stencil(f)` always means the exact expansion of `f` in the **lowest public primitive functions**.

Examples:

- `stencil(Gausslet(...))` expands over `Gaussian`s.
- `stencil(BoundaryGausslet(...))` expands over `HalfLineGaussian`s, with `Distorted(...)` wrappers if mapped.
- `stencil(RadialGausslet(...))` expands over `HalfLineGaussian`, optional `XGaussian`, and `Distorted(...)` wrappers if mapped.

### Stencil objects

```julia
StencilTerm
FunctionStencil
```

Public behavior:

```julia
st = stencil(f)

length(st)
st[j]

terms(st)
coefficients(st)
primitives(st)

coefficient(term)
primitive(term)

st(x)
```

### Notes

- `FunctionStencil` is callable.
- `stencil(f)` may return a lightweight view into parent basis data rather than allocating fresh storage.
- The stencil is not a full internal construction trace.

---

## Basis-wide primitive layer

```julia
primitives(basis)
stencil_matrix(basis)
```

Contract:

```julia
basis[i](x) = sum(stencil_matrix(basis)[μ, i] * primitives(basis)[μ](x)
                  for μ in eachindex(primitives(basis)))
```

### Convention

- Rows correspond to `primitives(basis)`.
- Columns correspond to basis functions in the same order as `basis[i]`.

This is the basis-wide analogue of `stencil(basis[i])`.

---

## Primitive contraction helpers

```julia
contract_primitive_vector(basis, vμ)
contract_primitive_diagonal(basis, dμ)
contract_primitive_matrix(basis, Aμν)
```

Conventions:

```julia
C = stencil_matrix(basis)

contract_primitive_vector(basis, vμ)   == C' * vμ
contract_primitive_diagonal(basis, dμ) == C' * Diagonal(dμ) * C
contract_primitive_matrix(basis, Aμν)  == C' * Aμν * C
```

### Notes

- Primitive-space ordering is always defined by `primitives(basis)`.
- This is the main expert entry point for analytic or mixed analytic/quadrature workflows.
- For v0, these helpers replace a larger public backend abstraction.

---

## Quadrature

```julia
RadialQuadratureGrid(points, weights; mapping = nothing)

radial_quadrature(basis::RadialBasis;
                  accuracy = :high,
                  refine = nothing,
                  quadrature_rmax = nothing)

quadrature_points(grid)
quadrature_weights(grid)
```

### Notes

- Radial quadrature is a separate object from the radial basis.
- The quadrature grid follows the same mapping as the basis but is usually finer.
- This is not a DVR basis.
- The default call should choose a conservative quadrature cutoff and a default starting resolution automatically.
- `accuracy` selects one of the built-in quadrature-accuracy profiles `:medium`, `:high`, or `:veryhigh`, with `:high` as the default.
- `refine` is an optional starting refinement hint relative to the reference-coordinate basis spacing.
- `quadrature_rmax` is an optional explicit physical-space cutoff override for expert use.

---

## High-level radial operators

Approximation choices:

```julia
abstract type AbstractDiagonalApproximation end

IntegralDiagonal()
```

Matrix builders:

```julia
overlap_matrix(basis, grid)
kinetic_matrix(basis, grid)

nuclear_matrix(basis::RadialBasis, grid; Z::Real)
centrifugal_matrix(basis::RadialBasis, grid; l::Int)

multipole_matrix(basis::RadialBasis,
                 grid;
                 L::Int,
                 approximation::AbstractDiagonalApproximation = IntegralDiagonal())
```

High-level bundle:

```julia
RadialAtomicOperators

atomic_operators(basis::RadialBasis,
                 grid::RadialQuadratureGrid;
                 Z::Real,
                 lmax::Int = 0,
                 approximation::AbstractDiagonalApproximation = IntegralDiagonal())
```

Preferred access pattern:

```julia
ops.overlap
ops.kinetic
ops.nuclear

centrifugal(ops, l)
multipole(ops, L)
```

### Notes

- For v0, `multipole_matrix` means the supported two-index IDA-style radial
  multipole matrix. An exact non-diagonal electron-electron object is a
  separate future API question and should not be folded into the same function
  name now.
- Use fields for the fixed matrices.
- Use accessor functions for the angular-momentum-indexed pieces.
- Avoid raw public `ops.Vcentr[ℓ]` / `ops.Veel[L]` style indexing in the high-level API.

---

## Source-ready docstring drafts

### `FunctionStencil`

```julia
"""
    FunctionStencil

Exact public stencil of a 1D function in terms of the lowest-level primitive
function objects exposed by the library.

A `FunctionStencil` represents a linear combination

    f(x) = Σ_j c_j ϕ_j(x)

where `c_j` are scalar coefficients and `ϕ_j` are primitive function objects.
The stencil is itself callable and evaluates to the same function as the object
from which it was obtained.

This is the one public stencil level of the library. It is intended to answer
the question "what simpler functions is this function exactly made from?" It is
not a full construction trace of the basis-building algorithm.

Examples of primitive types appearing in a stencil include:
- `Gaussian` for ordinary full-line gausslets
- `HalfLineGaussian` for half-line and radial constructions
- `XGaussian` for the extra near-origin radial primitives
- `Distorted(...)` wrappers when a coordinate mapping is part of the primitive

Use:
- `coefficients(st)` to get the scalar coefficients
- `primitives(st)` to get the primitive function objects
- `terms(st)` to iterate over coefficient/primitive pairs

See also: `stencil`, `primitives`, `stencil_matrix`, `direct_value`
"""
```

### `stencil`

```julia
"""
    stencil(f::AbstractFunction1D) -> FunctionStencil

Return the exact public stencil of `f`.

For a `Gausslet`, the stencil is its contraction over ordinary `Gaussian`
primitives. For half-line or radial basis functions, the stencil is over the
corresponding half-line primitive layer, including any `XGaussian` terms and
any explicit `Distorted(...)` wrappers required by the coordinate mapping.

`stencil(f)` is intended for inspection, debugging, documentation, and
primitive-level integration workflows. It should not be interpreted as a full
history of how `f` was constructed internally.

The returned stencil may be a lightweight view into cached basis data rather
than a newly allocated object.

See also: `FunctionStencil`, `direct_value`, `primitives`, `stencil_matrix`
"""
```

### `RadialBasisSpec`

```julia
"""
    RadialBasisSpec(family; rmax, mapping, reference_spacing=1.0, tails=6,
                    odd_even_kmax=6, xgaussians=XGaussian[])

    RadialBasisSpec(family; count, mapping, reference_spacing=1.0, tails=6,
                    odd_even_kmax=6, xgaussians=XGaussian[])

Recipe for constructing a radial gausslet basis on `[0, ∞)` for the reduced
radial function `u(r) = r R(r)`.

The resulting basis is intended to satisfy the radial boundary condition at the
origin directly in the basis itself, so that each basis function vanishes at
`r = 0`.

Arguments
=========

- `family`:
  canonical gausslet family, such as `:G6`, `:G8`, or `:G10`, or an explicit
  `GaussletFamily` object.

- `mapping`:
  coordinate mapping controlling variable resolution in physical radius `r`.
  The basis is uniform in the reference coordinate and mapped into `r`.

Keyword arguments
=================

- `rmax`:
  preferred user-facing size control. Build enough functions to represent the
  interval from the origin out to approximately `rmax` in physical space.

- `count`:
  expert/testing constructor. Build exactly `count` radial basis functions.

- `reference_spacing`:
  spacing of the underlying uniform gausslet lattice in the mapping coordinate.

- `tails`:
  number of extra half-line tail functions used to restore completeness near
  the boundary before the radial projection/refinement steps.

- `odd_even_kmax`:
  size cutoff for the small near-origin even block used in the odd-even radial
  repair.

- `xgaussians`:
  optional extra `XGaussian` primitives added before the final
  orthogonalization/X-localization step. These do not define the basis by
  themselves; they are additional near-origin directions.

Guidance
========

Most users should prefer the `rmax` constructor.
The `count` constructor is mainly for exact reproduction, testing, and
benchmarking.

See also: `build_basis`, `RadialGausslet`, `radial_quadrature`,
          `basis_diagnostics`
"""
```

### `primitives`

```julia
"""
    primitives(basis)

Return the ordered primitive function layer underlying `basis`.

The returned vector defines the row ordering used by `stencil_matrix(basis)` and
the expected ordering for primitive-space vectors, diagonals, and matrices in
the contraction helpers.

For a radial basis, these primitives are the lowest public half-line functions,
such as `HalfLineGaussian`, `XGaussian`, and mapped/distorted versions thereof.

See also: `stencil_matrix`, `contract_primitive_vector`,
          `contract_primitive_diagonal`, `contract_primitive_matrix`
"""
```

### `stencil_matrix`

```julia
"""
    stencil_matrix(basis)

Return the exact primitive-to-basis coefficient matrix for `basis`.

If `C = stencil_matrix(basis)` and `P = primitives(basis)`, then the `i`-th
basis function is defined by the column `C[:, i]`:

    basis[i](x) = Σ_μ C[μ, i] * P[μ](x)

Rows correspond to `primitives(basis)` in their returned order.
Columns correspond to basis functions in the same order as `basis[i]`.

This matrix is the basis-wide analogue of `stencil(basis[i])`. It is the main
entry point for users who want primitive-level analytic or mixed analytic /
quadrature integral workflows.

Implementations may return a dense matrix, sparse matrix, or lightweight view,
provided the above contract is respected.

See also: `stencil`, `primitives`, `contract_primitive_matrix`
"""
```

### `contract_primitive_vector`

```julia
"""
    contract_primitive_vector(basis, vμ)

Contract a primitive-space vector into basis space.

If `C = stencil_matrix(basis)`, returns

    C' * vμ

where `vμ` is ordered according to `primitives(basis)`.

See also: `contract_primitive_diagonal`, `contract_primitive_matrix`
"""
```

### `contract_primitive_diagonal`

```julia
"""
    contract_primitive_diagonal(basis, dμ)

Contract a primitive-space diagonal operator into basis space.

If `C = stencil_matrix(basis)`, returns

    C' * Diagonal(dμ) * C

where `dμ` is ordered according to `primitives(basis)`.

This is especially useful for diagonal approximations or quadrature-based
operators already expressed on the primitive layer.

See also: `contract_primitive_matrix`, `contract_primitive_vector`
"""
```

### `contract_primitive_matrix`

```julia
"""
    contract_primitive_matrix(basis, Aμν)

Contract a primitive-space operator matrix into basis space.

If `C = stencil_matrix(basis)`, this returns

    C' * Aμν * C

where `Aμν` is the operator in the primitive basis ordered according to
`primitives(basis)`.

Use this when analytic or semi-analytic integrals are most naturally computed on
the primitive layer and then transformed into the final basis. This function is
intended as the core contraction helper for expert workflows; higher-level
operator builders may call it internally.

If `Aμν` is Hermitian, the result should be Hermitian up to floating-point
roundoff.

See also: `stencil_matrix`, `primitives`, `contract_primitive_diagonal`
"""
```

### `BoundaryGausslet`

```julia
"""
    BoundaryGausslet

Callable basis function on the half line `[0, ∞)` obtained from a `HalfLineBasis`.

A `BoundaryGausslet` is the half-line analogue of an ordinary full-line
`Gausslet`. It is usually created as `hb[i]`, not by direct construction.

Boundary gausslets are localized and orthonormal on the half line. Near the
edge they differ from the corresponding full-line gausslets just enough to
restore orthogonality and low-order completeness; away from the edge they
approach the original full-line shapes.

The exact public stencil of a `BoundaryGausslet` is over the half-line
primitive layer:
- `HalfLineGaussian`
- optionally `Distorted(HalfLineGaussian(...), mapping)` when a coordinate
  mapping is present

Use
    f(x)
    value(f, x)
    direct_value(f, x)
    derivative(f, x)
    center(f)
    reference_center(f)
    moment_center(f, grid)
    integral_weight(f)
    stencil(f)

Notes
=====
- `center(f)` is the physical-space X-center.
- `reference_center(f)` is the center on the underlying uniform reference
  coordinate before any mapping.
- For modes close to the boundary, `moment_center(f, grid)` need not coincide
  with `center(f)`.

See also: `HalfLineBasis`, `RadialGausslet`, `stencil`
"""
```

### `RadialGausslet`

```julia
"""
    RadialGausslet

Callable radial basis function `χ_a(r)` from a `RadialBasis`.

A `RadialGausslet` is used to expand the reduced radial function

    u(r) = r R(r)

rather than `R(r)` itself. Each `RadialGausslet` satisfies the radial boundary
condition at the origin,

    χ_a(0) = 0,

and is orthonormal on `[0, ∞)` with respect to the ordinary `dr` measure.

When a coordinate mapping is present, the physical function is the mapped
version of an underlying half-line function on a uniform reference coordinate.
`RadialGausslet`s are usually created as `rb[i]`, not by direct construction.

The exact public stencil of a `RadialGausslet` is over the radial primitive
layer:
- `HalfLineGaussian`
- optional `XGaussian`
- optional `Distorted(...)` wrappers when a coordinate mapping is present

Use
    f(r)
    value(f, r)
    direct_value(f, r)
    derivative(f, r)
    center(f)
    reference_center(f)
    moment_center(f, grid)
    integral_weight(f)
    stencil(f)

Notes
=====
- `center(f)` is the physical radial X-center.
- `reference_center(f)` is the center on the underlying uniform reference
  coordinate before mapping.
- `moment_center(f, grid)` is especially useful near the origin, where it may
  differ from `center(f)`.
- Spherical harmonics are not part of this object; the corresponding atomic
  one-electron function is built later from `χ_a(r) / r * Y_{ℓm}(Ω)`.

See also: `RadialBasis`, `BoundaryGausslet`, `RadialBasisSpec`,
          `radial_quadrature`
"""
```

### `basis[i]`

```julia
"""
    basis[i]

Return the `i`-th basis function as a lightweight callable view.

The concrete return type depends on the basis:
- `UniformBasis`  -> `Gausslet`
- `HalfLineBasis` -> `BoundaryGausslet`
- `RadialBasis`   -> `RadialGausslet`

The returned object behaves like an ordinary function object but typically
shares cached data with its parent basis. Treat it as logically immutable.

Ordering
========
The function order returned by `basis[i]` is the canonical basis order used by:
- `centers(basis)`
- `integral_weights(basis)`
- columns of `stencil_matrix(basis)`

In particular, if `C = stencil_matrix(basis)` and `P = primitives(basis)`, then
the stencil of `basis[i]` is the `i`-th column of `C` expressed against `P`.

Examples
========
    f = rb[4]
    y = f(0.3)
    st = stencil(f)

See also: `stencil`, `primitives`, `stencil_matrix`
"""
Base.getindex(basis, i::Integer)
```

### `radial_quadrature`

```julia
"""
    radial_quadrature(basis::RadialBasis; accuracy=:high, refine=nothing, quadrature_rmax=nothing)

Build a fine quadrature grid matched to a radial basis.

The quadrature grid follows the same coordinate mapping as `basis` but is
typically much finer than the basis itself. This keeps quadrature accuracy as a
separate control knob from basis size.

Keyword arguments
=================
- `accuracy`:
  one of `:medium`, `:high`, or `:veryhigh`. The default is `:high`.
- `refine`:
  optional starting refinement factor on the uniform reference coordinate.
  `refine = 8` means roughly eight quadrature subintervals per basis spacing
  before any internal refinement.
- `quadrature_rmax`:
  optional explicit physical-space cutoff override for the quadrature grid.

With no keywords, the routine should choose a conservative cutoff and starting
resolution from the basis itself.

Use this grid for:
- one-body radial integrals (`overlap`, `kinetic`, `nuclear`, `centrifugal`)
- `moment_center(f, grid)` and related diagnostics
- radial IDA / multipole integrals

The returned grid is a true quadrature object, not a DVR basis. Its points and
weights may share the basis mapping, but they do not define the basis itself.

See also: `RadialQuadratureGrid`, `RadialBasis`, `atomic_operators`,
          `moment_center`
"""
```

### `reference_center`

```julia
"""
    reference_center(f)

Return the center of `f` on the underlying uniform reference coordinate.

For unmapped uniform functions this coincides with `center(f)`. For mapped
half-line or radial functions it returns the center before the coordinate
mapping is applied.

See also: `center`, `moment_center`
"""
```

### `moment_center`

```julia
"""
    moment_center(f, grid)

Return the first-moment center of `f` evaluated on `grid`:

    ∫ x f(x) dx / ∫ f(x) dx

or, in the radial case, the corresponding integral over `r`.

For full-line gausslets this should agree with `center(f)` to high accuracy.
For boundary and radial functions, differences near the edge or origin are
expected and are part of the method diagnostics.

See also: `center`, `reference_center`, `basis_diagnostics`
"""
```

---

## README opening draft

````markdown
# GaussletBases.jl

GaussletBases.jl is a small multilevel library for localized gausslet bases.

It supports three levels of use:

1. callable function objects
2. basis objects
3. ready-to-use operator matrices

The same library also exposes the exact primitive stencil behind each basis
function, so expert users can mix analytic primitive integrals with numerical
quadrature when needed.

## Why gausslets?

Gausslets are smooth, local, orthogonal basis functions built as compact sums
of Gaussians. For radial work, the basis is kept separate from a finer
quadrature grid, so singular Coulomb terms can be integrated accurately
without turning the basis itself into a DVR.

## Quick start

### Make and inspect a gausslet

```julia
using GaussletBases

g = Gausslet(:G10; center = 0.0, spacing = 1.0)

g(0.2)
direct_value(g, 0.2)
derivative(g, 0.2)

st = stencil(g)
coefficients(st)
primitives(st)

center(g)
integral_weight(g)
```

### Build a radial basis

```julia
map = AsinhMapping(c = 0.15, s = 0.15)

spec = RadialBasisSpec(:G10;
    rmax = 20.0,
    mapping = map,
    reference_spacing = 1.0,
    tails = 6,
    odd_even_kmax = 6,
    xgaussians = [XGaussian(alpha = 0.0936),
                  XGaussian(alpha = 0.0236)],
)

rb = build_basis(spec)

f = rb[4]

f(0.3)
center(f)
reference_center(f)

grid = radial_quadrature(rb)

moment_center(f, grid)
basis_diagnostics(rb, grid)
```

### Build radial atomic operators

```julia
ops = atomic_operators(rb, grid; Z = 2, lmax = 3)

ops.overlap
ops.kinetic
ops.nuclear

centrifugal(ops, 0)
multipole(ops, 0)
multipole(ops, 1)
```

### Work directly at the primitive level

```julia
P = primitives(rb)
C = stencil_matrix(rb)

Amunu = Matrix{Float64}(I, length(P), length(P))
A = contract_primitive_matrix(rb, Amunu)
```

This is the expert entry point for analytic or mixed analytic/quadrature
workflows.

## Design

GaussletBases.jl uses callable function objects as the main low-level interface:

```julia
f(x)
```

For a named form, use:

```julia
value(f, x)
```

To inspect the exact defining decomposition, use:

```julia
stencil(f)
```

For basis-wide primitive access, use:

```julia
primitives(basis)
stencil_matrix(basis)
```

For radial work, basis and quadrature are separate objects:

```julia
rb   = build_basis(spec)
grid = radial_quadrature(rb)
```

## Scope of v0

Included in the first public version:

- canonical gausslet families
- callable 1D gausslets
- half-line and radial bases
- coordinate mappings
- one public stencil layer
- primitive contraction helpers
- radial quadrature
- radial atomic operator builders

Not on the front page yet:

- PGDG internals
- hybrid Gaussian add-ons
- angular/radial many-body workflow layers
- HF / DMRG / MPO tooling

## Core names

```julia
GaussletFamily, Gausslet,
Gaussian, HalfLineGaussian, XGaussian, Distorted,

IdentityMapping, AsinhMapping, uofx, xofu, dudx, du2dx2,

UniformBasisSpec, HalfLineBasisSpec, RadialBasisSpec, build_basis,

stencil, primitives, stencil_matrix,
contract_primitive_vector,
contract_primitive_diagonal,
contract_primitive_matrix,

radial_quadrature, atomic_operators,

center, reference_center, moment_center,
integral_weight, basis_diagnostics
```
````

---

## Open questions still intentionally left open

These should stay open until code feedback makes them easier to judge:

1. Whether `Distorted` should be exported in v0 or only appear through `stencil(f)`.
2. Whether `XGaussian` should stay origin-centered in v0 or accept a center parameter from the start.
3. The internal storage type of `FunctionStencil` and `stencil_matrix`.
4. Which advanced family constructors belong in v0 public scope versus advanced docs only.
5. How much PGDG-specific structure should be public in the first PGDG extension.

---

## Suggested next step

Use this file as the single canonical design snapshot. Future revisions should update this file first, then adjust README/examples/docstrings to match.
