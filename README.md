# Gausslets.jl

Gausslets.jl is a small library for localized gausslet function objects.

This repository currently implements only the first narrow slice of the design:

- callable primitive function objects
- callable `Gausslet` objects built from exact Gaussian stencils
- concrete uniform / half-line / radial basis specs and basis objects
- a shared basis-wide primitive layer via `primitives(basis)` and `stencil_matrix(basis)`
- primitive contraction helpers
- radial quadrature, moment centers, and modest basis diagnostics
- one public stencil layer
- coordinate mappings

One-electron matrices, multipole matrices, and the high-level atomic operator
bundle are planned but are not implemented in this pass.

## Quick Start

```julia
using Gausslets

g = Gausslet(:G10; center = 0.0, spacing = 1.0)

x = 0.2
g(x)
value(g, x)
direct_value(g, x)
derivative(g, x)
center(g)
integral_weight(g)
```

## Inspect The Exact Stencil

`stencil(g)` returns the exact public Gaussian stencil of a gausslet.

```julia
st = stencil(g)

coefficients(st)
primitives(st)
st(x)
```

The stencil itself is callable and evaluates the same function as `g`.

## Build A Uniform Basis

```julia
ub = build_basis(UniformBasisSpec(:G10;
    xmin = -2.0,
    xmax = 2.0,
    spacing = 1.0,
))

length(ub)
f = ub[3]

f isa Gausslet
center(f)
centers(ub)
```

## Build A Radial Basis

```julia
map = AsinhMapping(c = 0.15, s = 0.15)

rb = build_basis(RadialBasisSpec(:G10;
    count = 6,
    mapping = map,
    reference_spacing = 1.0,
    tails = 3,
    odd_even_kmax = 2,
    xgaussians = [XGaussian(alpha = 0.2)],
))

length(rb)
rf = rb[2]

rf(0.2)
reference_center(rf)
center(rf)
stencil(rf)
```

The `count` constructor builds an exact number of radial basis functions. The
`rmax` constructor chooses enough functions to cover approximately up to the
requested range after mapping.

## Inspect The Shared Primitive Layer

```julia
P = primitives(rb)
C = stencil_matrix(rb)

x = 0.2
sum(C[mu, 2] * P[mu](x) for mu in eachindex(P))
rf(x)
```

Rows of `C` follow `primitives(rb)`. Columns of `C` follow `rb[i]`.

## Contract Primitive-Space Data

```julia
using LinearAlgebra

Amunu = Matrix{Float64}(I, length(P), length(P))
A = contract_primitive_matrix(rb, Amunu)

size(A)
```

The contraction helpers use the exact ordering returned by `primitives(basis)`.

## Radial Quadrature And Diagnostics

```julia
grid = radial_quadrature(rb; refine = 8)

quadrature_points(grid)
quadrature_weights(grid)

moment_center(rf, grid)

diag = basis_diagnostics(rb, grid)
diag.overlap_error
diag.D
```

`RadialQuadratureGrid` is separate from the basis. It follows the same mapping
as the radial basis, but is usually finer.

## Coordinate Mappings

The mapping object is the forward coordinate map:

```julia
map = AsinhMapping(c = 0.15, s = 0.15)
# Equivalent:
# map = AsinhMapping(a = 1.0, s = 0.15)

u = map(3.0)
x = xofu(map, u)

uofx(map, 3.0)
dudx(map, 3.0)
du2dx2(map, 3.0)
```

`AsinhMapping` in this implementation includes the constant-density tail term
controlled by `tail_spacing`. Its public constructor meanings are:

- `a` is the direct parameter in `asinh(x / a) / s`
- `c` is the derived near-origin control with `c = a * s`
