# GaussletBases.jl

GaussletBases.jl is a small library for localized gausslet function objects.

This repository currently implements a narrow callable/basis/quadrature/operator
slice of the design:

- callable primitive function objects
- callable `Gausslet` objects built from exact Gaussian stencils
- concrete uniform / half-line / radial basis specs and basis objects
- a shared basis-wide primitive layer via `primitives(basis)` and `stencil_matrix(basis)`
- primitive contraction helpers
- radial quadrature, moment centers, and modest basis diagnostics
- radial one-body matrix builders
- a two-index IDA-style radial `multipole_matrix`
- a high-level `RadialAtomicOperators` bundle
- one public stencil layer
- coordinate mappings

Exact four-index radial electron-electron APIs, PGDG, and downstream HF/DMRG
layers are deferred from v0.

## Quick Start

```julia
using GaussletBases

g = Gausslet(:G10; center = 0.0, spacing = 1.0)

x = 0.2
g(x)
value(g, x)
direct_value(g, x)
derivative(g, x)
center(g)
integral_weight(g)
```

Runnable versions of the README flow are in `examples/`.

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

The internal construction-grid spacing is a build-time control on
`build_basis`, not part of `RadialBasisSpec`. For example:

```julia
rb_fixed = build_basis(RadialBasisSpec(:G10;
    count = 6,
    mapping = map,
    reference_spacing = 1.0,
    tails = 3,
    odd_even_kmax = 2,
    xgaussians = [XGaussian(alpha = 0.2)],
); grid_h = 0.01, refine_grid_h = false)
```

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
grid = radial_quadrature(rb; refine = 24, rmax = 12.0)

quadrature_points(grid)
quadrature_weights(grid)

moment_center(rf, grid)

diag = basis_diagnostics(rb, grid)
diag.overlap_error
diag.D
```

`RadialQuadratureGrid` is separate from the basis. It follows the same mapping
as the radial basis, but is usually finer.

## Radial Operators

```julia
ops = atomic_operators(rb, grid; Z = 2.0, lmax = 2)

ops.overlap
ops.kinetic
ops.nuclear

centrifugal(ops, 2)
multipole(ops, 1)
```

The one-body matrices use the supplied quadrature grid directly. In this v0
slice, `multipole_matrix` and `multipole(ops, L)` mean the supported two-index
IDA-style radial multipole matrices, not an exact four-index electron-electron
object. `atomic_operators(...; lmax)` precomputes `centrifugal(ops, l)` for
`l = 0:lmax` and `multipole(ops, L)` for `L = 0:(2 * lmax)`.

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
