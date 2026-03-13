# Gausslets.jl

Gausslets.jl is a small library for localized gausslet function objects.

This repository currently implements only the first narrow slice of the design:

- callable primitive function objects
- callable `Gausslet` objects built from exact Gaussian stencils
- one public stencil layer
- coordinate mappings

Radial basis construction, quadrature, diagnostics, and operator bundles are
planned but are not implemented in this pass.

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

## Coordinate Mappings

The mapping object is the forward coordinate map:

```julia
map = AsinhMapping(c = 0.15, s = 0.15)

u = map(3.0)
x = xofu(map, u)

uofx(map, 3.0)
dudx(map, 3.0)
du2dx2(map, 3.0)
```

`AsinhMapping` in this implementation includes the constant-density tail term
controlled by `tail_spacing`.
