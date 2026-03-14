# GaussletBases.jl

GaussletBases.jl is a Julia package for working with gausslet basis functions.
Gausslets are smooth, localized functions built from Gaussians. They are meant
to combine some of the appealing features of localized basis sets and
grid-like methods: locality, systematic placement, and a direct connection to
simple underlying building blocks.

This package is currently focused on one-dimensional, half-line, and reduced
radial gausslet bases. In practice, that makes it most useful for atomic and
related electronic-structure experiments where one wants localized radial bases,
explicit quadrature, and direct access to operator matrices.

If you know Gaussian basis sets, grids, DVRs, Hartree-Fock, DMRG, PySCF, or
block2, the rough idea is that gausslets sit in that neighborhood while keeping
locality front and center.

## Why Gausslets Are Interesting

Gausslets are interesting because they try to give you some of the good parts
of two different worlds at once.

- Like localized basis functions, they are spatially local and can be adapted to
  nonuniform coordinates.
- Like grid-based approaches, they can be placed in a regular reference
  coordinate and paired naturally with quadrature.
- Because they are built from Gaussian pieces, there is also a clean analytic
  layer underneath the basis functions themselves.

This is especially appealing in radial electronic-structure work, where one
often wants basis functions that are local, well behaved near the origin, and
compatible with mapped coordinates and operator construction.

## What v0 Currently Includes

The current package lets you:

- make and evaluate individual gausslets
- build uniform, half-line, and radial gausslet bases
- inspect the exact Gaussian expansion used to define a function
- build explicit radial quadrature grids and basic diagnostics
- construct radial overlap, kinetic, nuclear, centrifugal, and IDA-style
  multipole matrices
- bundle those radial operators into a `RadialAtomicOperators` object

The package does not yet include an exact non-diagonal electron-electron API,
PGDG, hybrid Gaussian extras, or Python / Fortran interop layers.

## First Thing To Try

If you just want to see what a gausslet looks like in code, start by making one
gausslet and evaluating it at a point:

```julia
using GaussletBases

g = Gausslet(:G10; center = 0.0, spacing = 1.0)
x = 0.2

g(x)
value(g, x)
center(g)
integral_weight(g)
```

You can also inspect the exact Gaussian expansion used to build that function:

```julia
st = stencil(g)

st(x)
coefficients(st)
primitives(st)
```

Runnable versions of this first example live in `examples/01_first_gausslet.jl`.

## A Small Radial Example

The current package is most useful once you move to radial bases. The example
below builds a small mapped radial basis for the reduced radial function
`u(r) = r R(r)`, constructs an explicit quadrature grid, and checks a basic
diagnostic:

```julia
using GaussletBases

map = AsinhMapping(c = 0.15, s = 0.15)
spec = RadialBasisSpec(:G10;
    count = 6,
    mapping = map,
    reference_spacing = 1.0,
    tails = 3,
    odd_even_kmax = 2,
    xgaussians = [XGaussian(alpha = 0.2)],
)

rb = build_basis(spec)
f = rb[2]
grid = radial_quadrature(rb; refine = 24, rmax = 12.0)
diag = basis_diagnostics(rb, grid)

f(0.2)
reference_center(f)
center(f)
moment_center(f, grid)
diag.overlap_error
```

Two practical points are worth noting here.

- `radial_quadrature` takes an explicit `rmax`. That outer radius is a physical
  modeling choice, not something the package guesses for you.
- `build_basis(spec; grid_h = ..., refine_grid_h = ...)` lets you control the
  internal construction grid if you need to, but most users can start with the
  default adaptive behavior.

Runnable version: `examples/02_radial_basis.jl`.

## Construct Radial One-Body Operators

Once you have a radial basis and an explicit quadrature grid, you can build the
basic radial one-body operators used in atomic-style calculations:

```julia
using GaussletBases

map = AsinhMapping(c = 0.15, s = 0.15)
rb = build_basis(RadialBasisSpec(:G10;
    count = 6,
    mapping = map,
    reference_spacing = 1.0,
    tails = 3,
    odd_even_kmax = 2,
    xgaussians = [XGaussian(alpha = 0.2)],
))

grid = radial_quadrature(rb; refine = 24, rmax = 12.0)
ops = atomic_operators(rb, grid; Z = 2.0, lmax = 2)

ops.overlap
ops.kinetic
ops.nuclear

centrifugal(ops, 2)
multipole(ops, 1)
```

In this v0 package, `multipole_matrix` and `multipole(ops, L)` mean the
two-index radial multipole matrix used in the supported IDA-style workflow.
They are not an exact four-index electron-electron object.

Runnable version: `examples/03_radial_operators.jl`.

## A Note On Terminology

One package-specific word appears often enough to explain once up front:

- `stencil(f)` returns the exact Gaussian expansion used to build the function
  `f`. In the code this expansion is called the stencil.

For an entire basis, you can also inspect the common set of underlying building
blocks and the coefficients that combine them into each basis function:

```julia
P = primitives(rb)
C = stencil_matrix(rb)

sum(C[mu, 2] * P[mu](0.2) for mu in eachindex(P))
rb[2](0.2)
```

This is mainly expert-facing machinery, but it is useful if you want direct
access to the Gaussian pieces underlying the basis.

## Expert-Facing API Pointers

If you already know what you want to inspect, the main public entry points are:

- individual functions: `Gausslet`, `Gaussian`, `HalfLineGaussian`,
  `XGaussian`, `Distorted`
- basis recipes and bases: `UniformBasisSpec`, `HalfLineBasisSpec`,
  `RadialBasisSpec`, `build_basis`, `basis[i]`, `centers`,
  `reference_centers`, `integral_weights`
- basis-wide primitive layer: `primitives(basis)`, `stencil_matrix(basis)`,
  `contract_primitive_vector`, `contract_primitive_diagonal`,
  `contract_primitive_matrix`
- mappings: `IdentityMapping`, `AsinhMapping`, `uofx`, `xofu`, `dudx`,
  `du2dx2`
- quadrature and diagnostics: `radial_quadrature`, `moment_center`,
  `basis_diagnostics`
- radial operators: `overlap_matrix`, `kinetic_matrix`, `nuclear_matrix`,
  `centrifugal_matrix`, `multipole_matrix`, `atomic_operators`

## What Is Deferred From v0

The current release candidate does not yet include:

- an exact non-diagonal electron-electron API
- PGDG
- hybrid Gaussian extras
- Python / Fortran interop or export layers

## Examples

Short runnable examples are provided in:

- `examples/01_first_gausslet.jl`
- `examples/02_radial_basis.jl`
- `examples/03_radial_operators.jl`
