# GaussletBases.jl

GaussletBases.jl is a Julia package for building and using gausslet basis functions.

Gausslets are localized basis functions built from sums of Gaussians. They are meant to combine some of the useful features of ordinary basis sets and grid methods: they are smooth and orthonormal like a basis, but they are also local in space and can be paired with explicit quadrature grids when you need accurate integrals.

The current package is especially focused on 1D, half-line, and radial gausslet bases. Right now it is most useful for people who want to explore radial basis sets for atoms and related model problems in a way that keeps the basis compact and the radial operator construction explicit.

## Installation

At the moment, install GaussletBases directly from GitHub:

```julia
using Pkg
Pkg.add(url = "https://github.com/srwhite59/GaussletBases.jl")
```

Then load it in Julia with:

```julia
using GaussletBases
```

## What this package currently does

Version 0 includes:

* ordinary 1D gausslets
* uniform, half-line, and radial gausslet bases
* explicit coordinate mappings for variable radial resolution
* an explicit radial quadrature grid, separate from the basis itself
* basis diagnostics such as overlap error and center-mismatch measures
* radial one-body matrices (`overlap`, `kinetic`, `nuclear`, `centrifugal`)
* the current two-index IDA-style radial multipole matrices
* a high-level `RadialAtomicOperators` bundle for radial operator work

The package is intentionally still small. It is a basis/quadrature/operator library, not yet a full electronic-structure workflow package.

## First thing to try

Here is the smallest possible example: make one gausslet, evaluate it, and inspect how it is built from Gaussians.

```julia
using GaussletBases

g = Gausslet(:G10; center = 0.0, spacing = 1.0)

x = 0.2
g(x)
value(g, x)
derivative(g, x)
center(g)
integral_weight(g)
```

If you want to see the exact Gaussian expansion used to build `g`, use:

```julia
st = stencil(g)

coefficients(st)
primitives(st)
st(x)
```

### A note on the word `stencil`

The word `stencil` is not standard chemistry language, so here is what it means in this package.

`stencil(f)` returns the exact Gaussian expansion used to build the function `f`.

You do not need this for everyday use. It is mainly there when you want to inspect the underlying Gaussian building blocks or work directly with the common Gaussian layer shared by a whole basis.

Runnable versions of the examples in this README are in the `examples/` directory.

## Recommended atomic starting point

For atom-centered radial calculations, a good package-level starting point is:

* `s = 0.2`
* `c = s / (2Z)`
* `tail_spacing = 10.0`
* `tails = 6`
* `odd_even_kmax = 6`
* start with `rmax = 30.0` bohr for first-row atoms

In code, that looks like this:

```julia
Z = 2.0
s = 0.2

map = AsinhMapping(c = s / (2Z), s = s)

spec = RadialBasisSpec(:G10;
    rmax = 30.0,
    mapping = map,
    reference_spacing = 1.0,
    tails = 6,
    odd_even_kmax = 6,
    xgaussians = XGaussian[],
)

rb = build_basis(spec)
```

This is not the only reasonable choice. Larger values of `s`, such as `0.3` or even `0.5`, can still work surprisingly well in some situations. But `s = 0.2` is the standard starting recommendation in the package documentation.

The `tail_spacing = 10.0` part is built into the default `AsinhMapping` constructor. It keeps the outer spacing from growing without bound.

For more discussion of how to choose these parameters, see [`docs/recommended_atomic_setup.md`](docs/recommended_atomic_setup.md).

## A first radial basis example

Once you have a radial basis, you can inspect one basis function and build a matching quadrature grid:

```julia
f = rb[4]

f(0.3)
reference_center(f)
center(f)

grid = radial_quadrature(rb)

moment_center(f, grid)

diag = basis_diagnostics(rb)
diag.overlap_error
diag.D
```

Here `basis_diagnostics(rb)` chooses its own conservative integration grid internally, so it is a good default check while you are learning the basis.

The basis and the quadrature grid are separate on purpose. The basis is the compact variational object you expand in. The quadrature grid is the finer grid used to evaluate radial integrals accurately.

For most users, `radial_quadrature(rb)` is the right starting point. It chooses a conservative cutoff and starting resolution automatically. The optional `quadrature_rmax` and `refine` keywords are there mainly for manual control and benchmarking.

## Hydrogen ground-state example

A very good first scientific check is hydrogen. It is simple enough to understand immediately, but it already tests whether the radial basis and one-body operators are working together correctly.

The repository includes a runnable hydrogen example in:

* `examples/04_hydrogen_ground_state.jl`

The basic workflow is:

```julia
using LinearAlgebra
using GaussletBases

Z = 1.0
s = 0.2

map = AsinhMapping(c = s / (2Z), s = s)

rb = build_basis(RadialBasisSpec(:G10;
    rmax = 30.0,
    mapping = map,
    reference_spacing = 1.0,
    tails = 6,
    odd_even_kmax = 6,
    xgaussians = XGaussian[],
))

grid = radial_quadrature(rb)

S = overlap_matrix(rb, grid)
H = kinetic_matrix(rb, grid) +
    nuclear_matrix(rb, grid; Z = Z) +
    centrifugal_matrix(rb, grid; l = 0)

eig = eigen(Hermitian(H), Hermitian(S))
E0 = minimum(real(eig.values))
E0
```

The exact nonrelativistic ground-state energy of hydrogen is `-0.5 Ha`, so this is an easy place to begin checking basis and quadrature choices.

## Radial one-body operators

For one-electron or mean-field work, the package can build radial one-body matrices directly from a basis and an explicit quadrature grid:

```julia
S = overlap_matrix(rb, grid)
T = kinetic_matrix(rb, grid)
V = nuclear_matrix(rb, grid; Z = 2.0)
C2 = centrifugal_matrix(rb, grid; l = 2)
```

These matrices use the supplied quadrature grid directly. The package does not silently create a second hidden grid for them.

## Radial operator bundle

If you want a higher-level container for radial atomic work, use `atomic_operators`:

```julia
ops = atomic_operators(rb, grid; Z = 2.0, lmax = 2)

ops.overlap
ops.kinetic
ops.nuclear

centrifugal(ops, 2)
multipole(ops, 1)
```

In the current v0 package, `multipole_matrix` and `multipole(ops, L)` mean the supported two-index IDA-style radial multipole matrices. They are not exact four-index electron-electron tensors.

With `atomic_operators(...; lmax = 2)`, the bundle precomputes:

* `centrifugal(ops, l)` for `l = 0:2`
* `multipole(ops, L)` for `L = 0:4`

## Looking under the hood

Most users do not need this at first, but the package also lets you inspect the common Gaussian building blocks behind an entire basis.

```julia
P = primitives(rb)
C = stencil_matrix(rb)

x = 0.3
sum(C[mu, 4] * P[mu](x) for mu in eachindex(P))
rb[4](x)
```

Here:

* `P` is the shared list of underlying Gaussian-type building blocks
* `C` is the coefficient matrix that combines them into each basis function

This is useful if you want to work directly at the common Gaussian layer or contract primitive-space data into the final basis.

## Examples included in the repository

The repository currently includes these runnable examples:

* `examples/01_first_gausslet.jl`
* `examples/02_radial_basis.jl`
* `examples/03_radial_operators.jl`
* `examples/04_hydrogen_ground_state.jl`

If you are new to the package, start with them in that order.

## More documentation in this repository

If you want a little more guidance than the README can comfortably hold, these two pages are the best next stops:

* [`docs/recommended_atomic_setup.md`](docs/recommended_atomic_setup.md)
* [`docs/first_radial_workflow.md`](docs/first_radial_workflow.md)
* [`docs/terminology.md`](docs/terminology.md)

The first page is about parameter choices. The second is about the normal workflow of a first radial calculation. The third explains the small amount of package-specific vocabulary used in the code and docs.

## What is not in v0 yet

The current release intentionally does not yet include:

* an exact non-diagonal electron-electron API
* PGDG
* hybrid Gaussian add-ons
* Python or Fortran interop layers
* full HF / DMRG workflow tooling

Those are natural next extensions, but they are not the goal of the present release.
