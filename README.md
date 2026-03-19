# GaussletBases.jl

[![Docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://srwhite59.github.io/GaussletBases.jl/dev/)

Gausslets are localized basis functions built from short linear combinations of
Gaussians.

They are interesting because they try to keep two useful features at once:

- the analytic convenience and explicit primitive structure of Gaussian basis
  functions
- the locality and near-orthogonality that make grid-like bases attractive

GaussletBases.jl is a Julia package for building those basis functions,
constructing the quadrature grids that go with them, and forming the current
one-body and IDA-style operators on top of them.

Today the package can already do three scientifically useful things:

- a mature **radial / atomic workflow**
- a real but newer **ordinary Cartesian mapped/hybrid workflow**
- a separate **advanced research line** for contraction, hierarchy, and
  prototype PGDG-related work

If you are new, start with the radial path.

Documentation: <https://srwhite59.github.io/GaussletBases.jl/dev/>

## Who this package is for

GaussletBases is most useful today for people who want to:

- explore ordinary gausslets as localized basis functions built from explicit Gaussian primitive layers
- explore radial gausslet bases for atoms and related model problems
- inspect the underlying Gaussian layer behind a basis
- study how primitive Gaussian layers can be contracted into more useful localized functions
- experiment with simple one-dimensional hierarchy and contraction ideas without yet committing to a full molecular workflow

It is **not** yet a complete electronic-structure workflow package. It is a
basis, quadrature, and operator package, with a newer ordinary/cartesian
workflow and a separate advanced line for contraction and hierarchy in 1D.

## Installation

At present, install the package directly from GitHub:

```julia
using Pkg
Pkg.add(url = "https://github.com/srwhite59/GaussletBases.jl")
```

Then load it with:

```julia
using GaussletBases
```

## A first useful calculation

For a first atom-centered calculation, the recommended starting point is:

- `s = 0.2`
- `c = s / (2Z)`
- `rmax = 30.0` bohr

Here `s` roughly controls the overall radial spacing, while `c` roughly
controls how much resolution is concentrated near the nucleus. The fuller
setup discussion lives in the rendered manual at
[Recommended atomic setup](https://srwhite59.github.io/GaussletBases.jl/dev/howto/recommended_atomic_setup/).

The call

```julia
AsinhMapping(c = ..., s = ...)
```

uses ordinary Julia keyword arguments to the mapping constructor. There is no
special package-specific syntax there.

This first README example carries through to a real physical result:

```julia
using LinearAlgebra
using GaussletBases

Z = 1.0
s = 0.2
c = s / (2Z)

map = AsinhMapping(c = c, s = s)

rb = build_basis(RadialBasisSpec(:G10;
    rmax = 30.0,
    mapping = map,
))

grid = radial_quadrature(rb)

H = kinetic_matrix(rb, grid) +
    nuclear_matrix(rb, grid; Z = Z) +
    centrifugal_matrix(rb, grid; l = 0)

E0 = minimum(real(eigen(Hermitian(H)).values))
println("Lowest hydrogen energy: ", E0)
```

The exact nonrelativistic ground-state energy is `-0.5 Ha`, so this is the
cleanest first scientific check of the radial basis and quadrature together.

If you want the same workflow explained more slowly, with diagnostics and
setup discussion, go next to:

- [First radial workflow](https://srwhite59.github.io/GaussletBases.jl/dev/manual/)
- [Recommended atomic setup](https://srwhite59.github.io/GaussletBases.jl/dev/howto/recommended_atomic_setup/)

## What you usually do next

Once you have a first hydrogen result, the normal next steps are:

```julia
diag = basis_diagnostics(rb)
grid = radial_quadrature(rb)
ops = atomic_operators(rb, grid; Z = Z, lmax = 2)
```

Here:
- `diag` checks whether the basis is behaving well numerically
- `grid` is the separate quadrature grid used for integrals
- `ops` bundles the basic radial matrices for atomic-style work

That basis/quadrature separation is one of the central ideas in the package.

## Best first path through the repository

If you are new to the package, a good first path is:

1. read the example above
2. run the first four examples
3. read the radial quickstart
4. use the docs map to choose where to go next

From the repository root:

```bash
julia --project=. examples/01_first_gausslet.jl
julia --project=. examples/02_radial_basis.jl
julia --project=. examples/03_radial_operators.jl
julia --project=. examples/04_hydrogen_ground_state.jl
```

Then continue with:

- [Documentation home](https://srwhite59.github.io/GaussletBases.jl/dev/)
- [Manual](https://srwhite59.github.io/GaussletBases.jl/dev/manual/)
- [Example guide](https://srwhite59.github.io/GaussletBases.jl/dev/howto/example_guide/)

If you want the next atomic step after those four, run:

```bash
julia --project=. examples/15_atomic_hydrogen_ylm.jl
```

which adds the explicit `(l,m)` angular channels on top of the same radial
substrate.

## Documentation map

The main entry pages are:

- [Documentation home](https://srwhite59.github.io/GaussletBases.jl/dev/)  
  The rendered docs home page.
- [Manual](https://srwhite59.github.io/GaussletBases.jl/dev/manual/)  
  The best first read after this README.
- [Examples](https://srwhite59.github.io/GaussletBases.jl/dev/examples/)  
  The curated runnable-example entry point.
- [Reference](https://srwhite59.github.io/GaussletBases.jl/dev/reference/)  
  Curated API reference for the main exported entry points.
- [Current atomic branch](https://srwhite59.github.io/GaussletBases.jl/dev/explanations/current_atomic_branch/)
  The user-facing atomic status path.
- [Current ordinary branch](https://srwhite59.github.io/GaussletBases.jl/dev/explanations/current_ordinary_branch/)
  The user-facing ordinary-branch status path.
- [Developer Notes](https://srwhite59.github.io/GaussletBases.jl/dev/developer/)
  Lower-priority architecture and supporting-note entry points.

The narrower notes remain in the repository, but they should be read as
supporting notes after the current branch pages are clear.

## Current scope and limits

What is already useful today:

- ordinary 1D gausslet objects and explicit Gaussian constructions
- radial gausslet bases and radial one-body operators
- the first explicit one-electron `(l,m)` atomic layer
- the first static He / IDA-style interacting atomic ingredients
- the first direct / exchange / Fock / minimal-UHF atomic line in the current
  IDA model
- the experimental ordinary Cartesian mapped/hybrid one-body route in the
  friendlier regime
- explicit quadrature and diagnostics
- primitive-layer matrix construction and contraction
- the first 1D partition/hierarchy/contraction experiments

What is not yet here:

- exact non-diagonal electron-electron operators
- the fuller interacting spherical-angular atomic layer beyond the present
  one-electron `(l,m)` path
- a broad general HF workflow beyond the present minimal atomic IDA UHF kernel
- 2D or 3D nested workflows
- named Gaussian chemistry basis libraries
- larger solver layers such as DMRG or related workflows
- Python and Fortran interoperability layers

That is deliberate. The package is still settling its scientific structure before it grows into those directions.

## Acknowledgments

Development of this package was accelerated substantially with the help of
OpenAI Codex-style interactive coding assistance. Scientific direction, design
choices, and final review remained author-driven.
