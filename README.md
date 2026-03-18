# GaussletBases.jl

GaussletBases.jl is a Julia package for working with gausslet basis functions.

Its broad foundation is the ordinary one-dimensional gausslet story:

- explicit Gaussian constructions
- mapped and distorted coordinates
- basis functions built from a visible primitive layer
- matrix construction through that shared primitive layer

Today its most mature public-facing path is **atom-centered radial calculations**. The repository also contains:

- a small explicit atomic line built on top of the radial substrate
- an experimental ordinary Cartesian branch built around mapped/hybrid
  ordinary gausslets
- an advanced 1D primitive / contraction / hierarchy line

If you are new, start with the radial path.

## Who this package is for

GaussletBases is most useful today for people who want to:

- explore ordinary gausslets as localized basis functions built from explicit Gaussian primitive layers
- explore radial gausslet bases for atoms and related model problems
- inspect the underlying Gaussian layer behind a basis
- study how primitive Gaussian layers can be contracted into more useful localized functions
- experiment with simple one-dimensional hierarchy and contraction ideas without yet committing to a full molecular workflow

It is **not** yet a complete electronic-structure workflow package. It is a basis, quadrature, and operator package, with a newer experimental line for contraction and hierarchy in 1D.

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

## Current picture

The package currently has three visible layers:

- a mature numerical radial path
- a small atomic radial-plus-angular / IDA / minimal-UHF line on top of that
  radial substrate
- an experimental ordinary Cartesian branch that is currently most promising
  in the friendlier hybrid/core-supported regime

It is **not** yet a broad electronic-structure workflow package.

## Best first path through the repository

If you are new to the package, a good first path is:

1. run the first four examples
2. read the radial quickstart
3. read the recommended atomic setup note
4. only then move on to primitive layers, contraction, and hierarchy

From the repository root:

```bash
julia --project=. examples/01_first_gausslet.jl
julia --project=. examples/02_radial_basis.jl
julia --project=. examples/03_radial_operators.jl
julia --project=. examples/04_hydrogen_ground_state.jl
```

Those four examples take you from one gausslet, to a radial basis, to radial
operators, to a real hydrogen ground-state calculation.

If you want the next atomic step after those four, run:

```bash
julia --project=. examples/15_atomic_hydrogen_ylm.jl
```

which adds the explicit `(l,m)` angular channels on top of the same radial
substrate.

## A good first atom-centered setup

For an atom of nuclear charge `Z`, the package documentation recommends starting with:

- `s = 0.2`
- `c = s / (2Z)`
- `tails = 6`
- `odd_even_kmax = 6`
- `rmax = 30.0` bohr as a good first-row starting point

In code:

```julia
using GaussletBases

Z = 2.0
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
```

This is the recommended package starting point, not the only possible one. In practice, somewhat larger values of `s`, such as `0.3` to `0.5`, can still work surprisingly well in some situations. But if you want one clean starting recipe, use `s = 0.2`.

## What you usually do next

Once you have a radial basis, the normal next steps are:

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

## Documentation map

The main entry pages are:

- [`docs/index.md`](docs/index.md)  
  The current docs map.
- [`docs/first_radial_workflow.md`](docs/first_radial_workflow.md)  
  The best first read after this README.
- [`docs/recommended_atomic_setup.md`](docs/recommended_atomic_setup.md)
  Practical starting parameters for the radial atomic path.
- [`docs/current_atomic_branch.md`](docs/current_atomic_branch.md)
  The shortest current status page for the atomic line.
- [`docs/current_ordinary_branch.md`](docs/current_ordinary_branch.md)
  The shortest current status page for the ordinary Cartesian line.
- [`docs/example_guide.md`](docs/example_guide.md)  
  The example reading/running guide.
- [`docs/intermediate_primitive_layer.md`](docs/intermediate_primitive_layer.md)  
  The first advanced note about primitive layers and contraction.

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
