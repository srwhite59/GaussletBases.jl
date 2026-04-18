# GaussletBases.jl

[![Docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://srwhite59.github.io/GaussletBases.jl/dev/)

Gausslets are localized, orthogonal basis functions constructed from short
linear combinations of Gaussians.

They were developed to combine several properties that are rarely available in
one basis at the same time: orthogonality, locality, smoothness, systematic
improvability by refining the spacing, and a moment structure that supports
diagonal or integral-diagonal approximations for the Coulomb interaction.

The practical attraction of gausslets is that they can produce compact
second-quantized Hamiltonians. Compared with standard Gaussian orbital bases,
they are much more local and naturally orthogonal; compared with ordinary grid
or DVR constructions, they retain an explicit Gaussian primitive structure and
fit naturally with variable-resolution mappings and hybrid Gaussian
augmentations near nuclei.

GaussletBases.jl is a Julia package for building those basis functions,
constructing the quadrature grids that go with them, and forming the current
one-body and integral-diagonal approximation (IDA)-style operators on top of
them.

Today the package has four real current surfaces:

- a mature **radial / atomic workflow**
- a real **ordinary Cartesian workflow** with exact basis-to-basis Cartesian
  `cross_overlap`, `basis_projector`, and `transfer_orbitals` primitives
- a real **nested Cartesian / diatomic workflow surface** with one-center
  atomic fixed-block routes, bond-aligned diatomic nested source/fixed-block
  construction, geometry diagnostics, and geometry payload export
- an experimental **angular and advanced research track** for injected angular
  bases, contraction, hierarchy, and prototype PGDG-related work

If you are new, start with the radial path. It is also the most recent
published line in the repo: the radial gausslet paper is now posted on arXiv,
<https://doi.org/10.48550/arXiv.2603.22646>.

Documentation: <https://srwhite59.github.io/GaussletBases.jl/dev/>

The rendered docs are now organized into a small set of primary sections:
Manual, Algorithms, Examples, Reference, and Developer Notes.

## Who this package is for

GaussletBases is most useful today for people who want to:

- explore ordinary gausslets as localized basis functions built from explicit Gaussian primitive layers
- explore radial gausslet bases for atoms and related model problems
  The most recent posted paper on that line is
  <https://doi.org/10.48550/arXiv.2603.22646>.
- experiment with the manuscript-facing angular injected-basis research track
  while treating it as experimental rather than stable workflow surface
- inspect the underlying Gaussian layer behind a basis
- study how primitive Gaussian layers can be contracted into more useful localized functions
- experiment with simple one-dimensional hierarchy and contraction ideas without yet committing to a full molecular workflow

It is **not** yet a complete electronic-structure workflow package. It is a
basis, quadrature, and operator package, with a mature radial/atomic line, a
real newer ordinary/cartesian line, a real nested Cartesian/diatomic surface,
and narrower experimental solver-facing or hierarchy-facing research tracks.

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
controls how much resolution is concentrated near the nucleus. `rmax` is the
center of the last retained radial gausslet, so it is the main user-facing
scientific extent. The library chooses the separate internal build and
quadrature extents automatically. The fuller setup discussion lives in the
rendered manual at
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

In the normal user workflow, `radial_quadrature(rb)` should be the default.
The expert keyword `quadrature_rmax` still exists for compatibility, but the
library normally owns the internal quadrature extent.

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

From a fresh checkout, instantiate once before running examples:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

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
- [Algorithms](https://srwhite59.github.io/GaussletBases.jl/dev/algorithms/)  
  Basis-construction and operator-construction recipes with pseudocode, code
  pointers, and paper references.
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

## Top-level repo docs

At the repository root, the current documentation authority is now split
deliberately:

- `README.md`: onboarding and first trustworthy repo overview
- `DESIGN.md`: stable repo-wide design/contracts note and doc-authority map
- `STATUS.md`: current capability / trust matrix
- `ROADMAP.md`: strategic next-pressure note, not a schedule

For branch-specific current status, use:

- [Current atomic branch](https://srwhite59.github.io/GaussletBases.jl/dev/explanations/current_atomic_branch/)
- [Current ordinary branch](https://srwhite59.github.io/GaussletBases.jl/dev/explanations/current_ordinary_branch/)

For path-specific construction details, use the rendered
[Algorithms](https://srwhite59.github.io/GaussletBases.jl/dev/algorithms/)
pages rather than older flat note-history files.

## Current scope and limits

What is already useful today:

- ordinary 1D gausslet objects and explicit Gaussian constructions
- the mature radial/atomic line: radial bases, radial one-body operators, the
  current `(l,m)` atomic layer, static He / IDA-style interacting ingredients,
  direct / exchange / Fock helpers, a minimal UHF kernel, and dense/sliced
  export for the current density-density model
- exact Cartesian basis-to-basis cross overlap, projector, and transfer on the
  current working representation families
- one-center nested Cartesian shell-sequence, fixed-block, and QW consumer
  routes
- bond-aligned diatomic nested fixed-source / fixed-block / diagnostics /
  geometry payload support, plus the current bond-aligned diatomic QW workflow
- the newer ordinary Cartesian mapped/hybrid and Qiu-White residual-Gaussian
  routes in their supported regimes
- the experimental homonuclear chain / square-lattice nested producer routes
- the experimental angular research track: shell-local injected angular basis
  construction, shell-to-atom assembly, one-electron benchmark, HF-style
  benchmark, small-ED benchmark, and a direct in-memory HFDMRG payload path
- explicit quadrature and diagnostics
- primitive-layer matrix construction and contraction
- the current contraction / hierarchy research line

What is not yet here:

- broad exact four-index electron-electron workflows
- a broad stabilized angular atomic workflow beyond the current experimental
  benchmark ladder and direct HF payload handoff
- a broad general HF workflow beyond the present minimal atomic IDA UHF kernel
- a broad stabilized molecule-scale nested workflow beyond the current
  one-center and bond-aligned diatomic surfaces
- named Gaussian chemistry basis libraries
- true many-body DMRG or related workflow layers
- Python and Fortran interoperability layers

That is deliberate. The package is still settling its scientific structure before it grows into those directions.

## Acknowledgments

This repository draws on the author's substantial collection of earlier Julia
modules. Work on the repository itself, as a complete rewrite and
consolidation, began on March 12, 2026. It was made public on March 19, 2026,
with development accelerated substantially by OpenAI Codex-style interactive
coding assistance. Scientific direction, design choices, and final review
remain author-driven.
