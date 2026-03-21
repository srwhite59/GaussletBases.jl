# GaussletBases.jl

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

`GaussletBases.jl` is a Julia package for building those basis functions,
constructing the quadrature grids that go with them, and forming the current
one-body operators and integral-diagonal approximation (IDA)-style
interaction objects on top of them.

## Start here

If you are new to the package, use this short path:

1. [Manual](manual/index.md)
2. [Examples](examples/index.md)
3. [First radial workflow](tutorials/first_radial_workflow.md)
4. [Recommended atomic setup](howto/recommended_atomic_setup.md)

That is still the best front door because the radial line is the current
mature numerical workflow.

## Primary documents

The docs site is intentionally organized around this home page plus five
primary clickable documents rather than a large visible page tree:

- [Manual](manual/index.md)  
  The user-facing guide to what the package does today, where to start, and
  which workflow to follow next.
- [Algorithms](algorithms/index.md)  
  Implementation-facing basis-construction and operator-construction recipes,
  with pseudocode, code pointers, and references.
- [Examples](examples/index.md)  
  The curated runnable-example entry point.
- [Reference](reference/index.md)  
  Curated API reference built from real docstrings.
- [Developer Notes](developer/index.md)  
  Lower-priority architecture, supporting-note, and development material.

If you only want the practical reading order, go straight to the
[Manual](manual/index.md). If you already know what object or function you want
to call, use the [Reference](reference/index.md).

## Manual first, Reference second

The Manual, Algorithms, Examples, and Reference now divide the work in the
standard Julia-package way:

- the [Manual](manual/index.md) explains workflows, branch status, and what to
  run first
- the [Algorithms](algorithms/index.md) record the exact construction order for
  important basis and operator paths
- the [Examples](examples/index.md) show which scripts to run and in what order
- the [Reference](reference/index.md) answers API questions about the main
  exported entry points
- the [Developer Notes](developer/index.md) preserve lower-priority design and
  history material without competing with the user docs

That keeps the visible surface small while still preserving the full
development record underneath.

## Current scope

Today the package is strongest in three areas:

- a mature radial / atomic workflow
- a real but newer ordinary Cartesian mapped/hybrid workflow
- a separate advanced/research line for contraction, hierarchy, and supporting
  PGDG-related work

It is not yet a broad electronic-structure workflow package. The package is a
basis, quadrature, and operator package, with narrow solver-facing layers only
where they are already scientifically justified.

## Build the docs locally

From the repository root:

```bash
julia --project=docs -e 'using Pkg; Pkg.instantiate()'
julia --project=docs docs/make.jl
```
