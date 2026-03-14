# GaussletBases.jl

GaussletBases.jl is for numerical electronic-structure work with gausslets:
smooth, localized functions built from Gaussians. They are designed to live in
the same scientific conversation as Gaussian bases, grids, and DVR-like ideas,
while keeping locality and simple underlying building blocks in view.

The current package is especially focused on radial gausslet bases for atomic
and related problems. It lets you build mapped radial bases, explicit
quadrature grids, and radial operator matrices, while still making it possible
to inspect the Gaussian ingredients from which the basis functions are made.

The gausslet approach was introduced in 2017 and
developed further in later work on multisliced, hybrid gausslet/Gaussian, and
nested gausslet bases.

One reason gausslets attract interest is that they connect naturally to
DVR-like thinking and to the possibility of accurate, nearly diagonal
approximations to Coulomb interactions in many-body calculations. This package
does not yet provide the full exact four-index electron-electron machinery, but
it does provide the radial basis, quadrature, and operator tools that support
that direction of work.

## What Gausslets Are

Gausslets are smooth localized functions assembled from Gaussian pieces. They
are meant to give you a basis that is local in space, systematic in how it is
placed, and still closely tied to simple analytic building blocks.

In this package, gausslets appear in three closely related settings:

- individual one-dimensional gausslet functions
- half-line and radial basis sets built from them
- radial operator constructions built on top of those bases

## What This Package Is For

Typical uses of GaussletBases.jl are:

- building radial gausslet bases for atomic and related model problems
- constructing explicit quadrature grids matched to those bases
- forming radial one-body operator matrices
- experimenting with localized bases for electronic-structure and many-body
  workflows

## Why Gausslets Are Interesting

Gausslets are interesting because they try to combine some useful features of
localized basis sets and grid-based methods.

- Like localized basis functions, they are spatially local and can be adapted to
  nonuniform coordinates.
- Like grid-based approaches, they can be placed in a regular reference
  coordinate and paired naturally with quadrature.
- Because they are built from Gaussian pieces, it is possible to inspect the
  underlying Gaussian expansion directly when you need it.
- Their connection to DVR-like ideas is part of why they are attractive for
  many-body calculations, especially when one wants accurate nearly diagonal
  approximations to Coulomb interactions.

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

For most users, the natural next step is to build a small radial basis. The
example below constructs a mapped basis for the reduced radial function
`u(r) = r R(r)`, builds an explicit quadrature grid, and checks a simple
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

This is useful if you want to work directly with the Gaussian pieces underlying
the basis, rather than only with the higher-level basis functions.

## If You Want More Direct Control

If you already know the parts of the package you want to inspect, these are the
main public entry points:

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

The current v0 package does not yet include:

- an exact non-diagonal electron-electron API
- PGDG
- hybrid Gaussian extras
- Python / Fortran interop or export layers

## Examples

Short runnable examples are provided in:

- `examples/01_first_gausslet.jl`
- `examples/02_radial_basis.jl`
- `examples/03_radial_operators.jl`

## Background and References

For the scientific background on gausslets and related basis constructions, see:

- Steven R. White, *Hybrid grid/basis set discretizations of the Schrödinger equation*, J. Chem. Phys. **147**, 244102 (2017). [https://doi.org/10.1063/1.5007066](https://doi.org/10.1063/1.5007066)
  Original introduction of gausslets, including their construction and the diagonal approximation for two-electron Coulomb terms.

- Steven R. White and E. Miles Stoudenmire, *Multisliced gausslet basis sets for electronic structure*, Phys. Rev. B **99**, 081110 (2019). [https://doi.org/10.1103/PhysRevB.99.081110](https://doi.org/10.1103/PhysRevB.99.081110)
  Extends gausslets to efficient three-dimensional electronic-structure calculations using multislicing.

- Yiheng Qiu and Steven R. White, *Hybrid gausslet/Gaussian basis sets*, J. Chem. Phys. **155**, 184107 (2021). [https://doi.org/10.1063/5.0068887](https://doi.org/10.1063/5.0068887)
  Introduces hybrid bases combining gausslets with standard Gaussian functions to improve near-nuclear accuracy.

- Steven R. White and Michael J. Lindsey, *Nested gausslet basis sets*, J. Chem. Phys. **159**, 234112 (2023). [https://doi.org/10.1063/5.0180092](https://doi.org/10.1063/5.0180092)
  Introduces nested gausslet constructions, pure Gaussian distorted gausslet bases, and related extensions for larger-Z atoms and higher accuracy.
