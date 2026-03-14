# GaussletBases.jl

GaussletBases.jl is a Julia package for radial and one-dimensional gausslet
bases in numerical electronic-structure calculations.

Gausslets are orthonormal, infinitely smooth, highly localized functions built
as sums of Gaussians. When laid out on a regular reference grid, they also
satisfy moment conditions that make them integrate low-order polynomials like
delta functions centered on the grid points. That combination—orthogonality,
smoothness, locality, quadrature behavior, and an explicit Gaussian
representation—is what makes gausslets distinct. The key benefit is that they
allow the use of accurate two-index two-electron interactions, rather than the
usual expensive four-index form.

For 3D, ordinary gausslets use products of 1D functions along the coordinate axes, f(x) g(y) h(z).
The most recent development in gausslets is the introduction of radial
gausslets for atomic and related problems, and this package is built around
that development.  For ordinary gausslet bases, most integration for Hamiltonian
terms can be done analytically. For radial gausslets, an explicit quadrature
grid is currently a better choice, and GaussletBases.jl lets you construct
coordinate-mapped radial gausslet bases, build explicit radial quadrature
grids, and form radial operator matrices. It also keeps the exact Gaussian
expansion of each function accessible, so you can inspect the primitive
Gaussian layer directly when you need it.  In addition it allows the
construction and use of ordinary gausslets, but currently without Hamiltonian
generation.

Gausslets were introduced in 2017 and developed further in later work on
multisliced gausslets, hybrid gausslet/Gaussian bases, and nested gausslet
bases. A short list of references appears at the end of this README.

## What Gausslets Are

Gausslets are wavelet-like basis functions represented exactly as linear
combinations of Gaussians. In practice, that means you get several useful
things at once:

- an orthonormal basis with strong spatial locality
- infinite smoothness, rather than piecewise or cusp-like behavior
- moment properties that support delta-function-like quadrature behavior
- direct access to the underlying Gaussian expansion

In this package, gausslets appear in three closely related settings:

- individual one-dimensional gausslet functions
- half-line and radial basis sets built from them
- radial quadrature and operator constructions built on top of those bases

## What This Package Is For

Typical uses of GaussletBases.jl are:

- building radial gausslet bases for atomic and related model problems
- constructing explicit quadrature grids matched to those bases
- forming radial overlap, kinetic, nuclear, centrifugal, and multipole matrices
- inspecting the Gaussian stencil behind a basis function or an entire basis
- experimenting with localized bases for electronic-structure and many-body
  workflows

## Why Gausslets Are Interesting

Gausslets are interesting because they combine properties that usually do not
come together in one basis.

- They are orthonormal and highly local at the same time.
- They are infinitely smooth.
- Their moment conditions make quadrature and diagonal approximations natural.
- Because they are built from Gaussian pieces, there is an explicit analytic
  layer underneath the basis.
- In many-body settings, these features make accurate diagonal
  approximations to Coulomb interactions possible.

This is especially useful in electronic-structure work, where one wants
localized basis functions, controlled mapped coordinates, explicit quadrature,
and direct access to operator matrices.

## What v0 Currently Includes

The current package lets you:

- make and evaluate individual gausslets
- build uniform, half-line, and radial gausslet bases
- inspect the exact Gaussian expansion used to define a function
- build explicit radial quadrature grids and basic diagnostics
- construct radial overlap, kinetic, nuclear, centrifugal, and IDA-style
  two-index multipole matrices
- bundle those radial operators into a `RadialAtomicOperators` object

The present v0 scope is intentionally centered on the radial basis,
quadrature, and operator layer.

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

In this package, `multipole_matrix` and `multipole(ops, L)` refer to the
supported two-index radial multipole matrices used in the IDA-style workflow.

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

Features deferred from the current v0 package include:

- pure Gaussian distorted gausslet (PGDG) extensions
- hybrid gausslet/Gaussian extras beyond the present radial layer
- Python / Fortran interop or export layers

## Examples

Short runnable examples are provided in:

- `examples/01_first_gausslet.jl`
- `examples/02_radial_basis.jl`
- `examples/03_radial_operators.jl`

## Background and References

For the scientific background on gausslets and related basis constructions, see:

- Steven R. White, *Hybrid grid/basis set discretizations of the Schrödinger equation*, J. Chem. Phys. **147**, 244102 (2017). [https://doi.org/10.1063/1.5007066](https://doi.org/10.1063/1.5007066)
  Original introduction of gausslets: orthogonal, infinitely smooth, local, polynomially complete functions built from sums of Gaussians, together with diagonal approximations for two-electron Coulomb terms.

- Steven R. White and E. Miles Stoudenmire, *Multisliced gausslet basis sets for electronic structure*, Phys. Rev. B **99**, 081110 (2019). [https://doi.org/10.1103/PhysRevB.99.081110](https://doi.org/10.1103/PhysRevB.99.081110)
  Extends gausslets to efficient three-dimensional electronic-structure calculations using multislicing.

- Yiheng Qiu and Steven R. White, *Hybrid gausslet/Gaussian basis sets*, J. Chem. Phys. **155**, 184107 (2021). [https://doi.org/10.1063/5.0068887](https://doi.org/10.1063/5.0068887)
  Introduces hybrid bases combining gausslets with standard Gaussian functions to improve near-nuclear accuracy while preserving orthonormality and the diagonal two-electron structure.

- Steven R. White and Michael J. Lindsey, *Nested gausslet basis sets*, J. Chem. Phys. **159**, 234112 (2023). [https://doi.org/10.1063/5.0180092](https://doi.org/10.1063/5.0180092)
  Introduces nested gausslet constructions and related extensions, and clarifies the relation between completeness, orthogonality, zero-moment conditions, and diagonal structure.
