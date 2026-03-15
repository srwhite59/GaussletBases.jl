# A first radial workflow

This page is meant to answer a simple question:

**If I want to do a first atom-centered radial calculation with GaussletBases, what do I actually do?**

The package is built in layers, but you do not need to learn all of them at once. For a first calculation, the usual sequence is:

1. choose a radial mapping and basis recipe
2. build the radial basis
3. build a matching quadrature grid
4. check a few diagnostics
5. build one-body operators
6. solve a simple test problem such as hydrogen

That is the whole basic workflow.

## 1. Choose a radial mapping and basis recipe

For an atom of nuclear charge `Z`, the package documentation recommends starting with:

- `s = 0.2`
- `c = s / (2Z)`
- `tail_spacing = 10.0`
- `tails = 6`
- `odd_even_kmax = 6`
- `reference_spacing = 1.0`
- `rmax = 30.0` bohr for a first-row atom

The first two numbers set the coordinate mapping, which controls how resolution is distributed in radius. The next two control the small near-origin and boundary corrections used in the radial basis construction. `rmax` is the physical outer radius you want to represent.

In code:

```julia
using GaussletBases

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

For a first serious calculation, this is a much better starting point than the tiny toy examples used in quick smoke tests.

## 2. Build a matching quadrature grid

The basis and the quadrature grid are separate in this package, on purpose.

The basis is the compact set of functions you expand in. The quadrature grid is the finer grid you use to evaluate radial integrals accurately.

```julia
grid = radial_quadrature(rb)
```

Here:

* `rb` is the radial basis
* the package chooses a conservative quadrature cutoff automatically
* the package also chooses a default starting resolution automatically

In many grid-based methods, the basis and the integration grid are tied together. In GaussletBases they are intentionally separate. That is one of the main ideas behind the package.

If you later want manual control, the package also offers `accuracy = :medium`, `:high`, and `:veryhigh`, with `:high` as the default. The advanced keywords underneath that are `quadrature_rmax` and `refine`.

## 3. Check the basis diagnostics

Before doing anything more ambitious, it is a good habit to inspect the diagnostics.

```julia
diag = basis_diagnostics(rb)

diag.overlap_error
diag.D
```

Here `basis_diagnostics(rb)` chooses its own conservative integration grid internally. That is usually the best first diagnostic path.

Two especially useful quantities are:

* `diag.overlap_error`: how far the numerically evaluated overlap matrix is from exact orthonormality
* `diag.D`: an aggregate measure of how much the nominal basis centers differ from the first-moment centers computed on the quadrature grid

You do not need `D` to be “mysteriously small” in some abstract sense. The point is simpler: if `D` looks unexpectedly large or unstable when you change parameters, that is a sign to inspect the setup more carefully.

## 4. Inspect one basis function if you want to understand the basis

You do not have to do this every time, but it is often helpful while learning.

```julia
f = rb[4]

f(0.3)
reference_center(f)
center(f)
moment_center(f, grid)
```

A useful mental picture is:

* `reference_center(f)` is the center before the coordinate mapping is applied
* `center(f)` is the mapped physical-space center
* `moment_center(f, grid)` is the first-moment center computed numerically on the quadrature grid

Near the origin, those do not have to agree perfectly, and that is part of why the diagnostics matter.

## 5. Build one-body operators

For a one-electron radial problem, the Hamiltonian comes from the overlap, kinetic, nuclear, and centrifugal pieces.

```julia
S = overlap_matrix(rb, grid)
T = kinetic_matrix(rb, grid)
V = nuclear_matrix(rb, grid; Z = Z)
C0 = centrifugal_matrix(rb, grid; l = 0)

H = T + V + C0
```

At this stage, everything is still one-electron and radial. That makes it a very clean place to learn how the package works.

## 6. Solve hydrogen as a first scientific test

Hydrogen is the best first test problem because it is simple, physically familiar, and free of electron-electron complications.

```julia
using LinearAlgebra

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

The exact nonrelativistic ground-state energy is `-0.5 Ha`, so this gives you a direct and easy-to-understand accuracy check.

A runnable version of this calculation is in:

* `examples/04_hydrogen_ground_state.jl`

## 7. What changes for multi-electron atomic work?

The basic structure stays the same:

1. build a radial basis
2. build a radial quadrature grid
3. check diagnostics
4. build one-body operators
5. build the current two-index radial multipole operators

For the current v0 package, a convenient higher-level entry point is:

```julia
ops = atomic_operators(rb, grid; Z = 2.0, lmax = 2)

ops.overlap
ops.kinetic
ops.nuclear

centrifugal(ops, 2)
multipole(ops, 1)
```

In the current package, `multipole_matrix` and `multipole(ops, L)` mean the supported two-index IDA-style radial multipole matrices. They are not exact four-index electron-electron tensors.

## 8. When should you try x-gaussians?

For a first one-electron calculation such as hydrogen, it is perfectly reasonable to start with:

```julia
xgaussians = XGaussian[]
```

If you are doing more serious radial interaction work and the near-origin diagnostics look unsatisfactory, then trying one or two x-gaussians is reasonable.

A simple one-function experiment is:

```julia
xgaussians = [XGaussian(alpha = 0.2)]
```

A more serious two-function starting point is something like:

```julia
xgaussians = [XGaussian(alpha = 0.1), XGaussian(alpha = 0.025)]
```

These are starting points, not magic numbers. If you are using x-gaussians, check the diagnostics rather than assuming one choice is always best.

## 9. Common beginner mistakes

Here are the most common mistakes new users make.

### Confusing the basis with the quadrature grid

They are different objects with different jobs. Build both.

### Treating a tiny smoke-test basis as a production recommendation

A very small basis may be fine for testing that code runs, but not for a serious atom calculation.

### Worrying about `stencil` too early

You do not need `stencil(f)` to start doing calculations. It is there when you want to inspect the exact Gaussian expansion behind a function.

### Choosing a cutoff that is too small

If the tail of the state matters and the basis `rmax` or an explicit `quadrature_rmax` is too small, the calculation will look worse for a physical reason, not because the basis idea failed.

## 10. Where to go next

After this page, the most useful next reads are:

* [`docs/recommended_atomic_setup.md`](recommended_atomic_setup.md)
* [`docs/terminology.md`](terminology.md)

And after that, the runnable examples in `examples/` are the best way to get comfortable with the package.
