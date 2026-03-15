# Quickstart: a first radial workflow

This page answers one practical question:

**If I want to do a first atom-centered radial calculation with GaussletBases, what do I actually do?**

The short answer is:

1. choose a mapping and a radial basis recipe
2. build the basis
3. run a quick diagnostic check
4. build a quadrature grid
5. form one-body operators
6. test the setup on hydrogen

That is the basic workflow.

If you are new to the package, this is the best page to read after the README.

## 1. Choose a radial recipe

For an atom of nuclear charge `Z`, the package documentation recommends starting with:

- `s = 0.2`
- `c = s / (2Z)`
- `tail_spacing = 10.0` through the default `AsinhMapping`
- `tails = 6`
- `odd_even_kmax = 6`
- `reference_spacing = 1.0`
- `rmax = 30.0` bohr as a good first-row starting point

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
```

This is the package’s standard starting recipe. It is not the only workable choice, but it is a good first choice if you want one recommendation rather than a large tuning discussion.

## 2. Build the basis

```julia
rb = build_basis(spec)
```

That gives you a compact radial basis for the reduced radial function `u(r) = r R(r)`.

At this point you already have the actual basis functions. For example:

```julia
f = rb[4]

f(0.3)
reference_center(f)
center(f)
```

If you are just learning the package, it is enough to know that:
- `reference_center(f)` is the center before the coordinate map is applied
- `center(f)` is the physical-space center after mapping

## 3. First check: diagnostics

The quickest first sanity check is:

```julia
diag = basis_diagnostics(rb)

diag.overlap_error
diag.D
```

This one-argument form is the simplest path for a new user. The package chooses its own conservative quadrature internally for the diagnostic.

The two most useful quantities at first are:

- `diag.overlap_error`: how close the numerically checked overlap is to exact orthonormality
- `diag.D`: an aggregate measure of how much the nominal basis centers differ from first-moment centers on the diagnostic grid

You do not need to obsess over `D` on the first day. The point is simply that it helps you notice when the near-origin behavior may need more attention.

## 4. Build an explicit quadrature grid

When you want to build operators or inspect moment centers directly, make the quadrature grid explicitly:

```julia
grid = radial_quadrature(rb)
```

That is the default recommended call for most users.

The package chooses a conservative quadrature profile automatically. If you later want more explicit control, the advanced form is:

```julia
grid = radial_quadrature(rb; accuracy = :high, quadrature_rmax = 30.0)
```

You can think of the quadrature grid as the numerical integration partner of the basis. The basis and the quadrature grid are intentionally separate objects.

## 5. Compare different center notions if you want intuition

Near the origin, several “center” notions need not agree exactly. That is normal, and it is part of why the diagnostics exist.

For one basis function:

```julia
reference_center(f)
center(f)
moment_center(f, grid)
```

The simplest interpretation is:

- `reference_center(f)`: center before mapping
- `center(f)`: mapped physical-space center
- `moment_center(f, grid)`: first-moment center computed numerically on the grid

If these begin to disagree in a way that looks unstable under reasonable changes of basis or quadrature, that is a sign to look at the setup more carefully.

## 6. Build one-body operators

For a radial one-electron problem, the basic matrices are:

```julia
S = overlap_matrix(rb, grid)
T = kinetic_matrix(rb, grid)
V = nuclear_matrix(rb, grid; Z = Z)
C0 = centrifugal_matrix(rb, grid; l = 0)
```

and the one-electron Hamiltonian is then:

```julia
H = T + V + C0
```

At this stage, everything is still one-electron and radial. That makes it the cleanest place to understand what the package is doing.

## 7. Solve hydrogen as the first real test

Hydrogen is the best first scientific example because it is simple, familiar, and free of electron-electron complications.

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

The exact nonrelativistic ground-state energy is `-0.5 Ha`, so this is a very direct way to check whether the basis and quadrature are working together well.

A runnable version is in:

- `examples/04_hydrogen_ground_state.jl`

## 8. What changes for multi-electron radial work?

The basic radial structure stays the same:

1. build a radial basis
2. check diagnostics
3. build a quadrature grid
4. build one-body operators
5. build the current two-index radial multipole operators

For that current v0 workflow, the package also provides:

```julia
ops = atomic_operators(rb, grid; Z = 2.0, lmax = 2)

ops.overlap
ops.kinetic
ops.nuclear

centrifugal(ops, 2)
multipole(ops, 1)
```

In the current package, `multipole_matrix` and `multipole(ops, L)` mean the supported two-index IDA-style radial multipole matrices. They are not exact four-index electron-electron tensors.

## 9. What about x-gaussians?

For a first one-electron example such as hydrogen, it is perfectly reasonable to start with:

```julia
xgaussians = XGaussian[]
```

If you later care more about near-origin behavior in the current radial interaction approximation, it can be useful to try one or two x-gaussians.

A simple one-function experiment is:

```julia
xgaussians = [XGaussian(alpha = 0.2)]
```

A more serious two-function starting point is something like:

```julia
xgaussians = [XGaussian(alpha = 0.1), XGaussian(alpha = 0.025)]
```

Those are only starting points. If you use x-gaussians, check the diagnostics rather than treating any one choice as sacred.

## 10. Common beginner mistakes

The most common beginner mistakes are:

### Confusing the basis with the quadrature grid

They are different objects with different jobs.

### Treating a tiny demo basis as a production recommendation

A tiny basis is useful for smoke tests. It is not automatically a good atom setup.

### Worrying about `stencil` too early

You do not need `stencil(f)` to start doing calculations. It is there when you want to inspect the exact Gaussian expansion behind a function.

### Choosing `rmax` too small

If the state has meaningful tail weight and the represented radius is too short, the result will look worse for a real physical reason.

## 11. Where to go next

After this page, the most useful next reads are:

- [`recommended_atomic_setup.md`](recommended_atomic_setup.md)
- [`example_guide.md`](example_guide.md)
- [`terminology.md`](terminology.md)

If you want the advanced contraction/hierarchy side of the repo after that, then read:

- [`intermediate_primitive_layer.md`](intermediate_primitive_layer.md)
