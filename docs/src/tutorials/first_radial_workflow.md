# First radial workflow

This is the best first tutorial for a new user.

The first radial path is:

1. choose a mapping and a radial basis recipe
2. build the basis
3. run a quick diagnostic check
4. build a quadrature grid
5. form one-body operators
6. test the setup on hydrogen

## A first working recipe

For a first atom-centered radial calculation, start with:

- `s = 0.2`
- `c = s / (2Z)`
- `rmax = 30.0` bohr

In code:

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
```

`AsinhMapping(c = c, s = s)` uses ordinary Julia keyword arguments. At this
level, the only two parameters you need to recognize are:

- `s`: roughly controls the overall radial spacing
- `c`: roughly controls how much extra resolution is concentrated near the nucleus

The fuller tuning discussion lives in
[Recommended atomic setup](../howto/recommended_atomic_setup.md).

## Diagnostics and quadrature

After building the basis, the first quick numerical check is:

```julia
diag = basis_diagnostics(rb)
diag.overlap_error
diag.D
```

Then build the explicit quadrature grid used for operator construction:

```julia
grid = radial_quadrature(rb)
```

One of the central ideas in the package is that the basis and the quadrature
grid are separate objects.

## First useful physical result: hydrogen

The cleanest first scientific check is the radial hydrogen ground state:

```julia
H = kinetic_matrix(rb, grid) +
    nuclear_matrix(rb, grid; Z = Z) +
    centrifugal_matrix(rb, grid; l = 0)

E0 = minimum(real(eigen(Hermitian(H)).values))
println("Lowest hydrogen energy: ", E0)
```

The exact nonrelativistic ground-state energy is `-0.5 Ha`, so this first
result checks the basis and the quadrature together.

The corresponding runnable example is `examples/04_hydrogen_ground_state.jl`.

## What you usually do next

Once hydrogen is working, the normal next steps are:

```julia
ops = atomic_operators(rb, grid; Z = Z, lmax = 2)
```

That moves you from a one-electron radial test into the current small atomic
operator line.

For more detail after this tutorial, continue with:

- [Recommended atomic setup](../howto/recommended_atomic_setup.md)
- [Example guide](../howto/example_guide.md)
- [Current atomic branch](../explanations/current_atomic_branch.md)
