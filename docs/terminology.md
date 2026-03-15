# Terminology guide

GaussletBases uses a small amount of package-specific language. This page gives short plain-language meanings for the terms that matter most.

The goal is not to force a special vocabulary on you. The goal is simply to make the package easier to read.

## Gausslet

A gausslet is a localized basis function built from a sum of Gaussians.

For most users, it is enough to think:

**a gausslet is a smooth localized basis function with an explicit Gaussian construction underneath it**

## Mapping

A mapping controls how resolution is distributed in physical space.

In radial work, this is what lets you place more resolution near the nucleus without forcing the whole outer region to use that same fine spacing.

A common example is:

```julia
map = AsinhMapping(c = s / (2Z), s = s)
```

## Radial basis

A radial basis is the set of basis functions used for the reduced radial function

```math
u(r) = r R(r).
```

In the package, the basis object is `RadialBasis`, and one function from it can be accessed as `rb[i]`.

## Quadrature grid

A quadrature grid is the numerical integration grid used to evaluate integrals.

One of the central ideas in GaussletBases is:

**the basis is not the quadrature grid**

You build a basis with `build_basis(...)`, and you build a quadrature grid separately with `radial_quadrature(...)`.

## Primitive

A primitive is one of the lower-level Gaussian-type building blocks from which a function or basis is assembled.

Examples in this package include:

- `Gaussian`
- `HalfLineGaussian`
- `XGaussian`
- mapped versions such as `Distorted(Gaussian(...), mapping)`

## Primitive set

A primitive set is the shared ordered collection of primitives behind a basis.

For a basis `b`, the preferred high-level call is:

```julia
P = primitive_set(b)
```

If you want the raw primitive list, you can also inspect:

```julia
primitives(b)
```

## Stencil

`stencil(f)` returns the exact Gaussian expansion used to build the function `f`.

This is not standard chemistry language, so if it sounds unfamiliar, that is normal.

In plain words, it means:

**the exact list of primitive terms and coefficients that define the function**

For a full basis, the same information appears as a coefficient matrix:

```julia
C = stencil_matrix(b)
```

## Contraction

Contraction means building a smaller or more useful basis by taking linear combinations of a larger primitive layer.

In the package, this idea appears through the contraction matrix returned by:

```julia
C = stencil_matrix(b)
```

and through helpers such as:

```julia
contract_primitive_matrix(b, A)
```

## `reference_center`, `center`, and `moment_center`

These names are similar, so it helps to separate them.

### `reference_center(f)`

The center before coordinate mapping.

### `center(f)`

The physical-space center after mapping.

### `moment_center(f, grid)`

The first-moment center computed numerically on a supplied quadrature grid.

For radial functions near the origin, these do not have to agree exactly.

## `D`

`D` is the package’s aggregate center-mismatch diagnostic.

Very roughly, it measures how much the nominal basis centers differ from first-moment centers on the quadrature grid.

You do not need to memorize its formula in order to use it well. In practice, it is a warning light for near-origin behavior.

## `multipole_matrix`

In the current v0 package, `multipole_matrix` means the supported two-index IDA-style radial multipole matrix.

It does **not** mean an exact four-index electron-electron tensor.

That larger exact object is a separate future extension.

## `atomic_operators`

`atomic_operators(rb, grid; Z, lmax)` is a convenience constructor for the basic radial operator bundle.

It gives you:

- `ops.overlap`
- `ops.kinetic`
- `ops.nuclear`
- `centrifugal(ops, l)`
- `multipole(ops, L)`

## Final advice

If you are new to the package, you do not need to master all of these terms at once.

A good order is:

1. learn what the radial basis is
2. learn why the quadrature grid is separate
3. learn the recommended atomic setup
4. learn hydrogen
5. only then worry about primitives, stencils, and contraction if you need them
